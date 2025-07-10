# MegaLMM HPC parallel script
# each run for each trait at each time point can be run in parallel to improve
# computational efficiency.

library(MegaLMM)
library(tibble)
library(tidyr)
library(dplyr)

args = commandArgs(trailingOnly = T)

# all_BLUPs_file
all_BLUPs_f = args[1]

# all_BLUPs_HS_file
all_BLUPs_HS_f = args[2]

# Kinship file
Kin_f = args[3]

# genotypes
geno_f = args[4]

# CV matrix
CV_mat_file = args[5]

# output folder
out_path = args[6]

# trait
trait = args[7]

# DAT
dat = as.numeric(args[8])

# runID
run = as.numeric(args[9])




## parallel access to files is staggered to prevent reading errors 
read_csv_retry <- function(file,rn=F, max_retries = 3, delay = 5) {
  attempts <- 0
  while (attempts < max_retries) {
    attempts <- attempts + 1
    tryCatch({
      if(rn){
        data <- read.csv(file,row.names = 1)
      }else{
        data <- read.csv(file)
      }
      return(data)
    }, error = function(e) {
      if (attempts == max_retries) stop("Max retries reached. Unable to read file.")
      Sys.sleep(delay) # Wait before retrying
    })
  }
}

read_table_retry <- function(file, max_retries = 3, delay = 5) {
  attempts <- 0
  while (attempts < max_retries) {
    attempts <- attempts + 1
    tryCatch({
      data <- read.table(file)
      return(data)
    }, error = function(e) {
      if (attempts == max_retries) stop("Max retries reached. Unable to read file.")
      Sys.sleep(delay) # Wait before retrying
    })
  }
}

if(!dir.exists(out_path)){
  dir.create(out_path)
}

if(!dir.exists(file.path(out_path,"Accuracy"))){
  dir.create(file.path(out_path,"Accuracy"))
  dir.create(file.path(out_path,"Predictions"))
}

all_BLUPs = read_csv_retry(all_BLUPs_f, rn = T)
all_BLUPs_HS = read_csv_retry(all_BLUPs_HS_f, rn = T)
K = read_csv_retry(Kin_f, rn = T)
CV_mat = read_csv_retry(CV_mat_file,rn=T)

names(all_BLUPs)[c(1,4)]=c("X","BLUP")

genotypes = sort(read_table_retry(geno_f)$X)

K_cut = K[genotypes,genotypes]


acc_frame = data.frame()
result_frame = data.frame()

if(trait == "RGB1_Plant_Avg_HEIGHT_MM"){
  foc_dat = dat+1
}else{
  foc_dat = dat
}

BLUPS_foc_spec <- all_BLUPs %>%
  filter(DAT==foc_dat) %>%
  filter(Trait==trait) %>%
  filter(X %in% genotypes) %>%
  select(X,BLUP)

BLUPS_HS_spec <- all_BLUPs_HS %>%
  filter(DAT==dat) %>%
  filter(X %in% genotypes)

HS_mat = pivot_wider(BLUPS_HS_spec,id_cols = c(X),names_from = Trait,values_from = BLUP) %>%
  select(where(~any(. !=1))) %>%
  column_to_rownames(var="X")

HS_mat = HS_mat[,-which(apply(HS_mat,2,var)==0)]

HS_mat_n_col = ncol(HS_mat)


Y = cbind(BLUPS_foc_spec$BLUP[match(BLUPS_foc_spec$X,genotypes)],HS_mat[genotypes,])
names(Y)[1]=trait


# Set up 5 fold CV 

k_fold = length(unique(CV_mat[,1]))


for(fold in 1:k_fold){
  fold_ID = fold
  
  test_ids = CV_mat[,run]==fold_ID
  Y_train = Y_testing = Y
  Y_train[test_ids,1] = NA
  Y_testing[!test_ids,1] = NA
  
  ## MegaLMM
  run_parameters = MegaLMM_control(
    h2_divisions = 20, 
    # Each variance component is allowed to explain between 0% and 100% of the
    # total variation. How many segments should the range [0,100) be divided 
    # into for each random effect?
    burn = 1000,  
    # number of burn in samples before saving posterior samples. I set this to 
    # zero and instead run the chain in small chunks, doing the burning manually, a
    # s described below.
    thin = 2, #,
    # during sampling, we'll save every 2nd sample to the posterior database.
    K = 100 # number of factors. With 19 traits, this is likely way higher than needed.
  )
  
  MegaLMM_state = setup_model_MegaLMM(
    Y = Y_train,  
    # The n x p trait matrix
    # formula = as.formula(paste0("~ ",fixed_effect_t," (1|Genotype)")),  
    formula = ~ (1|Genotype),  
    # This is syntax like lme4 for mixed effect models. 
    # We specify a fixed effect of population and a random effect for genotype (Line)
    data = data.frame(Genotype=rownames(Y)),         
    # the data.frame with information for constructing the model matrices
    relmat = list(Genotype = K_cut),
    # extra_regressions = list(X=X_fe),
    # A list of covariance matrices to link to the random effects in formula.
    # each grouping variable in formula can be linked to a covariance matrix.
    # If so, every level of the grouping variable must be in the rownames of K.
    # additional rows of K not present in data will still be predicted 
    # (and therefore will use memory and computational time!)
    run_parameters=run_parameters,
    # This list of control parameters created above
    run_ID = paste0(out_path,"/Runs/",trait,"_dat",dat,"_run",run,"_",sprintf('MegaLMM_fold_%02d',fold_ID))
    # A run identifier. The function will create a folder with this name 
    # and store lots of useful data inside it
  )
  
  Lambda_prior = list(
    sampler = sample_Lambda_prec_horseshoe, 
    # function that implements the horseshoe-based Lambda prior 
    # described in Runcie et al 2020. 
    #See code to see requirements for this function.
    # other options are:
    # ?sample_Lambda_prec_ARD,
    # ?sample_Lambda_prec_BayesC
    prop_0 = 0.1,    
    # prior guess at the number of non-zero loadings in the first and most important factor
    delta = list(shape = 3, scale = 1),    
    # parameters of the gamma distribution giving the expected change 
    # in proportion of non-zero loadings in each consecutive factor
    delta_iterations_factor = 100   
    # parameter that affects mixing of the MCMC sampler. This value is generally fine.
  )
  
  priors = MegaLMM_priors(
    tot_Y_var = list(V = 0.5,   nu = 5),      
    # Prior variance of trait residuals after accounting for fixed effects and factors
    # See MCMCglmm for meaning of V and nu
    tot_F_var = list(V = 18/20, nu = 20),     
    # Prior variance of factor traits. This is included to improve MCMC mixing, 
    # but can be turned off by setting nu very large
    h2_priors_resids_fun = function(h2s,n)  1,  
    # Function that returns the prior density for any value of the h2s vector 
    # (ie the vector of random effect proportional variances across all random effects. 
    # 1 means constant prior. 
    # n is the number of h2 divisions above (here=20)
    # 1-n*sum(h2s)/n linearly interpolates between 1 and 0, 
    # giving more weight to lower values
    h2_priors_factors_fun = function(h2s,n) 1, 
    # See above. 
    # sum(h2s) linearly interpolates between 0 and 1,
    # giving more weight to higher values
    # Another choice is one that gives 50% weight to h2==0: ifelse(h2s == 0,n,n/(n-1))
    Lambda_prior = Lambda_prior
    # from above
  )
  
  MegaLMM_state = set_priors_MegaLMM(MegaLMM_state,priors)
  
  
  MegaLMM_state = initialize_variables_MegaLMM(MegaLMM_state)
  
  estimate_memory_initialization_MegaLMM(MegaLMM_state)
  
  MegaLMM_state = initialize_MegaLMM(MegaLMM_state,verbose = T)
  
  MegaLMM_state$Posterior$posteriorSample_params = c('Lambda','F_h2','resid_h2','tot_Eta_prec')
  MegaLMM_state$Posterior$posteriorMean_params = 'Eta_mean'
  
  MegaLMM_state$Posterior$posteriorFunctions = list(
    #B2 = 'B2_F %*% Lambda + B2_R',
    U = 'U_F %*% Lambda + U_R' # total genetic value (factors + residuals)
    #G = 't(Lambda) %*% diag(F_h2[1,]) %*% Lambda + diag(resid_h2[1,]/tot_Eta_prec[1,])', # genetic covariance among traits
    #R = 't(Lambda) %*% diag(1-F_h2[1,]) %*% Lambda + diag((1-resid_h2[1,])/tot_Eta_prec[1,])'#, # Residual covariance among traits
    #h2 = '(colSums(F_h2[1,]*Lambda^2)+resid_h2[1,]/tot_Eta_prec[1,])/(colSums(Lambda^2)+1/tot_Eta_prec[1,])' # additive heritability of each trait
  )
  
  MegaLMM_state = clear_Posterior(MegaLMM_state) 
  
  n_iter = 1000
  for(k in 1:2) {
    print(sprintf('Sampling run %d',k))
    MegaLMM_state = sample_MegaLMM(MegaLMM_state,n_iter) 
    MegaLMM_state = save_posterior_chunk(MegaLMM_state)
    print(MegaLMM_state)
  }
  
  
  # Load individual samples saved on disk
  U_samples = load_posterior_param(MegaLMM_state,'U')
  U_hat = get_posterior_mean(U_samples)
  Eta_mean = load_posterior_param(MegaLMM_state,'Eta_mean')
  
  MegaLMM_Uhat_accuracy = cor(Y_testing[,1],U_hat[,1],use='p')
  MegaLMM_Eta_mean_accuracy = cor(Y_testing[,1],Eta_mean[,1],use='p')
  
  class = rep("Train",nrow(BLUPS_foc_spec))
  class[is.na(fold_ID_matrix[,1])]="Test"
  cur_frame = data.frame(Geno = genotypes,
                         BLUPs = BLUPS_foc_spec$BLUP,
                         U_hat = U_hat[,1],
                         Eta_mean = Eta_mean[,1],
                         Class = class,
                         Run = run,
                         Fold = fold,
                         DAT = dat,
                         Trait = trait,
                         n_WL = HS_mat_n_col)
  t_acc_frame = data.frame(U_hat_accuracy = MegaLMM_Uhat_accuracy,
                           Eta_mean_accuracy = MegaLMM_Eta_mean_accuracy, 
                           Run = run, 
                           Fold = fold,
                           DAT = dat,
                           Trait = trait)
  if(fold==1){
    result_frame = cur_frame
    acc_frame = t_acc_frame
  }else{
    result_frame = rbind(result_frame,cur_frame)
    acc_frame = rbind(acc_frame,t_acc_frame)
  }
}

write.csv(acc_frame,paste0(out_path,"/Accuracy/",trait,"_dat",dat,"_run",run,"_accuracy_frame_MegaLMM.csv"))
write.csv(result_frame,paste0(out_path,"/Predictions/",trait,"_dat",dat,"_run",run,"_prediction_frame_MegaLMM.csv"))