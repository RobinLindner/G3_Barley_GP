# SNP mapping 
snps_MorexV1 = read.table("../Data/Genotype/B1K_IPK/B1K_IPK_filtered_sites.txt",header = T,row.names = 1)

chr1 = read.csv("../Data/Genotype/Assembly/chrom1.coord.csv",row.names=1)

selected_SNP = snps_MorexV1[1,]

selection = chr1 %>%
  filter((S2<selected_SNP$Position & E2>selected_SNP$Position) | (S2>selected_SNP$Position & E2<selected_SNP$Position) ) %>%
  filter(Identity... <=100)



position = selected_SNP$Position
result = data.frame(SNP = NA,old_position=NA,new_position=NA,Seq_Identity=NA)
for(i in 1:nrow(selection)){
  per_base_comp_factor = selection$LEN1[i] / selection$LEN2[i]
  query_range=range(selection$S2[i],selection$E2[i])
  relative_pos = position - query_range[1]
  new_pos = selection$S1 + relative_po * per_base_comp_factor
  
  result
  
}

map_SNP_position <- function(selection,selected_SNP){
  position = selected_SNP$Position
  result = data.frame(SNP = NA,old_position=NA,new_position=NA,Seq_Identity=NA)
  for(i in 1:nrow(selection)){
    per_base_comp_factor = selection$LEN1[i] / selection$LEN2[i]
    query_range=range(selection$S2[i],selection$E2[i])
    relative_pos = position - query_range[1]
    new_pos = selection$S1 + relative_pos * per_base_comp_factor
    
    result[i,] = c(selected_SNP$Name,position,new_pos,selection[i,]$Identity...)
    
  }
  return(result)
}

map_SNP_positions <- function(chr_map,chr_SNPs){
  result = data.frame(SNP = NA,old_position=NA,new_position=NA,Seq_Identity=NA)
  k=1
  for(i in 1:nrow(chr_SNPs)){
    old_pos = chr_SNPs$Position[i]
    selection <- chr_map %>%
      filter((S2<old_pos & E2>old_pos) | (S2>old_pos & E2<old_pos) ) 
    
    if(nrow(selection)==0){next}
    for(j in 1:nrow(selection)){
      per_base_comp_factor = selection$LEN1[j] / selection$LEN2[j]
      query_range=range(selection$S2[j],selection$E2[j])
      relative_pos = old_pos - query_range[1]
      new_pos = selection$S1[j] + relative_pos * per_base_comp_factor
      
      result[k,] = c(chr_SNPs$Name[i],old_pos,round(new_pos),selection$Identity...[j])
      k = k+1
      
    }
      
  }
  result$old_position = as.numeric(result$old_position)
  result$new_position = as.numeric(result$new_position)
  result$Seq_Identity = as.numeric(result$Seq_Identity)
  return(result)
}

chr1_SNPs = snps_MorexV1%>%
  filter(Chromosome == 1)

res=map_SNP_positions(chr1,chr1_SNPs)

mappings = res %>%
  group_by(SNP) %>%
  filter(Seq_Identity == max(Seq_Identity)) %>%
  filter(!duplicated(SNP))
  
write.csv(mappings,"../Data/Genotype/Assembly/chr1_SNP_map")



for(i in 1:7){
  chr_t = read.csv(paste0("../Data/Genotype/Assembly/csv_files/chrom",i,".coords.csv"),row.names=1)
  chr_SNP <- snps_MorexV1 %>%
    filter(Chromosome==i)
  res = map_SNP_positions(chr_t,chr_SNP)
  mappings = res %>%
    group_by(SNP) %>%
    filter(Seq_Identity == max(Seq_Identity)) %>%
    filter(!duplicated(SNP)) %>%
    mutate(Chromosome = i)
  if(i == 1){
    full_remap_df = mappings
  }else{
    full_remap_df = rbind(full_remap_df,mappings)
  }
}

write.csv(full_remap_df,"../Data/Genotype/Assembly/SNP_remap.csv")

