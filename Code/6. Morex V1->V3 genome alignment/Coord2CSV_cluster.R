# genome alignment => SNP translation btw. MorexV1 and MorexV3
args = commandArgs(trailingOnly = T)

read.coords <- function(coords_file){
  varnames = c("S1","E1","S2","E2","LEN1","LEN2","Identity[%]","LenR","LenQ","CovR","CovQ","Tags")
  df = data.frame(matrix(NA,nrow=1,ncol=12))
  names(df)=varnames
  df_idx = 1
  lines = readLines(coords_file)
  sf=F
  for(line in lines){
    if(sf){
    elements = unlist(strsplit(gsub("\\|","",trimws(line)),"\\s+"))
    
    elements[12] = paste(elements[12:length(elements)],collapse = " ")
    df[df_idx,] = elements[1:12]
    df_idx = df_idx+1
    print(df_idx)
    }
    if(startsWith(line,"=")){
      sf=T
    }
  }
  return(df)
}

read.coords2 <- function(coords_file){
  varnames = c("S1","E1","S2","E2","LEN1","LEN2","Identity[%]","LenR","LenQ","CovR","CovQ","Tags")
  
  
  lines <- readLines(coords_file)
  
  lines <- lines[-c(1:5)]
  split_lines <- strsplit(trimws(lines), "\\|")
  
  # Further split each sub-element by whitespace and flatten the structure
  df_list <- lapply(split_lines, function(x) unlist(strsplit(x, "\\s+")))
  
  # Convert list to data frame
  df <- as.data.frame(do.call(rbind, df_list))
  max_col=ncol(df)
  for(col_i in max_col:1){
    if(all(df[,col_i]=="")){
      df = df[,-col_i]
    }
  }
  names(df)[c(1:12)] = varnames
  
  return(df)
}

for(file in list.files(args[1])){
  filename = file.path(args[1],file)
  frame = read.coords2(filename)
  write.csv(frame,paste0(args[2],file,".csv"))
}

