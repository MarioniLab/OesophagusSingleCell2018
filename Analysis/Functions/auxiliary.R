#############################################
#### Script to store auxiliary functions ####
#############################################

# load data by tissue
loadData <- function(file_vec){
  
  # load normalized data
  normdata <- sapply(file_vec, readRDS)
  
  # add patients to column names
  ps <- sapply( strsplit(file_vec, "/"), grep, pattern="Patient", value=TRUE)
  for(i in 1:length(normdata)){
    colnames(normdata[[i]]) <- paste(ps[i], colnames(normdata[[i]]), sep="_")
  }
  
  if(length(file_vec)>1){
    # reduce to genes present in all data sets
    genes <- lapply(normdata, row.names)
    genes <- table(unlist(genes))
    genes <- names(genes)[genes>=length(normdata)]
    
    # reduce to genes expressed in at least one sample
    normdata_sub <- sapply(normdata, function(x, g) as.matrix(x[g, ]), g=genes)
    genes_sub <- genes[rowSums(do.call("cbind", normdata_sub))>0]
    normdata_sub <- sapply(normdata, function(x, g) as.matrix(x[g, ]), g=genes_sub)
    
    return(normdata_sub)
  }else{
    return(normdata)
  }
}

# filter contaminating cells
filterData <- function(normdata_list, cd45=TRUE, vim=TRUE){
  # prep
  select <- sapply(normdata_list, function(x) rep(TRUE, ncol(x)))
  
  # helper
  filterhelper <- function(geneid, normdata, select){
    if(geneid%in%row.names(normdata)){
      select <- select & (!normdata[geneid, ]>0)
    }else{
      print("gene not measured")
    }
    return(select)
  }
  
  # filter immune cells
  if(cd45) {
    select <- sapply(1:length(normdata_list), function(i) 
      filterhelper("ENSG00000081237", normdata_list[[i]], select[[i]]))
  }
  
  if(vim) {
    select <- sapply(1:length(normdata_list), function(i) 
      filterhelper("ENSG00000026025", normdata_list[[i]], select[[i]]))
  }
  normdata_list <- sapply(1:length(normdata_list), function(i) 
    return(normdata_list[[i]][,select[[i]]]))
  
  return(normdata_list)
}

# highly variable genes 
getHVGs <- function(normdata_list){
  
  # remove lowly expressed genes
  genes <- lapply(normdata_list, function(x) row.names(x)[rowMeans(x)>0.1])
  genes <- table(unlist(genes))
  genes <- names(genes)[genes==length(normdata_list)]
  normdata_list <- lapply(normdata_list, function(x, g) x[g,], g=genes)
  
  # find highly variable genes, single data set
  hvgHelper <- function(data){
    # fit variance mean relationship
    logdata <- log2(data + 1)
    varfit <- trendVar(logdata)
    decomp <- decomposeVar(logdata, varfit) 
    
    # plot fit
    #plot(varfit$mean, varfit$var)
    #curve(varfit$trend(x), col="red", lwd=2, add=TRUE)
    
    return( decomp )
  }
  par(mfrow=c(1,length(normdata_list)))
  HVG_list <- lapply(normdata_list, hvgHelper)
  
  # combine if more than one sample
  if(length(HVG_list)>1){
    # unpack arguments for combineVar
    HVG.df <- do.call(combineVar, HVG_list)
  }else{
    HVG.df <- HVG_list[[1]]
  }
  
  # get top 1000 most variable genes
  HVG.df <- HVG.df[order(HVG.df$bio, decreasing=TRUE), ]
  HVG <- rownames(HVG.df)[1:1000]
  return(HVG)
}

# batch correct by tissue
batchCorrect <- function(normdata_list, HVG){
  # batch correction based on top highly variable genes
  normdata_list_hvg <- sapply(normdata_list, function(x, hvg) log2(x[hvg, ] +1), hvg=HVG)
  
  # unpack arguments for mnnCorrect
  funcText <- paste0("mnnCorrect(", 
                     paste0("normdata_list_hvg[[", 1:length(normdata_list_hvg), "]]", collapse=", "), 
                     ", cos.norm.in=TRUE, cos.norm.out=TRUE, sigma=0.1)")
  bcdata <- eval( parse(text=funcText) )
  
  return(bcdata)
}