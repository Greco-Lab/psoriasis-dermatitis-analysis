#Add SVA's to phenodata tables

##arrange df vars by position
##'vars' must be a named vector, e.g. c("var.name"=1)
arrange.vars <- function(data, vars){
  ##stop if not a data.frame (but should work for matrices as well)
  stopifnot(is.data.frame(data))
  
  ##sort out inputs
  data.nms <- names(data)
  var.nr <- length(data.nms)
  var.nms <- names(vars)
  var.pos <- vars
  ##sanity checks
  stopifnot( !any(duplicated(var.nms)), 
             !any(duplicated(var.pos)) )
  stopifnot( is.character(var.nms), 
             is.numeric(var.pos) )
  stopifnot( all(var.nms %in% data.nms) )
  stopifnot( all(var.pos > 0), 
             all(var.pos <= var.nr) )
  
  ##prepare output
  out.vec <- character(var.nr)
  out.vec[var.pos] <- var.nms
  out.vec[-var.pos] <- data.nms[ !(data.nms %in% var.nms) ]
  stopifnot( length(out.vec)==var.nr )
  
  ##re-arrange vars by position
  data <- data[ , out.vec]
  return(data)
}

#list datasets
pheno = list.files("/Users/vh421689/Documents/eUtopia/phenodata/harmonized/")
processed = list.files("/Users/vh421689/Documents/eUtopia/processed/")
setdiff(pheno, processed)

#folder for finished phenodata
dir.create("/Users/vh421689/Documents/eUtopia/phenodata/SVA")


##add SVAs to harmonized phenodata tables

for(i in 1:length(processed)){
  load(paste0("/Users/vh421689/Documents/eUtopia/phenodata/harmonized/", processed[i], "/", processed[i],".RData"))
  print(paste("processing", phenodata$GSE[1]))
  print(dim(phenodata))
  phenodata = unique(phenodata)
  file_name = list.files(paste0("/Users/vh421689/Documents/eUtopia/processed/", processed[i]), pattern = "Phenotype_Table", full.names = T)
  SVA_df = read.delim(file_name)
  SVA_df2 = as.data.frame(SVA_df[,c(1,grep("sva", colnames(SVA_df)))])
  #dir.create(paste0("/Users/vh421689/Documents/eUtopia/phenodata/SVA_rerun/", processed[i]))
  if(dim(SVA_df2)[2]==1){
    print("no SVA's")
    phenodata = arrange.vars(phenodata, c("GSE" = 1))
    save(phenodata, file = paste0("/Users/vh421689/Documents/eUtopia/phenodata/SVA/", processed[i], "/", processed[i], ".RData"))
    write.table(x=phenodata, file = paste0("/Users/vh421689/Documents/eUtopia/phenodata/SVA/", processed[i], "/", processed[i], ".txt"), 
                quote = F, sep = "\t")
  }else if (dim(SVA_df2)[2]>1 & SVA_df2$GSM[1] == phenodata$GSM[1]){
    print(paste(processed[i], "has", (dim(SVA_df2)[2]-1)/2, "SVAd's"))
    phenodata = merge(phenodata, SVA_df2, by = "GSM", all = T)
    phenodata = arrange.vars(phenodata, c("GSE" = 1))
    save(phenodata, file = paste0("/Users/vh421689/Documents/eUtopia/phenodata/SVA/", processed[i], "/", processed[i], ".RData"))
    write.table(x=phenodata, file = paste0("/Users/vh421689/Documents/eUtopia/phenodata/SVA/", processed[i], "/", processed[i], ".txt"), 
                quote = F, sep = "\t")
  } else{
    print("SVA and phenodata not matching")
  }
  print(dim(phenodata))
  #View(phenodata)
  #phenodata = arrange.vars(phenodata, c("GSE" = 1))
  print(paste(processed[i], "processed"))
}
  



