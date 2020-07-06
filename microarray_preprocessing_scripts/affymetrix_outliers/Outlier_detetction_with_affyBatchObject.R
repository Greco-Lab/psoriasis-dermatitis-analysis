#outliers detection

convert_char_cols_to_factor <- function(phTable){
  coltypes <- unlist(lapply(phTable, class))
  print("coltypes")
  print(coltypes)
  print(table(coltypes))
  coltypes.charOnly.idx <- which(coltypes=="character")
  coltypes.charOnly.len <- length(coltypes.charOnly.idx)
  if(coltypes.charOnly.len>0){
    for(i in 1:coltypes.charOnly.len){
      phTable[, coltypes.charOnly.idx[i]] <- factor(phTable[, coltypes.charOnly.idx[i]])
    }
  }
  return(phTable)
}

View(utopia_affy)
all_gseid=as.vector(utopia_affy$GSE)
gseid = all_gseid[54]
load(paste("/Volumes/veera/GEO/pheno_checked/", gseid, "/", gseid, ".RData", sep=""))
View(geopheno_checked)

for (i in 1:length(gseid)){
  print(paste("Started prosessing series", gseid[i]))
  
  cdfLibName = "hgu95av2hsensgcdf" #switch to corresponding annotation file 
  #Read phenotype
  phTable <- read.delim(file= paste("/Volumes/veera/GEO/pheno_checked/", gseid[i], "/", gseid[i], ".txt", sep = ""))#, header=TRUE, sep="\t") #, stringsAsFactors=F, check.names=F, as.is=TRUE)
 
  #Clean phenotype
  source('~/Documents/R_scripts/CleanPh.R', echo=TRUE)
  phTable <- cleanPh(phTable)
  
  sampleColID <- 1
  fileNameColID <- 2
  
  ##Check sample duplicity and missing samples
  sampleIDs <- phTable[,sampleColID]
  if(any(duplicated(sampleIDs))){
    print(paste0("Sample IDs are not unique!\n\nPlease check the phenotype data and ensure that the selected Sample ID column has unique values!"))
    return(NULL)
  }else if(any(sampleIDs=="") || any(sampleIDs==" ") || any(is.na(sampleIDs))){
    print(paste0("Sample IDs contain BLANK and/or NA values!\n\nPlease check the phenotype data and ensure that the selected Sample ID column has complete information!"))
    return(NULL)
  }
  
  rownames(phTable) <- phTable[,sampleColID]
  pdFactor <- convert_char_cols_to_factor(phTable) #With Factors
  
  #Convert phenotype data frame to phenotype annotation object
  pheno <- new("AnnotatedDataFrame", data=pdFactor)
  
  #Run justRMA() to create expression set from affymetrix cel files
  raw_files_dir <- paste("/Volumes/veera/GEO/Pheno_preprosessed2/", gseid[i], "/", gseid[i], "_RAW/", sep = "")
  #raw_files_dir = "/Volumes/veera/GEO/Pheno_preprosessed2/GSE2737/GSE2737_RAW/"
  
  tmpFile <- file.path(raw_files_dir, phTable[1,fileNameColID])
  celHeader <- affyio::read.celfile.header(tmpFile)
  cdfName <- gsub("-|_", "", celHeader$cdfName)
  
  if(!grepl(cdfName, cdfLibName, ignore.case=T)){
    print("CDF annotation mismatch with the CEL files!")
    print(paste0("Need CDF file corresponding to '", cdfName, "'"))
    print(paste0("User provided CDF '", cdfLibName, "'"))
  }else{
    print("CDF annotation matched correctly!")
    print(paste0("CDF Name: ", cdfName, ", CDF Lib Name: ", cdfLibName))
  }
  
  require(hgu95av2hsensgcdf) #remember to change according to annotation file
  eset <- justRMA(filenames=phTable[,fileNameColID], celfile.path=raw_files_dir, phenoData=pheno, sampleNames=phTable[,sampleColID], normalize=TRUE, background=TRUE, cdfname=cdfLibName)
  
  norm.data <- exprs(eset)
  write.table(norm.data, file = paste("/Volumes/veera/GEO/pheno_checked/", gseid[i], "/", gseid[i], "Expression_Matrix_Normalized.txt", sep = ""), sep = "\t", quote = FALSE,row.names = TRUE)
  #write.table(norm.data, file = "/Volumes/veera/GEO/pheno_checked/GSE5667_GPL97/GSE5667_GPL97_Expression_Matrix_Normalized.txt", sep = "\t", quote = FALSE,row.names = TRUE)
  
  ##QC
  #Get AffyBatch Object
  tmpFile <- file.path(raw_files_dir, phTable[,fileNameColID])
  celHeader <- affyio::read.celfile.header(tmpFile)
  print("str(celHeader)")
  print(str(celHeader))
  cdfName <- gsub("-|_", "", celHeader$cdfName)
  print("cdfName:")
  print(cdfName)
  allInstalledPkgs <- rownames(installed.packages())
  cdfNameBioc <- paste0(tolower(cdfName), "cdf")
  print("Bioconductor cdfName:")
  print(cdfNameBioc)
  affCDF <- cdfLibName
  qcCDF <- affCDF
  warned <- FALSE
  if(!any(grepl(cdfNameBioc, allInstalledPkgs))){
    tryCatch(BiocInstaller::biocLite(cdfNameBioc, suppressUpdates=TRUE))
    allInstalledPkgs <- rownames(installed.packages())
    if(any(grepl(cdfNameBioc, allInstalledPkgs))){
      qcCDF <- cdfNameBioc
      library(package=qcCDF, character.only=TRUE)
    }
  }else{
    qcCDF <- cdfNameBioc
    library(package=qcCDF, character.only=TRUE)
  }
  if(!any(grepl(cdfName, affCDF, ignore.case=T))){
    stop(
      paste0("CDF annotation mismatch with the CEL files!\n\nNeed CDF file corresponding to '", cdfName, "'\n\nUser provided CDF '", affCDF, "'")
    )
  }
  
  print("qcCDF:")
  print(qcCDF)
  
  #Create affyBatch object for QC
  
  source('~/Documents/R_scripts/obtain_affyBatch.R', echo=TRUE)
  affyBatchObject = get_affyBatchObject(celDir = file.path(raw_files_dir), fileNames = as.character(phTable[,fileNameColID]), 
                                        pheno = pheno, cdfname = cdfLibName)
  

  ##Generate QC report
  # source('~/qc_report.R', echo=TRUE)
  # affy_qc_report(affyBatchObject, yaqc_QC=TRUE, image_QC=TRUE, NUSE_QC=TRUE, RLE_QC=TRUE, density_QC=TRUE,
  #                output_dir=paste("/Volumes/veera/GEO/pheno_checked/", gseid[i], "/", sep = ""),
  #                output_filename = paste(gseid[i], "_affy_qc.pdf", sep = ""), output_prefix="affyQC")

  print(paste(gseid[i], "prosessed"))
}

for (i in 1:length(gseid)){
  print(paste("started prosessing series", gseid[i]))
  fileFolder = paste("/Volumes/veera/GEO/Pheno_preprosessed2/", gseid[i], "/", gseid[i], "_RAW/", sep = "") #where to read the cel files
  folder = paste("/Volumes/veera/GEO/pheno_checked/", gseid[i], "/", sep = "")
  dir.create(paste(folder, "RData/", sep=""))
  dir.create(paste(folder, "figures/", sep=""))
  dir.create(paste(folder, "figures/RLE/", sep=""))
  dir.create(paste(folder, "figures/NUSE/", sep=""))
  dir.create(paste(folder, "figures/RNADeg/", sep=""))
  RDATAfileFolder = paste(folder, "RData/", sep="") #folder in which save all the intermediate and final results in RData file format
  RLE_fig_folder =  paste(folder, "figures/RLE/", sep="")
  NUSE_fig_folder = paste(folder, "figures/NUSE/", sep="")
  RNADeg_fig_folder = paste(folder, "figures/RNADeg/", sep="")
  no.files = length(list.files(fileFolder))
  print(no.files)
  
  #install.packages("doMC")
  
  library(affy)
  library(affyPLM)
  library(simpleaffy)
  library(doMC)
  #setwd(fileFolder)
  print("Reading cel files...")
  celFiles = list.files(fileFolder) #List of the names of cel files to be preprocessed
  celFiles = paste(fileFolder,celFiles,sep="")
  nFiles = length(celFiles) #number of cel files to analyse
  # RS =celFiles
  
  # The RLE, NUSE and RNA degradation indices will be evaluated on random partition of the pool of samples
  
  RS = list(celFiles)
  
  
  RLE_stat = c()
  RLE_values = c()
  NUSE_stat = c()
  NUSE_values = c()
  RNADeg = c()
  
  print("Starting the analyses...")
  
  
  PSetList = list()
  affyBList = list()
  for(i in 1:length(RS)){
    cat("Reading data...\n")
    #affyB = ReadAffy(filenames = RS[[i]])
    affyB = affyBatchObject
    
    ED=exprs(affyB)
    cat("AffyPLM...\n")
    
    Pset <- fitPLM(affyB)
    PSetList[[i]] = Pset
    affyBList[[i]] = affyB
    
    #   # Relative Log Expression (RLE) values. Specifically,
    #   #these RLE values are computed for each probeset by comparing the expression value
    #   #on each array against the median expression value for that probeset across all arrays.
    cat("RLE...\n")
    
    png(paste(RLE_fig_folder,"RLE_batch",i,".png",sep=""))
    rleP = RLE(Pset)
    dev.off()
    RLE_stat = cbind(RLE_stat,RLE(Pset,type="stat"))
    # RLE_values = cbind(RLE_values,RLE(Pset,type="value"))
    
    #   #Normalized Unscaled Standard Errors (NUSE) can also be used for assessing quality. In
    #   #this case, the standard error estimates obtained for each gene on each array from fitPLM
    #   #are taken and standardized across arrays so that the median standard error for that
    #   #genes is 1 across all arrays
    cat("NUSE...\n")
    
    png(paste(NUSE_fig_folder,"NUSE_batch",i,".png",sep=""))
    nuseP = NUSE(Pset)
    dev.off()
    NUSE_stat = cbind(NUSE_stat,NUSE(Pset,type="stat"))
    #NUSE_values = cbind(NUSE_values,NUSE(Pset,type="value"))
    
    cat("DEG...\n")
    
    #measure the difference in the signal at the 5' and 3' of a gene
    deg = AffyRNAdeg(affyB)
    degS = summaryAffyRNAdeg(deg)
    png(paste(RNADeg_fig_folder,"RNADeg_batch",i,".png",sep=""))
    plotAffyRNAdeg(deg)
    dev.off()
    
    RNADeg = cbind(RNADeg,degS)
    save(PSetList,affyBList,file=paste(RDATAfileFolder,"affy_pset_random_lists.RData",sep=""))
  }
  
  
  save(RNADeg,RLE_stat,NUSE_stat,file=paste(RDATAfileFolder,"outliers_stats.RData",sep=""))
  
  png(paste(RNADeg_fig_folder,"RNADeg_slope.png",sep=""))
  boxplot(RNADeg["slope",])
  dev.off()
  
  png(paste(RLE_fig_folder,"RLE_median.png",sep=""))
  boxplot(RLE_stat["median",])
  dev.off()
  
  png(paste(NUSE_fig_folder,"NUSE_median.png",sep=""))
  boxplot(NUSE_stat["median",])
  dev.off()
  
  print("Computing outliers...")
  
  RLEout <- names(boxplot.stats(RLE_stat["median",])$out)  # outlier values.
  NUSEout <- names(boxplot.stats(NUSE_stat["median",])$out)  # outlier values.
  deg_bstat = boxplot(RNADeg["slope",])$stats[4] #value of the 4% quantile
  DEGout = colnames(RNADeg)[which(RNADeg["slope",]>deg_bstat)]
  outs = union(DEGout,union(RLEout,NUSEout)) #names of all the outliers samples
  
  #M is a matrix of samples that contains a 1 or a 0 if the sample is considered outliers for one of the 6 measures
  M = matrix(0,nrow=length(outs),ncol=3)
  rownames(M) = outs
  colnames(M) = c("RLE","NUSE","DEG")
  M[outs %in% RLEout,"RLE"] = 1
  M[outs %in% NUSEout,"NUSE"] = 1
  M[outs %in% DEGout,"DEG"] = 1
  
  M = cbind(M,rowSums(M))
  colnames(M)[4] = "SUM"
  
  outs_for_at_least_one_method = outs
  outs = rownames(M)[M[,4]>1]
  save(M,outs,outs_for_at_least_one_method,file=paste(RDATAfileFolder,"outliers.RData",sep=""))
  
  print("End...")
}

