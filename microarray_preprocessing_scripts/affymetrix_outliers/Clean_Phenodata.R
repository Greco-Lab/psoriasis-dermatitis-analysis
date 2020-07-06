## Clean Phenotype table
cleanPh <- function(phTable){
  coltypes <- unlist(lapply(phTable, class))
  print("coltypes")
  print(coltypes)
  print(table(coltypes))
  coltypes.charOnly.idx <- which(coltypes=="character")
  coltypes.nonChar.idx <- which(!coltypes=="character")
  coltypes.charOnly.len <- length(coltypes.charOnly.idx)
  coltypes.nonChar.len <- length(coltypes.nonChar.idx)
  remInfo <- 0
  remStr <- ""
  if(coltypes.charOnly.len>0){
    phTable.charOnly <- phTable[, coltypes.charOnly.idx, drop=F]
    print("dim(phTable.charOnly)")
    print(dim(phTable.charOnly))
    numCheck <- unlist(lapply(phTable.charOnly, function(col){if(suppressWarnings(all(is.na(as.numeric(col))))){1}else{0}}))
    intCheck <- unlist(lapply(phTable.charOnly, function(col){if(suppressWarnings(all(is.na(as.integer(col))))){1}else{0}}))
    doubleCheck <- unlist(lapply(phTable.charOnly, function(col){if(suppressWarnings(all(is.na(as.double(col))))){1}else{0}}))
    allCheck <- numCheck+intCheck+doubleCheck
    if(all(allCheck==0)){
      checkFailed <- names(allCheck[allCheck==0])
      remInfo <- 1
      remStr <- paste0("Following columns are removed because they contain mixed character and numeric data types:\n[",paste0(checkFailed, collapse=", "), "]\n\n")
      if(coltypes.nonChar.len>0){
        phTable.nonChar <- phTable[, coltypes.nonChar.idx, drop=F]
        phTable.comb <- phTable.nonChar
      }else{
        remStr <- paste0(remStr, "No column survived filtering!!! Please define phenotype data columns with singular data type.")
        print(remStr)
        return(NULL)
      }
    }else{
      if(any(allCheck==0)){
        checkFailed <- names(allCheck[allCheck==0])
        phTable.charOnly <- phTable.charOnly[,-which(colnames(phTable.charOnly) %in% checkFailed)]
        remInfo <- 1
        remStr <- paste0("Following columns are removed because they contain mixed character and numeric data types:\n[",paste0(checkFailed, collapse=", "), "]\n\n")
      }
      print("str(phTable) -- before:")
      print(str(phTable))
      #phTable.charOnly <- as.data.frame(apply(phTable.charOnly, 2, function(x){sapply(x, function(y){gsub("[ -]", "_", y)})}), stringsAsFactors=F)
      phTable.charOnly <- as.data.frame(apply(phTable.charOnly, 2, function(x){res<-trimws(x); res<-gsub(" +", " ", res, perl=T); res<-gsub("[ -]", "_", res); return(res)}), stringsAsFactors=F)
      if(coltypes.nonChar.len>0){
        phTable.nonChar <- phTable[, coltypes.nonChar.idx, drop=F]
        print("dim(phTable.nonChar)")
        print(dim(phTable.nonChar))
        phTable.comb <- data.frame(phTable.charOnly, phTable.nonChar, stringsAsFactors=FALSE)
      }else{
        phTable.comb <- phTable.charOnly
      }
    }
    colOrgIdx <- sapply(colnames(phTable.comb), function(x){which(colnames(phTable) %in% x)})
    phTable <- phTable.comb[,names(colOrgIdx[order(colOrgIdx)]), drop=F]
  }
  
  print("str(phTable) -- check:")
  print(str(phTable))
  #Remove columns with single level data
  nrlevels <- apply(phTable, 2, function(x){length(levels(factor(x)))})
  nrlevels.singular <- which(nrlevels==1)
  
  if(length(nrlevels.singular)>0){
    remInfo <- 1
    remStr <- paste0(remStr, "Following columns are removed because they contain only single repeated value:\n[",paste0(names(nrlevels.singular), collapse=", "), "]")
    if(length(nrlevels.singular)==ncol(phTable)){
      remStr <- paste0(remStr, "\n\nNo column survived filtering!!! Please define phenotype data columns with singular data type.")
      print(remStr)
      return(NULL)
    }
    col2rem <- which(colnames(phTable) %in% names(nrlevels.singular))
    phTable <- phTable[,-col2rem, drop=F]
  }
  
  #Inform user with the columns removed from the data frame
  if(remInfo==1){
    print(remStr)
  }
  print("str(phTable) -- after:")
  print(str(phTable))
  return(phTable)
}

