#download rawdata and phenodata from GEO

library(GEOquery)
library(filesstrings)
install.packages("Rcpp")
install.packages("filesstrings")
installed.packages("strings")

#direction for data
dir.create("/home/MOKA/veera/GEO/Pheno_preprosessed/")

readLocalPheno <- function(path) {
  files = list.files(path)
  gsm = sapply(files, function(x) substr(x,0,regexpr('(_|\\.|-)', x)[1]-1))
  arrays = sapply(files, function(x) substr(x,regexpr('_', x)[1] + 1, regexpr('_[A-Z]',x)[1] - 1))
  slides = sapply(files, function(x) substr(x,regexpr('_[1-4]_', x)[1]+1, regexpr('_[1-4]\\.txt',x)[1]+1))
  data.frame(gsm, arrays, slides, files, row.names = gsm)
}

#vector of GEO indentificators
gseid = as.vector(PSO_AD_mastertable_only_GEO[[1]])

#folder for each dataset
dir.create("/home/MOKA/veera/GEO/Pheno_preprosessed_rerun2/")
rerun_folder = "/home/MOKA/veera/GEO/Pheno_preprosessed_rerun2/"
dir.create(rerun_folder)

gseid = as.vector(rerun_GEO[[1]])


for(i in 1:length(gseid)){
  tryCatch(
    expr = {
      #create a folder
      folder_name = paste(rerun_folder, gseid[i], sep = "")
      print(folder_name)
      dir.create(folder_name)
      
      #collect phenodata from GEO
      geodata=getGEO(GEO = gseid[i], destdir = folder_name, GSEMatrix = TRUE, getGPL = FALSE)
      geopheno = pData(phenoData(geodata[[1]]))
      #geopheno = rbind(pData(phenoData(geodata[[1]])),pData(phenoData(geodata[[2]])))

      ##download RAW data
      raw_data = getGEOSuppFiles(GEO = gseid[i], makeDirectory = TRUE, baseDir = folder_name)

      #create folder for RAW data
      raw_data_folder = paste(folder_name, "/", gseid[i], "_RAW", sep = "")
      dir.create(raw_data_folder)
      
      #open RAW data in created folder
      rawdata = paste(folder_name, "/", gseid[i], "/", gseid[i], "_RAW.tar", sep = "")
      untar(rawdata, exdir = raw_data_folder)
      list.files(raw_data_folder)
      
      ##if there are multiple files with same raw name
      #additional_folder = (paste(folder_name, "/additional_files", sep = ""))
      #dir.create(additional_folder)
      #move_these = list.files(raw_data_folder, pattern = "EXP")
      #file.move(paste(raw_data_folder, "/", move_these, sep = ""), additional_folder, overwrite = TRUE)
      
      #list.files(raw_data_folder)
      localpheno = readLocalPheno(path = raw_data_folder)
      geopheno = base::merge(geopheno, localpheno, by = "row.names")
      
      save(geopheno, file = paste(folder_name, "/", gseid[i], ".RData", sep = ""))
      write.table(x=geopheno, file = paste(folder_name, "/", gseid[i], ".txt", sep = ""), sep = "\t", quote = FALSE)
    },
    error = function(cond) {
      message(paste("Series did not work", gseid[i]))
      message(cond)
      return(NA)
    },
    warning = function(cond) {
      message(paste("Series caused a warning", gseid[i]))
      message(cond)
    },
    finally={
      message(paste("Processed series:", gseid[i]))
    }
  )
}