#download phenodata and rawdata from ArrayExpress

biocLite("ArrayExpress")
library(ArrayExpress)

#location for data
dir.create("/home/MOKA/veera/GEO/ArrayExpress")

#Array express identificators as vector
PSO_AD_mastertable_only_ArrayExpress = read.csv(paste("/home/MOKA/veera/GEO/PSO_AD_mastertable_only_ArrayExpress.txt", sep = ""))
arrayid = as.vector(PSO_AD_mastertable_only_ArrayExpress[[1]])


for(i in 1:length(arrayid)) {
  tryCatch(
    expr = {
      array_folder = paste("/home/MOKA/veera/GEO/ArrayExpress/", arrayid[i], sep = "")
      dir.create(array_folder)
      arraydata = ArrayExpress(arrayid[i], path = array_folder, save = TRUE, dataCols = NULL, drop = TRUE)
      
      #getAE(arrayid, path = array_folder, type = "full")
      
      arrayfile = paste(array_folder, "/", arrayid[i], ".sdrf.txt", sep = "")
      arraydf = read.csv(arrayfile, header = TRUE, sep = "", dec = ".")
      
      #create a new folder for phenodata (arraydf)
      r_data_folder = (paste(array_folder, "/", arrayid[i], "_R_data", sep = ""))
      dir.create(r_data_folder)
      save(arraydf, file = paste(r_data_folder, "/", arrayid[i], ".RData", sep = ""))
      write.table(x=arraydf, file = paste(r_data_folder, "/", arrayid[i], ".txt", sep = ""))
    },
    error = function(cond) {
      message(paste("Series did not work", arrayid[i]))
      message(cond)
      return(NA)
    },
    warning = function(cond) {
      message(paste("Series caused a warning", arrayid[i]))
      message(cond)
    },
    finally = {
      message(paste("Processed series:", arrayid[i]))
    }
  )
}

#check downloaded phenodata
View(arraydf)

all_array_id = list.files("/home/MOKA/veera/GEO/ArrayExpress/")

arrayid =all_array_id[7]
GSE= arrayid
GSM = NA
GSM = gsub(".f.*", "", GSM)
#arrayid = "E-MTAB-729"
file_name = as.vector(arraydf$REF.4)[1]
file_name = c(as.vector(arraydf$REF.4)[1:4], as.vector(arraydf$REF.5)[5:7])
grep("gz", file_name)
grep("gz", file_name2)
file_name = c(file_name[1:138], file_name2[139:204], file_name[205:210], file_name2[211:270], file_name[271:276], file_name2[277:336], file_name[337:360])
file_name = NA
diagnosis = as.vector(arraydf$Protocol.3)
diagnosis[5:7]="atopic_dermatitis"
diagnosis = gsub("pso.*","psoriatic_arthritis", diagnosis)
lesional = as.vector(arraydf$REF.3)
lesional = gsub("nor.*", "non_lesional", lesional)
lesional = gsub("L", "l", lesional)
treatment = NA
dose = NA
dose_unit = NA
time_point = NA

time_point_unit= NA
tissue = as.vector(arraydf$part.)
tissue = "skin"
anatomical_site = NA
platform = "GPL6947"
array = NA
slide = NA
dye = NA
slide = NA
gender = as.vector(arraydf$Protocol)
gender = gsub("m.*","m",gender )
gender = NA
age = arraydf$Characteristics.OrganismPart.
ethnicity = NA
patient_id = c(1:4, rep(5:7, each =2))
patient_id = as.vector(arraydf$Characteristics.age.)
patient_id = gsub(".*X", "", patient_id)
patient_id[which(patient_id =="")] =10
length(patient_id)

load(paste0("/home/MOKA/veera/GEO/ArrayExpress/", arrayid, "/", arrayid, "_R_data/", arrayid,".RData"))

#organize columns

arraydf_checked = data.frame(GSE, GSM, file_name, platform, diagnosis, lesional, treatment, dose, dose_unit, time_point, time_point_unit,
                               tissue, anatomical_site, array, dye, slide, 
                               gender, age, ethnicity, patient_id, row.names = 1:10, stringsAsFactors = FALSE)
# names(arraydf_checked) = c("GSE", "GSM", "file_name","platform", "diagnosis", "lesional", "treatment", "dose", "dose_unit",
#                              "time_point","time_point_unit", "tissue", "anatomical_site", "array", "dye", "slide",
#                              "gender", "age", "ethnicity")
arraydf_checked$arrayid = gsub("-","_", arraydf_checked$arrayid)

save(arraydf_checked, file = paste0("/home/MOKA/veera/GEO/ArrayExpress/", arrayid, "/", arrayid,".RData"))
write.table(x=arraydf_checked, file = paste0("/home/MOKA/veera/GEO/ArrayExpress/", arrayid, "/", arrayid,".txt"), sep = "\t", quote = FALSE)

