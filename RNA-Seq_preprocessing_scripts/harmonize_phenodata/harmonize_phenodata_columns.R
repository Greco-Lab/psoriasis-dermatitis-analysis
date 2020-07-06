#get similar columns for all phenodata files

# create folder for homogenized data
dir.create("/home/MOKA/veera/GEO/pheno_modified")
pheno_mod_folder = "/home/MOKA/veera/GEO/pheno_modified"

# list all of the unmodified phenodata identifications
gseid = list.files("/home/MOKA/veera/GEO/Pheno_preprosessed2/")


for (i in 1:length(gseid)){
  tryCatch(
    expr = {
      load(paste("/home/MOKA/veera/GEO/Pheno_preprosessed2/", gseid[i], "/", gseid[i], ".RData", sep = ""))
      geopheno = data.frame(lapply(geopheno, as.character), stringsAsFactors = FALSE)
      View(geopheno)
      
      GSE = gseid[i]
      
      GSM = geopheno[grep("GSM", geopheno)]
      if (ncol(GSM) > 1){
        GSM = GSM[1]
      }
      if (ncol(GSM) == 0){
        GSM = NA
      }
    
      lesional = geopheno[grep("lesional|Lesional|LS|NL", geopheno)]
      if (ncol(lesional) > 1){
        lesional = lesional[1]
      }
      if (ncol(lesional) == 0){
        lesional = NA
      }
      
      diagnosis = geopheno[grep("Pso|pso|AD|atopic|Atopic", geopheno)]
      if (ncol(diagnosis) > 1){
        diagnosis = diagnosis[1]
      }
      if (ncol(diagnosis) == 0){
        diagnosis = NA
      }
      
      treatment = geopheno[grep("treatment|Treatment", colnames(geopheno))] #NA #geopheno[grep("control", geopheno)[1]]
      if (ncol(treatment) > 1){
        treatment = treatment[1]
      }
      if (ncol(treatment) == 0){
        treatment = NA
      }
      
      dose = geopheno[grep("dose|Dose", colnames(geopheno))]
      if (ncol(dose) > 1){
        dose = dose[1]
      }
      if (ncol(dose) == 0){
        dose = NA
      }
      
      dose_unit = geopheno[grep("mg", geopheno)]
      if (ncol(dose_unit) > 1){
        dose_unit = dose_unit[1]
      }
      if (ncol(dose_unit) == 0){
        dose_unit = NA
      }
      
      time_point = geopheno[grep("time|Time", colnames(geopheno))]
      if (ncol(time_point) == 0){
        time_point = geopheno[grep("week|Week|month|Month|day|Day", geopheno)]
      }
      if (ncol(time_point) > 1){
        time_point = time_point[1]
      }
      if (ncol(time_point) == 0){
        time_point = NA
      }
      
      platform = geopheno[grep("GPL", geopheno)]
      if (ncol(platform) > 1){
        platform = platform[1]
      }
      if (ncol(platform) == 0){
        platform = NA
      }
      
      array = geopheno[grep("array", colnames(geopheno))]
      if (ncol(array) > 1){
        array = array[1]
      }
      if (ncol(array) == 0){
        array = NA
      }
      
      dye = geopheno[grep("cy|Cy|biotin|Biotin", geopheno)]
      if (ncol(dye) == 0){
        dye = geopheno[grep("dye", colnames(geopheno))]
      }
      if (ncol(dye) > 1){
        dye = dye[1]
      }
      if (ncol(dye) == 0){
        dye = NA
      }
      
      slide = geopheno[grep("slide", colnames(geopheno))]
      if (ncol(slide) > 1){
        slide = slide[1]
      }
      if (ncol(slide) == 0){
        slide = NA
      }
      
      file_name = geopheno[grep("files", colnames(geopheno))]
      if (ncol(file_name) > 1){
        file_name = file_name[1]
      }
      if (ncol(file_name) == 0){
        file_name = NA
      }
      
      tissue = geopheno[grep("skin|Skin|blood|Blood", geopheno)]
      if (ncol(tissue) > 1){
        tissue = tissue[1]
      }
      if (ncol(tissue) == 0){
        tissue = NA
      }
      
      anatomical_site = geopheno[grep("back|Back|trunk|Trunk|scalp|Scalp|buttock|Buttock|Arm|arm|Leg|leg", geopheno)]
      if (ncol(anatomical_site) > 1){
        anatomical_site = anatomical_site[1]
      }
      if (ncol(anatomical_site) == 0){
        anatomical_site = NA
      }
      
      gender = geopheno[grep("gender|Gender", colnames(geopheno))]
      if (ncol(gender) > 1){
        gender = gender[1]
      }
      if (ncol(gender) == 0){
        gender = NA
      }
      
      age = geopheno[grep("age", colnames(geopheno))]
      if (ncol(age) > 1){
        age = age[1]
      }
      if (ncol(age) == 0){
        age = NA
      }
      
      ethnicity = geopheno[grep("Ethnicity|ethnicity", colnames(geopheno))]
      if (ncol(ethnicity) == 0){
        ethnicity = geopheno[grep("caucasian|Caucasian|White|white|black|Black|asian|Asian", geopheno)]
      }
      if (ncol(ethnicity) > 1){
        ethnicity = ethnicity[1]
      }
      if (ncol(ethnicity) == 0){
        ethnicity = NA
      }
      
      SCORAD = geopheno[grep("SCORAD|Scorad|corad", colnames(geopheno))]
      if (ncol(ethnicity) > 1){
        ethnicity = ethnicity[1]
      }
      if (ncol(ethnicity) == 0){
        ethnicity = NA
      }
      
      geopheno_modified = data.frame(GSE,GSM, file_name, diagnosis, lesional, treatment, dose, dose_unit, time_point,
                                     tissue, anatomical_site, platform, array, dye, slide, 
                                     gender, age, ethnicity, SCORAD, row.names = 1:nrow(geopheno), stringsAsFactors = FALSE)
      names(geopheno_modified) = c("GSE", "GSM", "file_name", "diagnosis", "lesional", "treatment", "dose", "dose_unit",
                                   "time_point", "tissue", "anatomical_site", "platform", "array", "dye", "slide",
                                   "gender", "age", "ethnicity", "SCORAD")
      
      View(geopheno_modified)
      
      dir.create(paste(pheno_mod_folder, "/", gseid[i], sep = ""))
      mod_folder = paste(pheno_mod_folder, "/", gseid[i], sep = "")
      save(geopheno_modified, file = paste(mod_folder, "/", gseid[i], ".RData", sep = ""))
      write.table(x=geopheno_modified, file = paste(mod_folder, "/", gseid[i], ".txt", sep = ""), sep = "\t", quote = FALSE)
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
