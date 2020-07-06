#Harmonize phenodata and go through the original data to check that nothing was left out

library(dplyr)
library(plyr)

View(geopheno_checked)

# Create folder for modified and checked phenodata tables
dir.create("/home/MOKA/veera/GEO/pheno_checked")
checked_folder = "/home/MOKA/veera/GEO/pheno_checked"

#list data identificators
modified_geopheno = list.files("/home/MOKA/veera/GEO/pheno_modified/")
gseid = modified_geopheno

#original geopheno tables
load(paste("/home/MOKA/veera/GEO/Pheno_preprosessed2/", gseid, "/", gseid, ".RData", sep = ""))
geopheno = data.frame(lapply(geopheno, as.character), stringsAsFactors = FALSE)

# geophenotables with homogenized columns
load(paste("/home/MOKA/veera/GEO/pheno_checked/", gseid, "/", gseid, ".RData", sep = ""))
geopheno_checked = geopheno_modified


#because of the differencies in the original phenodata tables homogenization was done column by column


geopheno_checked$diagnosis
geopheno_checked$diagnosis = geopheno_modified$diagnosis
geopheno_checked$diagnosis = gsub(".*AD.*", "atopic_dermatitis", geopheno_checked$diagnosis)
geopheno_checked$diagnosis = gsub(".*Pso.*", "psoriasis", geopheno_checked$diagnosis)
geopheno_checked$diagnosis = gsub(".*Heal.*", "ctrl", geopheno_checked$diagnosis)
geopheno_checked$diagnosis = gsub(".*Art.*", "psoriatic_arthritis", geopheno_checked$diagnosis)
geopheno_checked$diagnosis = "atopic_dermatitis"
geopheno_checked$diagnosis = "psoriasis"


#healthy = ctrl
#blood = non_lesional
geopheno_checked$lesional
geopheno_checked$lesional = geopheno$source_name_ch1
geopheno_checked$lesional = "lesional"
geopheno_checked$lesional = "non_lesional"
geopheno_checked$lesional = gsub(".*LS.*", "lesional", geopheno_checked$lesional)
geopheno_checked$lesional = gsub(".*NonAD.*", "ctrl", geopheno_checked$lesional)
geopheno_checked$lesional = gsub(".*NL.*", "non_lesional", geopheno_checked$lesional)
geopheno_checked$lesional = gsub(".*gen.*", "NA", geopheno_checked$lesional)
geopheno_checked$lesional = geopheno_checked$diagnosis
geopheno_checked$lesional = NA
geopheno_checked$lesional = "lesional"
#eopheno_checked$lesional = "non_lesional"
geopheno_checked$lesional = gsub(".*non.*", "ctrl", geopheno_checked$lesional)


geopheno_checked$treatment
geopheno_checked$treatment= NA
#geopheno_checked$treatment = gsub(".*P.*", "placebo", geopheno_checked$treatment)

geopheno_checked$dose
#geopheno_checked$dose= geopheno$combined
#geopheno_checked$dose = 10
#geopheno_checked$dose = gsub(".*80.*", "80", geopheno_checked$dose)
#geopheno_checked$dose = gsub(".*P.*", "", geopheno_checked$dose)

geopheno_checked$dose_unit
geopheno_checked$dose_unit = "mg"
#geopheno_checked$dose_unit = gsub(".*10.*", "mg", geopheno_checked$dose_unit)

geopheno_checked$time_point
geopheno_checked$time_point = gsub(".*24.*", 1, geopheno_checked$time_point)
geopheno_checked$time_point = gsub(".*wk10.*", "70", geopheno_checked$time_point)
geopheno_checked$time_point = gsub(".*12_wk.*", "84_day", geopheno_checked$time_point)
geopheno_checked$time_point = gsub(".*16_wk.*", "112_day", geopheno_checked$time_point)
geopheno_checked$time_point = gsub(".*wk0.*", "0", geopheno_checked$time_point)
geopheno_checked$time_point = gsub(".*WK 1.*", "7_day", geopheno_checked$time_point)
geopheno_checked$time_point = gsub(".*wk2.*", "14", geopheno_checked$time_point)
geopheno_checked$time_point = gsub(".*3_wk.*", "21_day", geopheno_checked$time_point)
geopheno_checked$time_point = gsub(".*4_wk.*", "28_day", geopheno_checked$time_point)
geopheno_checked$time_point = gsub(".*6_wk.*", "42_day", geopheno_checked$time_point)
geopheno_checked$time_point = gsub(".*8_wk.*", "56_day", geopheno_checked$time_point)

geopheno_checked$time_point_unit = gsub(".*_.*", "day", geopheno_checked$time_point_unit)

geopheno_checked$tissue
geopheno_checked$tissue = geopheno$title
geopheno_checked$tissue = gsub(".*ski.*", "skin", geopheno_checked$tissue)
geopheno_checked$tissue = gsub(".*Uri.*", "urine", geopheno_checked$tissue)
geopheno_checked$tissue = gsub(".*synfluid.*", "synovial_fluid", geopheno_checked$tissue)
geopheno_checked$tissue = "blood"

geopheno_checked$anatomical_site 
geopheno_checked$anatomical_site = NA
#geopheno_checked$anatomical_site = geopheno_modified$anatomical_site
geopheno_checked$anatomical_site = gsub("_.*", "", geopheno_checked$anatomical_site)
#geopheno_checked$anatomical_site = "trunk_or_extremities"

geopheno_checked$array
geopheno_checked$array = NA

geopheno_checked$dye
geopheno_checked$dye = "biotin"
geopheno_checked$dye = "SYBR_green"
geopheno_checked$dye = "Cy3_and_Cy5"

geopheno_checked$slide
geopheno_checked$slide = NA

geopheno_checked$gender
geopheno_checked$gender = geopheno$Sex.ch1
geopheno_checked$gender = gsub(".*F.*", "f", geopheno_checked$gender)
geopheno_checked$gender = gsub(".*M.*", "m", geopheno_checked$gender)
#geopheno_checked$gender[189:208] = NA
#geopheno_checked$gender = geopheno_modified$gender

geopheno_checked$age
#geopheno_checked$age = geopheno$characteristics_ch1.2
geopheno_checked$age = gsub("y", "", geopheno_checked$age)
#geopheno_checked$age= "<6"

geopheno_checked$ethnicity
#geopheno_checked$ethnicity = gsub(".*African.*", "african_american", geopheno_checked$ethnicity)

# add id for each patient
patient_id = geopheno$source_name_ch1
patient_id = 1:15
patient_id = NA
patient_id = gsub(".*_", "", patient_id)
patient_id = gsub("(^|[^0-9])0+", "\\1", patient_id, perl = TRUE)
geopheno$individual.identifier.ch1 = gsub("C", 2, geopheno$individual.identifier.ch1)
length(unique(geopheno_checked$patient_id))
patient = unique(geopheno_checked$patient_id)
n=1:59
for(i in 1:length(n)){
  geopheno_checked$patient_id = gsub(patient[i], n[i], geopheno_checked$patient_id)
}

geopheno_checked = cbind(geopheno_checked, patient_id)

#missing data marked as NA
geopheno_checked$dose[which(geopheno_checked$dose =="")] =NA
geopheno_checked$dose_unit[which(geopheno_checked$dose_unit =="")]=NA
geopheno_checked$treatment[which(geopheno_checked$treatment =="")] =NA
geopheno_checked$time_point[which(geopheno_checked$time_point =="")] =NA
geopheno_checked$time_point_unit[which(geopheno_checked$time_point_unit =="")] =NA
geopheno_checked$anatomical_site[which(geopheno_checked$anatomical_site =="")] =NA
geopheno_checked$SCORAD[which(geopheno_checked$SCORAD =="")] =NA

View(geopheno_checked)

#save harmonized data (geopheno_checked)
dir.create(paste(checked_folder, "/", gseid, sep = ""))
check_folder = paste(checked_folder, "/", gseid, sep = "")
save(geopheno_checked, file = paste(check_folder, "/", gseid, ".RData", sep = ""))
write.table(x=geopheno_checked, file = paste(check_folder, "/", gseid, ".txt", sep = ""), sep = "\t", quote = FALSE)

