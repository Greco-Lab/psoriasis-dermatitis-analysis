#remove outliers based on results on "outlier_detection_with_affyBatchObject.R"

library(imager)

## create table where outlier data is saved
colnames(utopia_affy) = c("GSE","platform","annotation_file","sample.n","outliers_in_3","outliers_in_3.n","outliers_in_3.precentage","outliers_in_2","outliers_in_2.n","outliers_in_2.precentage")
write.table(x=utopia_affy, file = "/home/MOKA/veera/GEO/utopia_affy.txt" , sep = "\t", quote = FALSE)
save(utopia_affy, file = "/home/MOKA/veera/GEO/utopia_affy.RData")
comments = utopia_affy$outliers_in_3
comments = ""
utopia_affy = add_column(utopia_affy, comments)

#list all dataset identificators and go through every dataset one by one
all_gseid = as.vector(utopia_affy[,1])
gseid = all_gseid[1]
load(paste("/home/MOKA/veera/GEO/pheno_checked/", gseid, "/", gseid, ".RData", sep = ""))
utopia_affy$sample.n[grep(gseid, utopia_affy$GSE)] = dim(geopheno_checked)[1]

#load outlier data
load(paste("/home/MOKA/veera/GEO/pheno_checked/", gseid, "/RData/outliers.RData" , sep = ""))

M
colSums(M)
utopia_affy$all_outs.n[grep(gseid, utopia_affy$GSE)] = dim(M)[1] # same as length(outs_for_at_least_one_method)
# outs #filenames with two or more outs
# outs_for_at_least_one_method #filenames with at least one out
sum(M[,4]==3)
utopia_affy$outliers_in_3.n[grep(gseid, utopia_affy$GSE)] = sum(M[,4]==3)

M =as.data.frame(M)
three_outs = M[grep(3, M[,4]),] # matrix with three outs
two_outs = M[grep(2, M[,4]),] # matrix with two outs

## NUSE
NUSE_median = load.image(paste("/home/MOKA/veera/GEO/pheno_checked/", gseid, "/figures/NUSE/NUSE_median.png", sep=""))
plot(NUSE_median)
#list.files(paste("/home/MOKA/veera/GEO/pheno_checked/", gseid, "/figures/NUSE/", sep = ""))
NUSE_batch = load.image(paste("/home/MOKA/veera/GEO/pheno_checked/", gseid, "/figures/NUSE/NUSE_batch1.png", sep=""))
plot(NUSE_batch)
geopheno_checked[grep(paste(row.names(M[grep(1, M[,2]),]), collapse="|"),geopheno_checked$file_name),]
# NUSEout
# grep(paste(NUSEout, collapse = "|"), geopheno_checked$file_name)
# grep("GSM2859165", geopheno_checked$GSM)
# geopheno_checked

## RLE
RLE_median = load.image(paste("/home/MOKA/veera/GEO/pheno_checked/", gseid, "/figures/RLE/RLE_median.png", sep=""))
plot(RLE_median)
RLE_batch = load.image(paste("/home/MOKA/veera/GEO/pheno_checked/", gseid, "/figures/RLE/RLE_batch1.png", sep=""))
plot(RLE_batch)
geopheno_checked[grep(paste(row.names(M[grep(1, M[,1]),]), collapse="|"),geopheno_checked$file_name),]
# RLE_stat
# grep("GSM688772", geopheno_checked$GSM)
# RLEout


## RNADeg
RNADeg_slope = load.image(paste("/home/MOKA/veera/GEO/pheno_checked/", gseid, "/figures/RNADeg/RNADeg_slope.png", sep=""))
plot(RNADeg_slope)
RNADeg_batch = load.image(paste("/home/MOKA/veera/GEO/pheno_checked/", gseid, "/figures/RNADeg/RNADeg_batch1.png", sep=""))
plot(RNADeg_batch)
# RNADegouts = M[grep(1, M[,3]),]
# RNADeg[,grep("GSM3240039", colnames(RNADeg))]
geopheno_checked[grep(paste(row.names(M[grep(1, M[,3]),]), collapse="|"),geopheno_checked$file_name),]


## list the rows that should be removed from phenodata
#create folder for phenodata tables with outs removed
dir.create(paste("/home/MOKA/veera/GEO/pheno_checked/", gseid, "/outs_removed", sep = ""))
outs_removed_folder = paste("/home/MOKA/veera/GEO/pheno_checked/", gseid, "/outs_removed", sep ="")

# outliers in all three methods
# precentage of all of the samples
outliers_in_3.precentage = round(100 * length(row.names(three_outs))/dim(geopheno_checked)[1], digits = 1)
outliers_in_3.precentage
utopia_affy$outliers_in_3.precentage[grep(gseid, utopia_affy$GSE)] = outliers_in_3.precentage

# remove_rows = c()
remove_rows = grep(paste(row.names(three_outs), collapse =  "|"), geopheno_checked$file_name)
utopia_affy$outliers_in_3[grep(gseid, utopia_affy$GSE)] = paste(geopheno_checked$GSM[remove_rows], collapse = ", ")
geopheno_checked[remove_rows,]
outs_in_three=row.names(three_outs)

# remove list from above and save phenodata table
geopheno_3_outs_removed = geopheno_checked[-remove_rows,]
save(geopheno_3_outs_removed, file = paste(outs_removed_folder, "/", gseid, "_3_outs_removed.RData", sep = ""))
write.table(x=geopheno_3_outs_removed, file = paste(outs_removed_folder, "/", gseid, "_3_outs_removed.txt", sep = ""), sep = "\t", quote = FALSE)

# # combine geopheno and M 
# pheno_outs = geopheno_checked[grep(paste(outs, collapse="|"), geopheno_checked$file_name),]
# M_modified = add_column(as.data.frame(M), rownames(M), .before = "RLE")
# colnames(M_modified) = c("file_name", "RLE", "NUSE", "DEG", "SUM")
# pheno_outs = merge(pheno_outs, M_modified, by="file_name", .after = "file_name") # data frame with two or more outs
# pheno_outs_deg = pheno_outs[grep(1, pheno_outs$DEG),] # data frame with only the ones with DEG out with at least one more out
# View(pheno_outs)


## RNADeg + out
over_two_outs = rbind(three_outs, two_outs)
deg_and_three = over_two_outs[grep(1, over_two_outs[,3]),]
#deg_and_three = over_two_outs
utopia_affy$outliers_in_2.n[grep(gseid, utopia_affy$GSE)] = dim(deg_and_three)[1]
remove_rows_2 = grep(paste(row.names(deg_and_three), collapse = "|"), geopheno_checked$file_name)
utopia_affy$outliers_in_2[grep(gseid, utopia_affy$GSE)] = paste(geopheno_checked$GSM[remove_rows_2], collapse = ", ")
geopheno_checked[remove_rows_2,]

remove_rows_2 = c(1,3)

#precentage of all of the samples
outliers_in_2.precentage = round(100 * length(remove_rows_2)/dim(geopheno_checked)[1], digits = 1)
outliers_in_2.precentage
utopia_affy$outliers_in_2.precentage[grep(gseid, utopia_affy$GSE)] = outliers_in_2.precentage

# remove samples that are outliers in two methods including RNADeg and save phenodata
geopheno_2_outs_removed = geopheno_checked[-remove_rows_2,]
save(geopheno_2_outs_removed, file = paste(outs_removed_folder, "/", gseid, "_2_outs_removed.RData", sep = ""))
write.table(x=geopheno_2_outs_removed, file = paste(outs_removed_folder, "/", gseid, "_2_outs_removed.txt", sep = ""), sep = "\t", quote = FALSE)


# save utopia_affy table, where all outliers are saved for each dataset
#utopia_affy$comments[grep(gseid, utopia_affy$GSE)] = ""
View(utopia_affy)
write.table(x=utopia_affy, file = "/home/MOKA/veera/GEO/utopia_affy.txt" , sep = "\t", quote = FALSE)
save(utopia_affy, file = "/home/MOKA/veera/GEO/utopia_affy.RData")

#if no outliers were removed
utopia_affy$outliers_in_2.precentage[grep(gseid, utopia_affy$GSE)] = 0
utopia_affy$outliers_in_2.n[grep(gseid, utopia_affy$GSE)] = 0
utopia_affy$comments[grep(gseid, utopia_affy$GSE)] = "no outs removed"
