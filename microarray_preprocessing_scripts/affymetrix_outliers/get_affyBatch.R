#get affyBatch object

get_affyBatchObject <- function(celDir, fileNames, pheno, cdfname){
  celFilePaths <- file.path(celDir, fileNames)
  err <- 0
  tryCatch(simpleaffy::setQCEnvironment(cdfname), error=function(e){err<<-1})
  if(!err){
    affyBatchObject <- affy::read.affybatch(filenames=celFilePaths, phenoData=pheno, cdfname=cdfname)
  }else{
    headdetails <- affyio::read.celfile.header(celFilePaths[1])
    ref.cdfName <- headdetails[[1]]
    dim.intensity <- headdetails[[2]]
    rm.mask <- FALSE
    rm.outliers <- FALSE
    rm.extra <- FALSE
    verbose <- FALSE
    exprs <- affyio::read_abatch(celFilePaths, rm.mask, rm.outliers, rm.extra, ref.cdfName, dim.intensity[c(1, 2)], verbose)
    colnames(exprs) <- sampleNames(pheno)
    scandates <- sapply(seq_len(length(celFilePaths)), function(i) {
      sdate <- affyio::read.celfile.header(celFilePaths[i], info = "full")[["ScanDate"]]
      if (is.null(sdate) || length(sdate) == 0)
        NA_character_
      else sdate
    })
    protocol <- new("AnnotatedDataFrame", data = data.frame(ScanDate = scandates,
                                                            row.names = sampleNames(pheno), stringsAsFactors = FALSE),
                    dimLabels = c("sampleNames", "sampleColumns"))
    #protocol <- new("AnnotatedDataFrame")
    description <- new("MIAME")
    notes <- ""
    if (is.null(cdfname)) cdfname <- ref.cdfName
    affyBatchObject <- new("AffyBatch", exprs = exprs, cdfName = cdfname,
                           phenoData = pheno, nrow = dim.intensity[2], ncol = dim.intensity[1],
                           annotation = cleancdfname(cdfname, addcdf = FALSE),
                           protocolData = protocol, description = description,
                           notes = notes)
  }
  return(affyBatchObject)
}
