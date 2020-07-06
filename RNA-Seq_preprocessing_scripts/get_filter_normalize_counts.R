### Get, filter and normalise counts.

get_raw_counts <- function(bam.files, annotation, paired.end, ncores=2) {
  
  if(is.null(bam.files)){stop("Error: please provide one or more .bam files!")}
  if(is.null(annotation)){stop("Error: please provide an annotation file in .gtf or .gff format!")}
  if(!class(paired.end)=="logical"){stop("Error: indicate whether the reads are paired end")}
  if(!class(ncores)=="numeric"){stop("Error: please input a numeric value!")}
  
  raw.counts1 <- Rsubread::featureCounts(bam.files, annot.ext=annotation, isPairedEnd=paired.end, isGTFAnnotationFile = TRUE, nthreads = ncores)
  
  total.counts.matrix <- raw.counts1[[1]]
  
  genelengths <- raw.counts1$annotation$Length
  names(genelengths) <- raw.counts1$annotation$GeneID
  
  
  return(list(total.counts.matrix, genelengths))
}




filter_low_counts <- function(counts.matrix, conditions, method = "cpm", normalized=FALSE, depth=NULL, cpm=1, p.adj = "fdr"){
  
  if(is.null(counts.matrix)){stop("Error: please provide a numeric count matrix!")}
  if(is.null(conditions)){stop("Error: please provide a factor or a vector indicating the conditions!")}
  if(!method %in% c("cpm", "wilcoxon", "proportion")) {stop("Error: Please type in one of valid methods!")}
  
  
  if (method=="cpm"){
    filtered.counts = NOISeq::filtered.data(counts.matrix, factor = conditions, norm = normalized, method = 1, cv.cutoff = 100, cpm = cpm, p.adj = p.adj)
    
  }else if(method=="wilcoxon"){
    filtered.counts = NOISeq::filtered.data(counts.matrix, factor = conditions, norm = normalized, method = 2, cv.cutoff = 100, p.adj = p.adj)
    
  }else if(method=="proportion"){
    
    if(is.null(depth)){stop("Error: indicate a numeric vector indicating per sample library depths")}
    if(!class(depth)=="numeric"){stop("Error: please provide the depth argument with a numeric vector!")}
    ### Compute librarary depth
    
    
    filtered.counts = NOISeq::filtered.data(counts.matrix, factor = conditions, norm = normalized, depth = depth, method = 3, cv.cutoff = 100, cpm = cpm, p.adj = p.adj)
    
  }
  return(filtered.counts)
}





normalize_counts <- function(filtered.counts, method="UQUA", gene.lengths=NULL, length.correction=FALSE){
  
  if(is.null(filtered.counts)){stop("Error: please provide a matrix of raw or filtered counts")}
  if(!method %in% c("RPKM", "UQUA", "TMM", "CPM")) {stop("Error: please type in one of the allowed methods!")}
  if(method == "RPKM" & is.null(gene.lengths)){stop("Error: please provide a vector of transcripts lengths!")}
  
  if (method=="RPKM"){
    normalized.counts = NOISeq::rpkm(filtered.counts, long = gene.lengths, k = 0, lc = 1)
  }else if(method=="UQUA"){
    if(length.correction==TRUE){
      normalized.counts = NOISeq::uqua(filtered.counts, long = gene.lengths, lc = 0, k = 0)
    }else{
      normalized.counts = NOISeq::uqua(filtered.counts, lc = 0, k = 0)
    }
    
  }else if(method=="TMM"){
    if(length.correction==TRUE){
      normalized.counts = NOISeq::tmm(filtered.counts, long = gene.lengths, lc = 0)
    }else{
      normalized.counts = NOISeq::tmm(filtered.counts, long = 1000, lc = 0)
    }
    
  }else if(method=="CPM"){
    gene.lengths=1000
    normalized.counts = edgeR::cpm(filtered.counts, long = gene.lengths, k = 0, lc = 1)
  }
  
  
  return(normalized.counts)
}
