BICseq <- function(sample, reference, 
    seqNames = paste("chr", c(1:22, "X", "Y"), sep = "")){
	
	if(!file.exists(sample)){stop(paste("No such file: \""), sample,"\"",sep="")}
	if(!file.exists(reference)){stop(paste("No such file: \"",reference, "\"",sep=""))}

	header.sample = scanBamHeader(sample)
	target.sample = names(header.sample[[1]]$targets)
	ind.match = match(seqNames,target.sample)
	if(sum(is.na(ind.match))>0){
		ind.tmp = which(is.na(ind.match))
		stop(paste("No such chromosome \"",seqNames[ind.tmp[1]], "\" in the header of the BAM file \"", sample,"\"",sep=""))
		}
	header.reference = scanBamHeader(reference)
	target.reference = names(header.reference[[1]]$targets)
	ind.match = match(seqNames,target.reference)
	if(sum(is.na(ind.match))>0){
		ind.tmp = which(is.na(ind.match))
		stop(paste("No such chromosome \"",seqNames[ind.tmp[1]], "\" in the header of the BAM file \"", reference, "\"",sep=""))
		}

	return(new("BICseq", samp = sample, ref = reference, seqNames = seqNames))
}

Aligned <- function(sAlign, rAlign, name, chromCol = 2, locCol = 3){
    return(new("Aligned", sAlign = sAlign, rAlign = rAlign, name = name,
        chromCol = as.integer(chromCol), locCol = as.integer(locCol)))
}

## binSize - Size of windows (number of base pair) for counting reads
## lambda -  
## win_size: to determine whether read counts in a region is an outlier, its neighboring win-size windows will be
##     selected to calculate the quantile of reads in the windows. if the region has more than mul*q number of 
##     reads, it will be viewed as an outlier and the number of reads in this region will be set as mul*q. This number 
##     should be integer and larger than 1.
## quantile:  quantile to use (say 0.95th percentile) for p. This number should be between 0 and 1.
## mul: if a genomic position has more than mul*q reads, it will be viewed as an outlier, where q is the 
##    quantile calculated based on win_size and quantile.

segReads <- function(sample, reference, seqName, binSize, lambda, 
    winSize = 200, quant = 0.95, mult = 1, outFile = NULL){
    keep <- TRUE
    if(is.null(outFile)){
        keep <- FALSE
        outFile <- tempfile()
    }
    on.exit(ifelse(keep, break, unlink(outFile)))
    
    bins <-  .C("sort_rms_binning", as.integer(sample), length(sample), as.integer(reference), 
        length(reference), as.integer(binSize), as.integer(winSize), as.double(quant), 
        as.double(mult), as.character(outFile))
    #bins <- .C("binning", as.integer(sample), as.integer(reference),
    #    length(sample), length(reference), as.integer(binSize), 
    #    as.character(tempFile), PACKAGE = "BICseq")
        
    #calculates the log2 ratio for each bin
    bins <- read.delim(outFile, sep = "\t", header = FALSE, as.is = TRUE)
    #dat <- getRatios(bins) 
    bins <- GRanges(seqnames = Rle(seqName, nrow(bins)), 
        ranges = IRanges(bins[, 4], bins[, 5]), sampleReads = bins[, 1], 
        totalReads = bins[, 2], sampleFreq = bins[, 3])
    
    storage.mode(outFile) <- "character"
    storage.mode(lambda) <- "double"
    segs <- .Call("bic_seq_Rinterface", outFile, lambda, PACKAGE = "BICseq")
    segs <- GRanges(seqnames = Rle(seqName, nrow(segs)), 
        ranges = IRanges(segs[, 4], segs[, 5]), sampleReads = segs[, 1], 
        totalReads = segs[, 2], sampleFreq = segs[, 3])
    
    #segs <- cbind(segs[, c("loc.start", "loc.end")], 
    #    num.mark = paste(segs[, "sampleReads"], (as.numeric(segs[, "totalReads"]) -
    #        as.numeric(segs[, "sampleReads"])), sep = "/"),
    #    seg.mean = apply(segs, 1, FUN = function(x) 
    #    round(log2((as.numeric(x["sampleReads"])/
    #    (as.numeric(x["totalReads"]) - as.numeric(x["sampleReads"])))/
    #    (length(sample)/length(reference))), 4)))
    
    return(list(bin = bins, seg = segs))
}

getRatios <- function(bins, what = c("bin", "seg")){

    if(class(bins) == "GRangesList"){
        bins <- IRanges::unlist(bins)
    }
    
    #chrom <- gsub("chr", "", as.character(seqnames(bins)))
    chrom <- as.character(seqnames(bins))
    sampleTotal = sum(as.numeric(elementMetadata(bins)[, "sampleReads"]))
    controlTotal = sum(as.numeric(elementMetadata(bins)[, "totalReads"])) - sampleTotal
    #adjust <- log2(sum(as.numeric(elementMetadata(bins)[, "sampleReads"]))/
    #    (sum(as.numeric(elementMetadata(bins)[, "totalReads"]))
    #    -sum(as.numeric(elementMetadata(bins)[, "sampleReads"]))))
    adjust <- log2(sampleTotal/controlTotal)
    if(what == "bin"){
        ratios <- round(.C("getratios", 
            as.integer(as.numeric(elementMetadata(bins)[, "sampleReads"])), 
            as.integer(as.numeric(elementMetadata(bins)[, "totalReads"]) - 
            as.numeric(elementMetadata(bins)[, "sampleReads"])), 
            as.integer(length(bins)), as.double(adjust), 
            as.double(rep(0, length(bins))))[[5]], 4)
    }else{    
        ratios <- round(log2(as.numeric(elementMetadata(bins)[, "sampleReads"])/
            (as.numeric(elementMetadata(bins)[, "totalReads"]) -
             as.numeric(elementMetadata(bins)[, "sampleReads"]))) -adjust, 4)
            
            
    }
    
    return(data.frame(chrom = chrom, start = start(ranges(bins)), 
        end = end(ranges(bins)), ratio = ratios))
}

getSummary.BICseq <- function(bicsegment, correction=TRUE,cutoff=0.3){ ## cutoff is used to remove the regions with extreme initial copy ratio estimates and recalcualte the copy ratios

    bins = seg(bicsegment)
    if(class(bins) == "GRangesList"){
        bins <- IRanges::unlist(bins)
    }

  chrom <- as.character(seqnames(bins))
  sampleReads = as.numeric(elementMetadata(bins)[, "sampleReads"])
  controlReads = as.numeric(elementMetadata(bins)[, "totalReads"]) - sampleReads
  sampleTotal = sum(sampleReads)
  controlTotal = sum(controlReads)
  adjust <- log2(sampleTotal/controlTotal)
  ratios = round(log2(sampleReads/controlReads)-adjust,4)

  ind.rm = which(abs(ratios)>cutoff)
  if(length(ind.rm)>0) {
	sampleTotal.tmp = sampleTotal - sum(sampleReads[ind.rm])
	controlTotal.tmp = controlTotal - sum(controlReads[ind.rm])

	if(sampleTotal.tmp>=100 && controlTotal.tmp>=100){ ## make sure there are enough "total" reads
		sampleTotal = sampleTotal.tmp
		controlTotal = controlTotal.tmp
		}

	adjust <- log2(sampleTotal/controlTotal)
	ratios = round(log2(sampleReads/controlReads)-adjust,4)
	}

  start =  start(ranges(bins))
  end = end(ranges(bins))
  genomesize = sum(as.numeric(end - start + 1))

  log.pvalue = rep(1,length(start))
  prob = sampleTotal/(sampleTotal+controlTotal)
  probN = 1-prob

  for(k in c(1:length(log.pvalue))){
	ttreads = sampleReads[k] + controlReads[k]
	p = sampleReads[k]/ttreads
	if(p<=prob){
		p = (sampleReads[k]+0.5)/ttreads
		log.pvalue[k] = pnorm((p-prob)*sqrt(ttreads),sd=sqrt(prob*(1-prob)),log.p=TRUE)
		}else{
		p = (controlReads[k]+0.5)/ttreads
		log.pvalue[k] = pnorm((p-probN)*sqrt(ttreads),sd=sqrt(probN*(1-probN)),log.p=TRUE)
		}
	
	}
   log10.pvalue = log.pvalue/log(10)

   ##bonferroni correction
   if(correction!=FALSE){
	   log10.pvalue = log10.pvalue + log10(genomesize)
	   ind = which(log10.pvalue>=0.0)
	   log10.pvalue[ind] = 0.0
	}


   return(data.frame(chrom=chrom,start=start,end=end,sampleReads=sampleReads,referenceReads=controlReads, ratio=ratios,log10.pvalue=log10.pvalue))
}


segBAM <- function(seqs, BAMSample, BAMReference, binSize = 100, 
    lambda = 2, winSize = 200, quant = 0.95, mult = 1){
    
    segged <- lapply(seqs, segSeq, BAMSample = BAMSample, 
        BAMReference = BAMReference, binSize = binSize, lambda = lambda, 
        winSize = winSize, quant = quant, mult = mult) 
    gc()
    names(segged) <- seqs
    bins <- do.call("GRangesList", args = lapply(segged, FUN = function(x)
        return(x[["bin"]])))
    segs <- do.call("GRangesList", args = lapply(segged, FUN = function(x) 
        return(x[["seg"]])))
        
    return(list(bin = bins, seg = segs))
}  

## ... 
segSeq <- function(seqName, BAMSample, BAMReference, binSize = 100, 
    lambda = 2, winSize = 200, quant = 0.95, mult = 1, ...){
 
    param <- ScanBamParam(what = c("pos"), which = GRanges(seqName, 
        IRanges(1, 500000000)), scanBamFlag(isUnmappedQuery = FALSE))
    sample <- scanBam(BAMSample, ..., param = param)[[1]][["pos"]]
    ref <- scanBam(BAMReference, ..., param = param)[[1]][["pos"]]

    if(length(sample)==0){stop(paste("No read available for the chromosome \"", seqName, "\" in \"",BAMSample,"\"", sep=""))}
    if(length(ref)==0){stop(paste("No read available for the chromosome \"", seqName, "\" in \"",BAMReference,"\"", sep=""))}    

    segged <- segReads(sample, ref, seqName = seqName, binSize = binSize, 
        lambda = lambda, winSize = winSize, quant = quant, mult = mult)
        
    return(segged) 
}


plotBICseq <- function(bicsegment, sampleName, save = TRUE, 
    endrule = c("none", "median", "keep", "constant"), k = 3, indexOnly = FALSE,chromOrder=NULL,plotBin=TRUE){
  endrule <- match.arg(endrule)
  graphList <- list()
  
  #segments <- BICseq:::getRatios(seg(bicsegment), "seg")
  segments <- BICseq:::getSummary(bicsegment)
  segments <- segments[,c(1:3,6)]
  if(!is.null(chromOrder)){
	chroms.segments = as.character(segments[,"chrom"])
	chrom.order = as.character(chromOrder)
	ind = match(chrom.order,chroms.segments)
	if(sum(is.na(ind))){
		stop(paste("No such chromosome:",chrom.order[which(is.na(ind))[1]]))
		}
	}
  bins <- BICseq:::getRatios(bin(bicsegment), "bin")

  if(!is.null(chromOrder)){
	ind1 = match(as.character(bins[,"chrom"]),as.character(chromOrder))
	bins = bins[!is.na(ind1),]
	ind2 = match(segments[,"chrom"],as.character(chromOrder))
	segments = segments[!is.na(ind2),]
	}

  
  
  locations <- alignGenes(bins[, c("chrom", "start")], 
      indexOnly = indexOnly,chromOrder=chromOrder)
  adjustments <- getAdjustments(bins[, c("chrom", "start")], 
      indexOnly = indexOnly,chromOrder=chromOrder)
   
  if(save){
      graphList[[sampleName]] <- file.path(tempdir(),
                                      paste(sampleName, "Seg.png", sep = ""))
      png(filename = graphList[[sampleName]], width = 1200, height = 650)
  }
  max.locations = max(locations)
  graphics::plot(0, 0, type = "n", main = sampleName, xlab = "Chromsome",
         ylab = " Log2 ratio", ylim = c(-5, 5), axes = FALSE,
         xlim = c(0, max.locations + 10),cex=1.2)
  axis(2)
  box()
  highlightChrom(adjustments, -5, 5,max.locations)

  if(plotBin==TRUE){
	  if(endrule != "none"){
	      nas <- which(is.na(as.numeric(bins[, "ratio"])))
	      if(length(nas) != 0){
        	points(locations[-nas],
             	runmed(as.vector(bins[-nas, "ratio"]),
                                     k = k, endrule = endrule), cex = 0.4, pch = 16)
      		}else{
        	points(locations,
             	runmed(as.vector(bins[, "ratio"]),
                                     k = k, endrule = endrule), cex = 0.4, pch = 16)
      		}
  	}else{
      		points(locations,
             	as.numeric(as.numeric(as.vector(bins[, "ratio"]))),
             	cex = 0.4, pch = 16)
  	}
  }
  lines(c(min(locations), max(locations)), rep(0, 2), lwd = 2, col = "blue")
  lines(c(min(locations), max(locations)), rep(1, 2), col = "blue")
  lines(c(min(locations), max(locations)), rep(-1, 2), col = "blue")
  
  if(!indexOnly){
      if(is.vector(segments)){segments = t(segments)}                         
      drawSegs(segments, adjustments)
  }else{
      for(index in 1:nrow(segments)){
          positions <- range(locations[which(as.character(bins[, "chrom"]) == 
                               as.character(segments[index, "chrom"]) &
                             as.numeric(bins[, "start"]) < 
                               as.numeric(segments[index, "end"]) &
                             as.numeric(bins[, "end"]) >
                               as.numeric(segments[index, "start"]))])
          lines(c(positions[1], positions[2]), rep(segments[index, "ratio"], 2), 
            col = "red", lwd = 2)
      }
  }

  markChrom(adjustments, -5,max.locations)
  if(save){
      dev.off()
  } 
  if(save){
    return(graphList)
  }else{
    return(invisible())
  }
}

alignGenes <- function(positions, indexOnly = FALSE, chromOrder=NULL){
  if(indexOnly){
    return(1:nrow(positions))
  }else{
    adjustments <- getAdjustments(positions, indexOnly, chromOrder=chromOrder)
    for(chrom in names(adjustments)){
      positions[positions[, 1] == chrom, 2] <-
        as.numeric(positions[positions[, 1] == chrom, 2]) +
          as.numeric(adjustments[chrom])
    }
    return(as.numeric(positions[, 2]))
  }
}

# This function gets the values that can be used to adjust chromosomal
# locations to align genes on different chromosomes so that they  appear
# in sequence from chromosome one to Y along a single strand
getAdjustments <- function(positions, indexOnly = FALSE,chromOrder=chromOrder){
  if(indexOnly){
    temp <- split.data.frame(positions, factor(positions[, 1]))
    temp <- unlist(lapply(temp, FUN = function(x) nrow(x) + 1))
  }else{
    temp <- split.data.frame(positions, factor(positions[, 1]))
    temp <- unlist(lapply(temp, FUN = function(x) max(as.numeric(x[, 2]))))
  }
  # Chromosomes to 30 for other organisms (e. g. fish - 25)
  #chroms <- sort(as.numeric(names(temp)[names(temp) %in% 1:30]))
  chroms <- as.character(names(temp))
  if(!is.null(chromOrder)){
	  chrom.order = as.character(chromOrder)
	  ind = match(chrom.order,chroms)
	  if(sum(is.na(ind))>0){
		stop(paste("No such chromsome ",chrom.order[which(is.na(ind))[1]]))
		}
	  chroms = chrom.order	  
	}
  #if(any(names(temp) %in% "X")){
  #  chroms <- c(chroms, "X")
  #}
  #if(any(names(temp) %in% "Y")){
  #  chroms <- c(chroms, "Y")
  #}
  adjustments <- 0
  for(index in 1:length(chroms) - 1){
    adjustments <- c(adjustments, (adjustments[length(adjustments)]+
                                   temp[as.character(chroms[index])]))
  }
  names(adjustments) <- chroms
  return(adjustments)
}

highlightChrom <- function(adjustments_in, min, max,max.x){
  if(max.x < max(adjustments_in)){stop("Error in highlightChroom")}
  adjustments = c(adjustments_in,max.x)
  for(index in 1:(length(adjustments)-1)){
    if(index %% 2 == 1){
      polygon(c(adjustments[index], adjustments[index + 1],
                adjustments[index + 1], adjustments[index]),
              c(min, min, max, max), col = "gray", border = "white")
    }
  }
  return(invisible())
}

drawSegs <- function(segdata, adjustments, seqOnly = TRUE){
  drawSegLine <- function(segLocs){
    toAdd <- adjustments[as.character(as.vector(segLocs["chrom"]))]
    lines(c(as.numeric(segLocs["start"]) + toAdd,
            as.numeric(segLocs["end"]) + toAdd),
          rep(segLocs["ratio"], 2), col = "red", lwd = 2)
  }
  junck <- apply(segdata, 1, FUN = drawSegLine)
  return(invisible())
}

markChrom <- function(adjustments_in, min,max.x){
  chromLocs <- NULL
  chromNames <- NULL
  adjustments = c(adjustments_in,max.x)
  for(i in 1:length(adjustments) - 1){
    if(i %% 2 == 1){
      chromLocs <- c(chromLocs, mean(c(adjustments[i], adjustments[i + 1])))
      chromNames <- c(chromNames, names(adjustments)[i])
    }
  }
  text(chromLocs, rep(min - 0.125, length(chromLocs)), chromNames, cex = 1.2)
  return(invisible())
}


