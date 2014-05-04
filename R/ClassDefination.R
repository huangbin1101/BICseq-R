# A class that contains short reads of a sample and a reference
# from one chromosome
# 
# Copyright 2010, Jianhua Zhang, all rights reserved
#

setClass("Aligned", representation(sAlign = "character",
                                  rAlign = "character",
                                  name = "character",
                                  chromCol = "integer",
                                  locCol = "integer"))

setGeneric("segAligned", function(object, bin, lambda, coverage, ...)
           standardGeneric("segAligned"))
setMethod("segAligned", "Aligned", 
          function(object, bin, lambda, coverage, ...){
              if(tolower(coverage) %in% c("deep", "high")){
                  return(segBowtie(sAlign(object), rAlign(object), name(object), 
                      locCol(object), bin, lambda, ...))
              }else{
                  return(segBowtieLite(sAlign(object), rAlign(object), 
                      sampleName = name(object), chromCol = chromCol(object),
    		      locCol = locCol(object), binSize = bin, lambda = lambda, ...))
              }})
              
              
setGeneric("sAlign", function(object)
           standardGeneric("sAlign"))
setMethod("sAlign", "Aligned", 
          function(object) object@sAlign)
          
setGeneric("rAlign", function(object)
           standardGeneric("rAlign"))
setMethod("rAlign", "Aligned", 
          function(object) object@rAlign)   
          
setGeneric("name", function(object)
           standardGeneric("name"))
setMethod("name", "Aligned", 
          function(object) object@name)
          
setGeneric("chromCol", function(object)
           standardGeneric("chromCol"))
setMethod("chromCol", "Aligned", 
          function(object) object@chromCol)
          
setGeneric("locCol", function(object)
           standardGeneric("locCol"))
setMethod("locCol", "Aligned", 
          function(object) object@locCol)
              

setClass("BICseq", representation(samp = "character",
                                  ref = "character",
                                  seqNames = "character"))

setGeneric("getBICseg", function(object, bin = 100, lambda = 2, 
    winSize = 200,quant = 0.95, mult = 1)
           standardGeneric("getBICseg"))
setMethod("getBICseg", "BICseq", 
          function(object, bin = 100, lambda = 2, 
              winSize = 200, quant = 0.95, mult = 1){
              segged <- segBAM(seqNames(object), samp(object), ref(object), bin = bin, 
                  lambda = lambda, winSize = winSize, quant = quant, mult = mult)
              new("BICsegment", bin = segged[["bin"]], seg = segged[["seg"]])
              })

setGeneric("samp", function(object)
           standardGeneric("samp"))
setMethod("samp", "BICseq",
          function(object) object@samp)
          
setGeneric("ref", function(object)
           standardGeneric("ref"))
setMethod("ref", "BICseq",
          function(object) object@ref)
          
setGeneric("seqNames", function(object)
           standardGeneric("seqNames"))
setMethod("seqNames", "BICseq",
          function(object) object@seqNames)
 
setMethod("show", "BICseq",
          function(object) {
            cat("Object of BICseq\n")
            cat(paste("Sample BAM file: ", 
                 samp(object), 
                 "\n"), sep = "")            
            cat(paste("Reference BAM file: ", 
	                     ref(object), "\n"), sep = "")
	    cat("Sequence names:\n")
	    cat(paste(seqNames(object), collapse = ";"))
	    cat("\n")
          })
 
setClass("BICsegment", representation(bin = "GRangesList",
                                  seg = "GRangesList"))
setGeneric("bin", function(object)
           standardGeneric("bin"))
setMethod("bin", "BICsegment",
          function(object) object@bin)
          
setGeneric("seg", function(object)
           standardGeneric("seg"))
setMethod("seg", "BICsegment",
          function(object) object@seg)
                         
setGeneric("ratios", function(object)
           standardGeneric("ratios"))
setMethod("ratios", "BICsegment",
          function(object){
              list(bin = getRatios(bin(object), "bin"), seg = getRatios(seg(object), "seg"))
          })

setGeneric("getSummary",function(object, correction=TRUE) standardGeneric("getSummary"))
setMethod("getSummary","BICsegment",
	function(object){
		bicseg = getSummary.BICseq(object,correction=TRUE,cutoff=0.3)
		return(bicseg)
	})

                                  
setGeneric("plot", function(x, y, ...)
           standardGeneric("plot"))
setMethod("plot", "BICsegment",
          function(x, y, ...) {
             args <- list(...)
             args[["bicsegment"]] <- x
             do.call("plotBICseq", args = args)
             #plotBICseq(x, sampleName = args[["sampleName"]], save = args[["save"]], 
             #endrule = args[["endrule"]], k = args[["k"]], 
              #indexOnly = args[["indexOnly"]])
          })


setMethod("show", "BICsegment",
          function(object) {
            cat("Object of BICsegment\n")
            cat(paste("Number of bins:", 
                 nrow(bin(object)), 
                 "\n", sep = " "))
            print(bin(object)[1:5, ])
            cat("..........\n")
            
            cat(paste("Number of segments:", 
	                     nrow(seg(object)), "\n", sep = ""))
	    print(seg(object)[1:5, ])
            cat("..........\n")
          })
