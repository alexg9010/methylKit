
###  S3 Functions to be used in S4 methods

# This function will summarize overlaps of a methylRaw or methylRawList with a Biomarker GRanges
# and return a new GRanges or data.frame of the same length as Biomarker.
# The reason for this option is to preserve alignment, 
# so that multiple such objects can then be added together (since for all k, 
# the k-th entry of one such object refers to the same region as the k-th region of the other). 
.markerCounts <- function(object, regions,
                          strand.aware=TRUE,  
                          naToZero=TRUE){
  #require(GenomicRanges)
  
  # sort regions
  regions <- sortSeqlevels(regions)
  regions <- sort(regions,ignore.strand=TRUE)
  
  # overlap object with regions
  # convert object to GRanges
  if(!strand.aware){
    g.meth=as(object,"GRanges")
    strand(g.meth)="*"
    mat=IRanges::as.matrix( findOverlaps(regions,g.meth ) )
    #mat=matchMatrix( findOverlaps(regions,g.meth ) )
    
  }else{
    mat=IRanges::as.matrix( findOverlaps(regions,as(object,"GRanges")) )
    #mat=matchMatrix( findOverlaps(regions,as(object,"GRanges")) )
    
  }
  
  #require(data.table)
  # create a temporary data.table row ids from regions and counts from object
  dt=data.table(id = mat[, 1], getData(object)[mat[, 2], 5:ncol(object) ])
  
  # 
  # ## initialize some variables
  # coverage=.SD=numTs=id=numTs1=covered=NULL
  
  # use data.table to sum up counts per region
  dt=dt[,lapply(.SD,sum),by=id]
  
  # look for values with "name" in it, eg. "tx_name" or "name"
  # valuesList = names(values(regions))
  # nameid = valuesList[grep (valuesList, pattern="name")]
  
  n <- length(regions)
  dt <- merge(data.table(id=seq_len(n)),dt,all=TRUE)
  
  NAToZero = function(DT) {
    for (i in names(DT))
      DT[is.na(get(i)), (i):=0]
  }
  
  if(naToZero) NAToZero(dt)
  
  
  temp.df=as.data.frame(regions) # get regions to a dataframe
  #create a new methylBase object to return
  new.data=data.table(#id      =new.ids,
    chr     =temp.df[dt$id,"seqnames"],
    start   =temp.df[dt$id,"start"],
    end     =temp.df[dt$id,"end"],
    strand  =temp.df[dt$id,"strand"],
    as.data.frame(dt[,c(2:(ncol(dt))),with=FALSE]),stringsAsFactors=FALSE)
  
  new("methylRaw",new.data,sample.id=object@sample.id,
      assembly=object@assembly,context=object@context,
      resolution="region")
  
}



# End of S3 functions -----------------------------------------------------






#' Count overlap of object with Biomarker
#' 
#' This function will summarize overlaps of a methylRaw or methylRawList with a Biomarker GRanges
#' and return a new GRanges or data.frame of the same length as Biomarker.
#' The reason for this option is to preserve alignment, 
#' so that multiple such objects can then be added together (since for all k, 
#' the k-th entry of one such object refers to the same region as the k-th region of the other). 
#' 
#' 
#'
#' @param object either methylRaw or methylRawList,the genomic data that we are 
#' comparing to a set of biomarkers
#' @param regions GRanges the biomarker regions
#' @param strand.aware if set to TRUE only CpGs that match the strand of 
#' the region will be summarized. (default:FALSE)
#' @param naToZero should actual NAs be converted to zeros (default: TRUE)
#' @param as.granges return as GRanges instead of data.frame (default:TRUE)
#'
#' @return either GRanges or data.frame
#' @export
#'
#' @docType methods
setGeneric("markerCounts", 
           function(object, regions,
                    strand.aware=TRUE,  
                    naToZero=TRUE, as.granges=TRUE)
  standardGeneric("markerCounts") )


# GETs marker counts for given GRanges object
# RETURNS a new methylRaw object
setMethod("markerCounts", signature(object="methylRaw",regions="GRanges"),
          function(object, regions,
                   strand.aware=TRUE,  
                   naToZero=TRUE, as.granges=TRUE){
            
            res <- .markerCounts(object, regions,
                                 strand.aware,  
                                 naToZero) 
            
            ##not return object as methylRaw, since it is not campatible to other methods
            res <- getData(res)
            
            if(as.granges) { 
              return(makeGRangesFromDataFrame(res,keep.extra.columns = TRUE))
            } else {
              return(res)
            }
            
          }
)


# GETs marker counts for given GRanges object
# RETURNS a new methylRaw object
setMethod("markerCounts", signature(object="methylRaw",regions="GRangesList"),
          function(object, regions,
                   strand.aware=TRUE,  
                   naToZero=TRUE, as.granges=TRUE){
            
            # combine and sort GRanges from List
            regions <- unlist(regions)
            
            res <- .markerCounts(object, regions,
                                 strand.aware,  
                                 naToZero) 
            
            ##not return object as methylRaw, since it is not campatible to other methods
            res <- getData(res)
            
            if(as.granges) { 
              return(makeGRangesFromDataFrame(res,keep.extra.columns = TRUE))
            } else {
              return(res)
            }
            
          }
)

# GETs marker counts for given GRanges object
# RETURNS a new methylRaw object
setMethod("markerCounts", signature(object="methylRawList",regions="GRanges"),
          function(object, regions,
                   strand.aware=TRUE,  
                   naToZero=TRUE, as.granges=TRUE){
            
            outList=list()
            for(i in 1:length(object))
            {
              obj <- .markerCounts(object[[i]], regions,
                                   strand.aware,  
                                   naToZero) 
              outList[[i]] = obj
            }
            myobj=new("methylRawList", outList,treatment=object@treatment)
            
            res <- unite(myobj)
            
            ##not return object as methylBase, since it is not campatible to other methods
            res <- getData(res)
            
            if(as.granges) { 
              return(makeGRangesFromDataFrame(res,keep.extra.columns = TRUE))
            } else {
              return(res)
            }
          }
)

# GETs marker counts for given GRanges object
# RETURNS a new methylRaw object
setMethod("markerCounts", signature(object="methylRawList",regions="GRangesList"),
          function(object, regions,
                   strand.aware=TRUE,  
                   naToZero=TRUE, as.granges=TRUE){
            
            # combine and sort GRanges from List
            regions <- unlist(regions)
            
            outList=list()
            for(i in 1:length(object))
            {
              obj <- .markerCounts(object[[i]], regions,
                                   strand.aware,  
                                   naToZero) 
              outList[[i]] = obj
            }
            myobj=new("methylRawList", outList,treatment=object@treatment)
            
            res <- unite(myobj)
            
            ##not return object as methylBase, since it is not campatible to other methods
            res <- getData(res)
            
            if(as.granges) { 
              return(makeGRangesFromDataFrame(res,keep.extra.columns = TRUE))
            } else {
              return(res)
            }
          }
)
          
          
          
