
#' from_csv_to_count_matrices function takes csv and txt outputs of variant callers and smatools depth and returns a matrix x of mutation counts and a matrix d of related read depth
#' @param lsdt a mandatory data.frame of metadata containing the list of sample names (column Sample) and sampling dates (column Sampling.date)
#' @param path.to.csv a path indacting were .csv, .tsv and/or .txt files (output of variant caller and read depth) are located 
#' @param length.seq the length of the sequence (usually a genome) studied


from_csv_to_count_matrices <- function(lsdt, path.to.csv, length.seq) {

  nT = nrow(lsdt) # number of time points
  
  read.depth = matrix(NA, length.seq, nT)
  for (t in 1:nT) {
    read.depth[,t] =  read.table(paste0(path.to.csv, lsdt$Sample[t], 
                                        ".primertrimmed.rg.sorted.bam.samtoolsdepth.txt"), 
                                 header=FALSE)[,3] 
  }
  # EXTRACT x matrix ####
  # mut.t : list of mutations present in each sample and number of times they are read (mut.t)
  mut.t = vector('list', nT)
  for (t in 1:nT) {
    tmp = read.table(paste0(path.to.csv, lsdt$Sample[t], 
                            ".primertrimmed.rg.sorted.bam.tsv"), 
                     header=TRUE)
    # keep only substitutions
    tmp = tmp[tmp$ALT %in% c("A", "T", "C", "G"),]
    mut.t[[t]] = tmp$ALT_DP
    names(mut.t[[t]]) = paste0(tmp$REF, tmp$POS, tmp$ALT)
  }
  
  # retrieve all mutations present in the dataset (all samples pooled)
  tmp = unique(unlist(lapply(mut.t, names)))
  # order mutations according to their positino on the genome
  tmp = tmp[order(readr::parse_number(tmp))]
  mut = tmp
  
  # create matrices for mutation counts
  n = length(mut) # number of mutations, all samples pooled
  x = matrix(0, n, nT, dimnames = list(mut,as.character(as.Date(lsdt$`Sampling date`))))
  for (t in 1:nT) {
    idx = intersect(mut, names(mut.t[[t]]))
    x[idx,t] = mut.t[[t]][idx]
  }
  
  # create matrices of read depth
  d = read.depth[readr::parse_number(mut),]
  dimnames(d) = dimnames(x)
  
  # # check si aucun nombre de mutation supérieur à la profondeur de lecture ... 
  # sum(x/d>1, na.rm = TRUE)
  
  data = list(x=x, d=d)
  return(data)
  
}


#' reduce_data function takes a dataset in required format and returns a reduced dataset according to selected period of time and thresholds. 
#' @param data (a list composed of a matrix x of mutation counts and a matrix d of read depth of the same size (mandatory)
#' @param time.period (vector of character string, two values): start date and end date of the time period (default: NULL for the entire time period)
#' @param threshold.freq (single value): frequency threshold (default is 0.05)
#' @param prop (single value): proportion of samples required to meet threshold.freq (default is NULL for one sample among all)
#' @param threshold.depth (single value): threshold on read depth (default is 10)
#' @param reduce.profile (list): (default is NULL for no reduction on mutation profile). reduce.profile is a list of:
#'  lineages: vector of names of lineages required to meet reduce.profile$threshold 
#'  threshold: single value for the minimal probability that a mutation belong to chosen lineages
#'  profile.matrix: mutation x lineage matrix containing the probability for a mutation to belong to a lineage
  

reduce_data = function(data=data, time.period = NULL, threshold.freq = 0.05, prop = NULL, threshold.depth = 10, reduce.profile = NULL) {
  
  x = data$x
  d = data$d 
  
  if(!is.null(time.period)) {
    x = x[,as.Date(colnames(x)) >= as.Date(time.period[1]) & as.Date(colnames(x)) <= as.Date(time.period[2])]
    d = d[,as.Date(colnames(d)) >= as.Date(time.period[1]) & as.Date(colnames(d)) <= as.Date(time.period[2])]
  }

  # set a read depth below "threshold.depth" to zero along with associated mutation count
  d[d<=threshold.depth] = 0
  x[d==0] = 0

  # remove lines of zero mutation count
  idx = apply(x, 1, sum)!=0
  d = d[idx,]
  x = x[idx,]
  #check sum(d-x<0); which(d-x<0, arr.ind = TRUE)

  # remove mutations below minimal frequency "threshold.freq" and available data (read depth above zero) in a proportion "prop" of the samples 
  ## if "prop" is null, minimal frequency "threshold.freq" in at least one sample
  if (is.null(prop)) {
    # at least on sample with freq above threshold.freq
    idx = apply(x/d, 1, max, na.rm=TRUE) > threshold.freq
    x = x[idx,]
    d = d[idx,]

  } else {
    # otherwise a proportion prop of the samples meet threshold.freq ()
    idx = apply(x/d, 1, function(z) {
      tmp = sort(z, na.last = FALSE)
      tmp[is.na(tmp)] = 0 # NA samples (read depth at zero) are considered as not meeting criteria
      quantile(tmp, probs = 1-prop)
    }
    ) > threshold.freq
    x = x[idx,]
    d = d[idx,]
  }
  
  # remove mutations of constant frequency equal 1 for them to be not interesting and are not supported for the E phase of the EM algorithm 
  idx = apply(d-x, 1, sum) != 0
  x = x[idx,]
  d = d[idx,]

  # Remove lines of NA (depth equal zero) for all but 1 time point
  idx = apply(d!=0, 1, sum) > 1
  x = x[idx,]
  d = d[idx,]
  
  if(!is.null(reduce.profile)) {
    # keep only mutations with probability above "reduce.profile$threshold" to belong to at least one lineage among "reduce.profile$lineage" according to matrix reduce.profile$profile.matrix
    sel = names(which((apply(reduce.profile$mutation.profile[,reduce.profile$lineage], 1, max) > reduce.profile$threshold)))
    x = x[intersect(sel, rownames(x)),]
    d = d[intersect(sel, rownames(x)),]
  }

  return(list(x = x, d = d))
}







