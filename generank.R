# Functions for ranking a CNA dataset by CNA loss frequency and CNA gain frequency
# input is CNA data set:
#   Col 1: "Hugo_Symbol"
#   Col 2 - End: Sample names
#   Rows 1 - End: Feature names (Hugo Symbols for genes)
lossRank <- function(dataset){
  gene.names = dataset[,1]
  logical = (dataset < 0)
  loss_freq = data.frame()
  for (i in 1:nrow(logical)) {
    loss_freq[i,1] = gene.names[i]
    loss_freq[i,2] = sum(logical[i,2:ncol(logical)], na.rm = TRUE)
    loss_freq[i,3] = (ncol(logical)-1) - sum(logical[i,2:ncol(logical)], na.rm = TRUE)
    loss_freq[i,4] = sum(logical[i,2:ncol(logical)], na.rm = TRUE)/(ncol(logical)-1)
  }
  ordered = order(loss_freq[,4])
  x = c()
  y = c()
  z = c()
  a = c()
  b = c()
  for (i in 1:length(ordered)) {
    x[i] = i
    y[i] = loss_freq[ordered[i], 4]
    z[i] = gene.names[ordered[i]]
    a[i] = loss_freq[ordered[i], 2]
    b[i] = loss_freq[ordered[i], 3]
  }
  output = data.frame(z, x, a, b, y)
  colnames(output) = c("Hugo_Symbol", "Gene Rank", "total copy number loss", 
                       "total unaltered", "copy number loss freq")
  return(output)
}
gainRank <- function(dataset){
  gene.names = dataset[,1]
  logical = (dataset > 0)
  gain_freq = data.frame()
  for (i in 1:nrow(logical)) {
    gain_freq[i,1] = gene.names[i]
    gain_freq[i,2] = sum(logical[i,2:ncol(logical)], na.rm = TRUE)
    gain_freq[i,3] = (ncol(logical)-1) - sum(logical[i,2:ncol(logical)], na.rm = TRUE)
    gain_freq[i,4] = sum(logical[i,2:ncol(logical)], na.rm = TRUE)/(ncol(logical)-1)
  }
  ordered = order(gain_freq[,2])
  x = c()
  y = c()
  z = c()
  a = c()
  b = c()
  for (i in 1:length(ordered)) {
    x[i] = i
    y[i] = gain_freq[ordered[i], 4]
    z[i] = gene.names[ordered[i]]
    a[i] = gain_freq[ordered[i], 2]
    b[i] = gain_freq[ordered[i], 3]
  }
  output = data.frame(z, x, a, b, y)
  colnames(output) = c("Hugo_Symbol", "Gene Rank", "total copy number gain", 
                       "total unaltered", "copy number gain freq")
  return(output)
}
