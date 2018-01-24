if(!require(ape)){stop('Code requires library "ape".\n')}
if(!require(reldna)){stop('Code requires library "reldna".\n')}
if(!require(parallel)){stop('Code requires library "parallel".\n')}
if(!require(Biostrings)){stop('Code requires library "Biostrings".\n')}

get_forward_subseq <- function(seq, tss, boundaries) {
  return(as.character(
    DNAString(seq,
              start = tss-boundaries[1],
              nchar = sum(boundaries)+1)))
}

get_reverse_subseq <- function(seq, tss, boundaries) {
  # return(as.character(
  #   reverseComplement(
  #     DNAString(seq,
  #               start = tss-boundaries[2],
  #               nchar = sum(boundaries)+1))))
  #DNAString(rcshort, start = 501 - 10 - 1, nchar = sum(c(10,15))+1)
  return (as.character(
    DNAString(seq,
              start = tss - boundaries[1] - 1,
              nchar = sum(boundaries)+1)))
}
get_forward_substring <- function(seq, tss, boundaries) {
   return(seq[(tss-boundaries[1]):(tss+boundaries[2])])
}

get_reverse_substring <- function(seq, tss, boundaries) {
  res <- rev(chartr("ATGC", "TACG", seq[(tss-boundaries[2]):(tss+boundaries[1])]))
#  return(tail(res, n=1))
  return(res)
}

#' Title
#'
#' @param seq full sequence (chromosome)
#' @param interval requested interval
#' @param tss
#' @param strand DNA strand
#'
#' @return
#' @export
#'
#' @examples
dynchars<-function(seq, average_interval_size, tss, boundaries, strand=c('forward','reverse')) {

  if (missing(average_interval_size))
  stop("Need to specify average_interval_size")

  if(!is.character(seq))
    stop("Sequence must be a character vector containing A, C, G, T letters only")
  # if((tss - average_interval_size < 1) || (tss + average_interval_size > length(seq)))
  #   stop("Considered average_interval_size exceeds the size of the sequence")
  seq <- switch(strand, 
                forward = seq,
                reverse = as.character(reverseComplement(DNAString(seq))))
  seq <- unlist(strsplit(seq, ''))
  seq <- toupper(seq)

  boundaries <- boundaries + average_interval_size %/% 2

  strand <- match.arg(strand)
  seq <- switch(strand,
              forward=get_forward_substring(seq, tss, boundaries),
              reverse=get_reverse_substring(seq, tss, boundaries))

  a<-3.4*10^(-10)
  #I<-c(7.6, 4.8, 8.2, 4.1)*10^(-44)
  K<-c(227, 155, 220, 149)*10^(-20)
  V<-c(2.09, 1.43, 3.12, 2.12)*10^(-20)

  csA<-cumsum(seq=='A')
  csT<-cumsum(seq=='T')
  csG<-cumsum(seq=='G')
  csC<-cumsum(seq=='C')

  countA = csA[average_interval_size:length(csA)]-c(0, csA[1:(length(csA)-average_interval_size)])
  countT = csT[average_interval_size:length(csT)]-c(0, csT[1:(length(csT)-average_interval_size)])
  countG = csG[average_interval_size:length(csG)]-c(0, csG[1:(length(csG)-average_interval_size)])
  countC = csC[average_interval_size:length(csC)]-c(0, csC[1:(length(csC)-average_interval_size)])

  M <- cbind(countA, countT, countG, countC)/average_interval_size

  #Is <- as.numeric(M%*%I)#! numeric conversion
  Ks <- as.numeric(M%*%K)
  Vs <- as.numeric(M%*%V)

  E0 <- (8*(Ks*Vs)^0.5)* 6E23 / 4184
  d <- ((Ks*a^2)/Vs)^(0.5)/a;
  gc = M[,3] + M[,4]

  dynchars_return<-list(E0=E0, d=d, gc=gc)
  return(dynchars_return)
}

#' Title
#'
#' @param tss_position
#' @param extended_string -- vector of characters
#' @param ep_interval
#' @param zout
#' @param strand
#'
#' @return
#' @export
#'
#' @examples
calculate_EP_on_interval <- function(tss_position, extended_string, ep_interval=c(250, 150), zout=-480:239, strand=c('forward','reverse')) {
  #lseqspline1D(substr(e.coli_U00096.2, exp_tsss[i]-250, exp_tsss[i]+150), bound=c(50, 350), ref=251 )
  strand<- match.arg(strand)
  #extended_string <- unlist(strsplit(extended_string, ''))
  subseq <-switch(strand,
                  forward = get_forward_substring(extended_string, tss_position, ep_interval),
                  reverse = get_reverse_substring(extended_string, tss_position, ep_interval))
  subseq <- paste0(subseq, collapse = "")
  p <- lseqspline1D(
    subseq,
    bound=c(50, 350),
    ref=251)
  return(p$mpot[p$x %in% zout])
}

#' Title
#'
#' @param dnaSeq
#' @param dynamic_interval
#' @param tss_position
#' @param ep_interval
#' @param zout
#' @param strand
#'
#' @return
#' @export
#'
#' @examples
calculate_profile <- function(dnaSeq, average_interval_size, tss_position, dynamic_interval, ep_interval=c(267, 217), zout=-480:239, strand=c('forward','reverse')) {
  #ep_interval = c(163, 213) for reverse, for forward c(267, 217), zout - the same
  #seq, average_interval_size, tss, boundaries, strand=c('forward','reverse')
  #dnaSeq <- unlist(strsplit(dnaSeq, ''))
  dynRes <- dynchars(dnaSeq, average_interval_size, tss_position,  dynamic_interval, strand)
  epRes <- calculate_EP_on_interval(tss_position, dnaSeq, ep_interval, zout, strand)
  return(c(dynRes$E0, dynRes$d, epRes, dynRes$gc))
}

