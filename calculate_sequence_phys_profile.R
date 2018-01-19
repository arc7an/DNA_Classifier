if(!require(ape)){stop('Code requires library "ape".\n')}
if(!require(reldna)){stop('Code requires library "reldna".\n')}
if(!require(parallel)){stop('Code requires library "parallel".\n')}

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
dynchars<-function(seq, interval, tss, strand=c('forward','reverse')) {
  seq<-unlist(strsplit(seq, ''))

  if (missing(interval))
  stop("Need to specify interval")

  if(!is.character(seq))
    stop("Sequence must be a character vector containing A, C, G, T letters only")
  # if((tss - interval < 1) || (tss + interval > length(seq)))
  #   stop("Considered interval exceeds the size of the sequence")
  strand <- match.arg(strand)
  seq <- toupper(seq)
  seq <- switch(strand, 
              forward=seq,
              reverse=rev(seq))
  
  seq <- seq[(tss - 150):(tss + 250 - 1)]
  
  a<-3.4*10^(-10)
  #I<-c(7.6, 4.8, 8.2, 4.1)*10^(-44)
  K<-c(227, 155, 220, 149)*10^(-20)
  V<-c(2.09, 1.43, 3.12, 2.12)*10^(-20)

  csA<-cumsum(seq=='A')
  csT<-cumsum(seq=='T')
  csG<-cumsum(seq=='G')
  csC<-cumsum(seq=='C')

  countA = csA[interval:length(csA)]-c(0, csA[1:(length(csA)-interval)])
  countT = csT[interval:length(csT)]-c(0, csT[1:(length(csT)-interval)])
  countG = csG[interval:length(csG)]-c(0, csG[1:(length(csG)-interval)])
  countC = csC[interval:length(csC)]-c(0, csC[1:(length(csC)-interval)])

  
  M <- switch(strand,
              forward = cbind(countA, countT, countG, countC)/interval,
              reverse = cbind(countT, countA, countC, countG)/interval)

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
#' @param extended_string
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
  p <- lseqspline1D(
    substr(extended_string, tss_position-ep_interval[1], tss_position+ep_interval[2]),
    bound=c(50, 350),
    ref=251)
  return(
    switch(strand,
           forward=p$mpot[p$x %in% zout],
           reverse=rev(p$mpot[p$x %in% zout])))
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
calculate_profile <- function(dnaSeq, dynamic_interval, tss_position, ep_interval=c(267, 217), zout=-480:239, strand=c('forward','reverse')) {
  #ep_interval = c(163, 213) for reverse, for forward c(267, 217), zout - the same
  dynRes <- dynchars(dnaSeq, dynamic_interval, tss_position, strand)
  epRes <- calculate_EP_on_interval(tss_position, dnaSeq, ep_interval, zout, strand)
  return(c(dynRes$E0, dynRes$d, epRes, dynRes$gc))
}

