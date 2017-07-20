basic.stats2<-function (data, diploid = TRUE, digits = 4) 
{
  if (is.genind(data)) 
    data <- genind2hierfstat(data)
  loc.names <- names(data)[-1]
  if (length(table(data[, 1])) < 2) 
    data[dim(data)[1] + 1, 1] <- data[dim(data)[1], 1] + 
    1
  if (dim(data)[2] == 2) 
    data <- data.frame(data, dummy.loc = data[, 2])
  p <- getp(data[,-1], data[,1])
  n <- t(ind.count(data))
  if (diploid) {
    dum <- data[, -1]
    sHo<-calc.Hobs.pop(dum, data[,1])
    sHo<-matrix(unlist(sHo), nrow=dim(data[,-1])[2], byrow=FALSE)
    mHo <- apply(sHo, 1, mean, na.rm = TRUE)
  }
  else {
    sHo <- NA
    mHo <- NA
  }
  sp2 <- lapply(p, fun <- function(x) apply(x, 2, fun2 <- function(x) sum(x^2)))
  sp2 <- matrix(unlist(sp2), nrow = dim(data[, -1])[2], byrow = TRUE)
  # sp2<-calc.Hexp.pop(data[,-1], data[,1])
  # sp2<-matrix(unlist(sp2), nrow = dim(data[, -1])[2], byrow = TRUE)
  if (diploid) {
    # Hs <- (1 - sp2 - sHo/2/n)
    # Hs <- n/(n - 1) * Hs
    Hs <- n/(n - 1) * (1 - sp2 - sHo/2/n)
    Fis = 1 - sHo/Hs
  }
  else {
    Hs <- n/(n - 1) * (1 - sp2)
    Fis <- NA
  }
  np <- apply(n, 1, fun <- function(x) sum(!is.na(x))) #number of population for each locus
  mn <- apply(n, 1, fun <- function(x) {
    np <- sum(!is.na(x))
    np/sum(1/x[!is.na(x)])
  })
  msp2 <- apply(sp2, 1, mean, na.rm = TRUE)
  mp <- lapply(p, fun <- function(x) apply(x, 1, mean, na.rm = TRUE))
  mp2 <- unlist(lapply(mp, fun1 <- function(x) sum(x^2)))
  if (diploid) {
    mHs <- mn/(mn - 1) * (1 - msp2 - mHo/2/mn)
    Ht <- 1 - mp2 + mHs/mn/np - mHo/2/mn/np
    mFis = 1 - mHo/mHs
  }
  else {
    mHs <- mn/(mn - 1) * (1 - msp2)
    Ht <- 1 - mp2 + mHs/mn/np
    mFis <- NA
  }
  Dst <- Ht - mHs
  Dstp <- np/(np - 1) * Dst
  Htp = mHs + Dstp
  Fst = Dst/Ht
  Fstp = Dstp/Htp
  Dest <- Dstp/(1 - mHs)
  res <- data.frame(cbind(mHo, mHs, Ht, Dst, Htp, Dstp, Fst, 
                          Fstp, mFis, Dest))
  names(res) <- c("Ho", "Hs", "Ht", "Dst", "Htp", "Dstp", "Fst", 
                  "Fstp", "Fis", "Dest")
  if (diploid) {
    rownames(sHo) <- loc.names
    rownames(Fis) <- loc.names
  }
  is.na(res) <- do.call(cbind, lapply(res, is.infinite))
  overall <- apply(res, 2, mean, na.rm = TRUE)
  overall[7] <- overall[4]/overall[3]
  overall[8] <- overall[6]/overall[5]
  overall[9] <- 1 - overall[1]/overall[2]
  overall[10] <- overall[6]/(1 - overall[2])
  names(overall) <- names(res)
  if (!diploid) {
    overall[-2] <- NA
  }
  all.res <- list(n.ind.samp = n, pop.freq = lapply(p, round, 
                                                    digits), Ho = round(sHo, digits), Hs = round(Hs, digits), 
                  Fis = round(Fis, digits), perloc = round(res, digits), 
                  overall = round(overall, digits))
  class(all.res) <- "basic.stats"
  all.res
}