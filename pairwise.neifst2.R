pairwise.neifst2<-function (dat, diploid = TRUE) 
{
  dat <- dat[order(dat[, 1]), ]
  pops <- unique(dat[, 1])
  npop <- length(pops)
  fstmat <- matrix(nrow = npop, ncol = npop, dimnames = list(pops, 
                                                             pops))
  if (is.factor(dat[, 1])) {
    dat[, 1] <- as.numeric(dat[, 1])
    pops <- as.numeric(pops)
  }
  for (a in 2:npop) {
    for (b in 1:(a - 1)) {
      subdat <- dat[dat[, 1] == pops[a] | dat[, 1] == pops[b], 
                    ]
      fstmat[a, b] <- fstmat[b, a] <- basic.stats2(subdat, 
                                                  diploid = diploid)$overall[8]
    }
  }
  fstmat
}