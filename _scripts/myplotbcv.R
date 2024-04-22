myplotbcv <- function(dge, formula, rngbcv, syms, cutoffdisp1=1.5,
                      cutoffdisp2=1, cutoffA=7) {
  require(plotrix)
  plotBCV(dge, pch=".", cex=4, las=1, main=formula, ylim=rngbcv)
  A <- dge$AveLogCPM
  disp <- getDispersion(dge)
  f <- cut(dge$prior.df, breaks=c(0, 2, 4, 8, Inf))
  priordfcolors <- c("darkred", "darkorange", "darkgreen", "black")
  names(priordfcolors) <- levels(f)
  for (pd in unique(f)) {
    mask <- f == pd
    points(A[mask], sqrt(disp)[mask], pch=".", cex=4,
           col=priordfcolors[pd])
  }
  mask <- sqrt(disp) > cutoffdisp1 | (sqrt(disp) > cutoffdisp2 & A > cutoffA)
  pos <- setNames(thigmophobe(A[mask], sqrt(disp)[mask]), syms[mask])
  text(A[mask], sqrt(disp)[mask], labels=syms[mask], pos=pos)
  legend("right", names(priordfcolors), fill=priordfcolors, inset=0.01,
         title="Prior df")
}
