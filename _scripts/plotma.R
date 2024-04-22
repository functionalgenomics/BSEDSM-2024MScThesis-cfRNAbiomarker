plotma <- function(tt, FCcutoff=1.5, FDRcutoff=0.1, highlight=character(0),
                   labelTop=10) {
  labelTop <- min(length(highlight), labelTop)
  par(mar=c(4, 5, 1, 1))
  xrng <- c(floor(min(tt$AveExpr)), ceiling(max(tt$AveExpr)))
  yrng <- c(floor(min(tt$logFC)), ceiling(max(tt$logFC)))
  plot(tt$AveExpr, tt$logFC, pch=".", cex=3, col=grey(0.75),
       main="", xlab=expression(paste(log[2], " CPM Average Expression")), axes=FALSE,
       xlim=xrng, ylim=yrng,
       ylab=expression(paste(log[2], " fold change")), las=1, cex.lab=1.4)
  axis(1, at=seq(xrng[1], xrng[2], by=1), labels=TRUE, cex.axis=1.2)
  axis(2, at=seq(yrng[1], yrng[2], by=1), labels=TRUE, cex.axis=1.2, las=1)
  abline(h=0, col="gray", lwd=2, lty=1)
  abline(h=c(-log2(FCcutoff), log2(FCcutoff)), col="gray", lwd=2, lty=2)
  points(tt$AveExpr, tt$logFC, pch=".", cex=3, col=gray(0.75))
  mask <- rownames(tt) %in% highlight
  if (any(mask)) {
    points(tt$AveExpr[mask], tt$logFC[mask], pch=".", cex=3, col="darkred")
    idxtopgenesByAveExpr <- order(tt$AveExpr, decreasing=TRUE)
    idxtopgenesByAveExpr <- idxtopgenesByAveExpr[idxtopgenesByAveExpr %in% which(mask)]
    idxtopgenesByFC <- order(abs(tt$logFC), decreasing=TRUE)
    idxtopgenesByFC <- idxtopgenesByFC[idxtopgenesByFC %in% which(mask)]
    idxtopgenes <- unique(c(idxtopgenesByAveExpr[1:labelTop], idxtopgenesByFC[1:labelTop]))
    if (length(idxtopgenes) > 2)
      posg <- thigmophobe(tt[idxtopgenes, "AveExpr"], tt[idxtopgenes, "logFC"])
    else
      posg <- rep(3, length(idxtopgenes))
    names(posg) <- tt[idxtopgenes, "Symbol"]
    text(tt$AveExpr[idxtopgenes], tt$logFC[idxtopgenes], labels=names(posg),
         pos=posg, col="darkred", cex=0.5, font=2)
  }
  text(max(tt$AveExpr)-0.7, -log2(FCcutoff), sprintf("> %d%%\nfold change", (FCcutoff-1)*100),
       pos=1, cex=0.85)
  text(max(tt$AveExpr)-0.7, log2(FCcutoff), sprintf("> %d%%\nfold change", (FCcutoff-1)*100),
       pos=3, cex=0.85)
}
