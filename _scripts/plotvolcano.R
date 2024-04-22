plotvolcano <- function(tt, FCcutoff=1.5, FDRcutoff=0.1,
                        highlight=character(0), labelTop=10) {
  require(plotrix)
  labelTop <- min(length(highlight), labelTop)
  par(mar=c(4, 5, 1, 1))
  xrng <- c(floor(min(tt$logFC)), ceiling(max(tt$logFC)))
  yrng <- c(floor(min(-log10(tt$P.Value))), ceiling(max(-log10(tt$P.Value)))-0.5)
  plot(tt$logFC, -log10(tt$P.Value), pch=".", cex=3, col=grey(0.75),
       main="", xlab=expression(paste(log[2], " fold change")), axes=FALSE,
       xlim=xrng, ylim=yrng,
       ylab=expression(paste(-log[10], " raw P-value")), las=1, cex.lab=1.4)
  axis(1, at=seq(xrng[1], xrng[2], by=1), labels=TRUE, cex.axis=1.2)
  axis(2, at=seq(yrng[1], yrng[2], by=1), labels=TRUE, cex.axis=1.2, las=1)
  abline(v=c(-log2(FCcutoff), log2(FCcutoff)), col="gray", lwd=2, lty=2)
  abline(h=-log10(max(tt$P.Value[tt$adj.P.Val <= FDRcutoff])),
         col="gray", lwd=2, lty=2)
  points(tt$logFC, -log10(tt$P.Value), pch=".", cex=3, col=gray(0.75))
  mask <- rownames(tt) %in% highlight
  if (any(mask)) {
    points(tt$logFC[mask], -log10(tt$P.Value[mask]), pch=".", cex=3, col="darkred")
    idxtopgenesByP <- order(tt$P.Value, decreasing=FALSE)
    idxtopgenesByP <- idxtopgenesByP[idxtopgenesByP %in% which(mask)]
    idxtopgenesByFC <- order(abs(tt$logFC), decreasing=TRUE)
    idxtopgenesByFC <- idxtopgenesByFC[idxtopgenesByFC %in% which(mask)]
    idxtopgenes <- unique(c(idxtopgenesByP[1:labelTop], idxtopgenesByFC[1:labelTop]))
    if (length(idxtopgenes) > 2)
      posg <- thigmophobe(tt[idxtopgenes, "logFC"], -log10(tt[idxtopgenes, "P.Value"]))
    else
      posg <- rep(3, length(idxtopgenes))
    names(posg) <- tt[idxtopgenes, "Symbol"]
    text(tt$logFC[idxtopgenes], -log10(tt$P.Value[idxtopgenes]), labels=names(posg),
         pos=posg, col="darkred", cex=0.5, font=2)
  }
  text(max(tt$logFC)*0.95,
       -log10(max(tt$P.Value[tt$adj.P.Val <= FDRcutoff])),
       sprintf("%d%% FDR", 100*FDRcutoff), pos=1, cex=0.85)
  text(-log2(FCcutoff), 0.05, sprintf("> %d%%\nfold change", (FCcutoff-1)*100),
       pos=2, cex=0.85)
  text(log2(FCcutoff), 0.05, sprintf("> %d%%\nfold change", (FCcutoff-1)*100),
       pos=4, cex=0.85)
}
