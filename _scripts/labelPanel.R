## function to add a lable to a panel in a multipanel figure.
## parameters: lab - label, typically a letter "a", "b", "c", etc.
##             font - font type, font=2 is boldface, see 'font' in 'help(par)'
##             cex - character expansion factor; see 'cex' in 'help(text)'
##             offsetx - offset factor in the x-axis to adjust label position
##             offsety - offset factor in the y-axis to adjust label position
labelPanel <- function(lab, font=2, cex=2, offsetx=0.05, offsety=0.05) {
  par(xpd=TRUE)
  w <- par("usr")[2] - par("usr")[1]
  h <- par("usr")[4] - par("usr")[3]
  text(par("usr")[1]-w*offsetx, par("usr")[4]+h*offsety, lab, font=font, cex=cex)
  par(xpd=FALSE)
}
