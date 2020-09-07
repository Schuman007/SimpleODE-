f <- seq(0,1,length=100)
R0 <- -log(1-f)/f
R0[is.nan(R0)] <- 0

plot(f~R0,type='l',xlab=expression(R[0]),ylab="fraction infected",bty='l')

FinalSize <- cbind.data.frame(f, R0)
