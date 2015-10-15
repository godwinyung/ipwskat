SNPInfo   <- read.table("out.pos-1",header=T)
R <- 10000
regions.start  <- sample(min(SNPInfo$CHROM_POS):(max(SNPInfo$CHROM_POS)-3e3), R, replace=TRUE)
regions.pos    <- matrix(c(regions.start, regions.start+3e3), ncol=2)
regions        <- matrix(0, nrow=R, ncol=2)
for (r in 1:R) {
  temp <- which(SNPInfo$CHROM_POS >= regions.pos[r,1] & SNPInfo$CHROM_POS <= regions.pos[r,2])
  regions[r,] <- c(min(temp), max(temp))
}
write.table(regions,"regions3e3.txt",row.names=F,col.names=F)

