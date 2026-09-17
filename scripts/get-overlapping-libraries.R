# Get all library pairs with overlapping donors
library(data.table)
library(foreach)
library(doMC)
registerDoMC(cores=4)

args <- commandArgs(trailingOnly=TRUE)
infile <- args[1]
outfile <- args[2]

dat <- fread(infile)
dat[, donor := tstrsplit(specimenID, split='_')[1]]


libraries <- dat$libraryBatch

combos <- unique(CJ(LIB1=libraries, LIB2=libraries))[LIB1 != LIB2]

o <- foreach(i=1:nrow(combos), .combine='rbind') %dopar% {
    l1 <- combos[i,LIB1]
    l2 <- combos[i,LIB2]
    tmp <- dat[libraryBatch %in% c(l1,l2)][, .N, by=donor][N>1]
    if(nrow(tmp) > 0) { return(data.table(L1=min(c(l1,l2)), L2=max(c(l1,l2)))) }
}

o <- unique(o)

fwrite(o, file=outfile, quote=F, row.names=F, col.names=F, sep='\t')