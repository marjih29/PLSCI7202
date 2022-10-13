dat <- read.delim("A619final2.map.txt")
dat.sub <- dat[, 3:ncol(dat)]
tab <- table(as.matrixdat <- read.delim("F2_A619.map.txt")
dat.sub <- dat[, 3:ncol(dat)]
tab <- table(as.matrix(dat.sub))
write.csv(tab, file = "Table.originalA619.GC.csv")(dat.sub))
write.csv(tab, file = "Table.A619.GC.csv")

dat <- read.delim("F2_C49A.map.txt")
dat.sub <- dat[, 3:ncol(dat)]
tab <- table(as.matrix(dat.sub))
write.csv(tab, file = "Table.originalC49A.GC.csv")
