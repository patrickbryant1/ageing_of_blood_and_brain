library("IlluminaHumanMethylation450kanno.ilmn12.hg19")
## get the 450k annotation data
ann450k = getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)
colnames(ann450k)
ann450k$Name
write.csv(ann450k, 'anno450k.csv')
