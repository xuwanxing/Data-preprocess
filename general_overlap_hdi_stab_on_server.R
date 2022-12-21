library(dplyr)
library(data.table)
args <- commandArgs(TRUE)

MHT.correction <- function(df){
  ## Bonferonni
  assoc_Bonferonni <- p.adjust(df$pval, method = "bonferroni")
  df$Bonferonni_global <- assoc_Bonferonni
  ## BY
  assoc_BY <- p.adjust(df$pval, method="BY")
  df$BY_global <- assoc_BY
  ## BH 
  assoc_BH <- p.adjust(df$pval, method="BH")
  df$BH_global <- assoc_BH
  
  return(df)
}

hdi<-fread(args[1],sep="\t",header=T,check.names=F)
hdi.mht.cor<-MHT.correction(hdi)
hdi.mht.cor.sort<-hdi.mht.cor[order(hdi.mht.cor$BH_global),]

stab<-fread(args[2],sep="\t",check.names=F,header=T)
stab.filt <- stab[stab$taxa != "None",] 
colnames(stab.filt)[2]<-"gene"
colnames(stab.filt)[3]<-"taxa"

overlap_hdi_stabsel <- merge(hdi.mht.cor.sort,stab.filt, by = c("gene","taxa"))
overlap_hdi_stabsel <- overlap_hdi_stabsel[order(overlap_hdi_stabsel$BH_global),]
overlap_hdi_stabsel.BH.0.1 <- overlap_hdi_stabsel[overlap_hdi_stabsel$BH_global < 0.1,]

write.table(overlap_hdi_stabsel.BH.0.1, file=paste0(getwd(),"/",args[3],"_overlap_hdi_stab_BH.0.1.txt"), quote=F, sep="\t",col.names = NA, row.names = TRUE)
