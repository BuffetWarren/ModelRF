#Get command arguments - species prefix
library(protr)



args<-commandArgs(TRUE)

geneFile <- args[1]

proteinFile <- args[2]

featureFile <- args[3]

peptide_list_essential=readFASTA(proteinFile)
coding_list_essential=readFASTA(geneFile)

#peptide_list_essential <- peptide_list_essential[(sapply(peptide_list_essential, protcheck))]
peptide_list_essential <- peptide_list_essential[(sapply(peptide_list_essential, protcheck))]

#print(coding_list_essential)
print("Begin Extract feature to AA sequence")

print("####Generation of features composition from protein sequence###")
# ## Sub group AAC
a1 = t(sapply(peptide_list_essential, extractAAC))
colnames(a1) <- paste("AAC", colnames(a1), sep = "_")
featureFileSave=paste(featureFile,"Feature_collection_AAC.csv",sep="")
write.table(a1,featureFileSave,quote=FALSE,sep=",")
file.copy(from = featureFileSave, to = paste("subgroup/",featureFileSave))
#file.remove(paste("subgroup/",featureFileSave))

## Sub group Geary
a2 = t(sapply(peptide_list_essential, extractGeary, nlag=30))
colnames(a2) <- paste("Geary", colnames(a2), sep = "_")
featureFileSave=paste(featureFile,"Feature_collection_Geary.csv",sep="")
write.table(a2,featureFileSave,quote=FALSE,sep=",")
file.copy(from = featureFileSave, to = paste("subgroup/",featureFileSave))

## APAAC
# a3 = t(sapply(peptide_list_essential, extractPAAC,lambda=30))
# colnames(a3) <- paste("APAAC", colnames(a3), sep = "_")
# featureFileSave=paste(featureFile,"Feature_collection_APAAC.csv",sep="")
# write.table(a3,featureFileSave,quote=FALSE,sep=",") 
#file.copy(from = featureFileSave, to = paste("subgroup/",featureFileSave))

## PAAC
a4 = t(sapply(peptide_list_essential, extractPAAC,lambda=30))
colnames(a4) <- paste("PAAC", colnames(a4), sep = "_")
featureFileSave=paste(featureFile,"Feature_collection_PAAC.csv",sep="")
write.table(a4,featureFileSave,quote=FALSE,sep=",")
file.copy(from = featureFileSave, to = paste("subgroup/",featureFileSave))

## CTDC
a5 = t(sapply(peptide_list_essential, extractCTDC))
colnames(a5) <- paste("CTDC", colnames(a5), sep = "_")
featureFileSave=paste(featureFile,"Feature_collection_CTDC.csv",sep="")
write.table(a5,featureFileSave,quote=FALSE,sep=",")
file.copy(from = featureFileSave, to = paste("subgroup/",featureFileSave))

# ## CTDD
# # # a6 = t(sapply(peptide_list_essential, extractCTDD))
# # # colnames(a6) <- paste("CTDD", colnames(a6), sep = "_")
# # # featureFileSave=paste(featureFile,"Feature_collection_CTDD.csv",sep="")
# # # write.table(a6,featureFileSave,quote=FALSE,sep=",")
# # # file.copy(from = featureFileSave, to = paste("subgroup/",featureFileSave))
# #  
## CTDT
a7 = t(sapply(peptide_list_essential, extractCTDT))
colnames(a7) <- paste("CTDT", colnames(a7), sep = "_")
featureFileSave=paste(featureFile,"Feature_collection_CTDT.csv",sep="")
write.table(a7,featureFileSave,quote=FALSE,sep=",")
file.copy(from = featureFileSave, to = paste("subgroup/",featureFileSave))

## SOCN
a8 = t(sapply(peptide_list_essential, extractSOCN, nlag=30))
colnames(a8) <- paste("SOCN", colnames(a8), sep = "_")
featureFileSave=paste(featureFile,"Feature_collection_SOCN.csv",sep="")
write.table(a8,featureFileSave,quote=FALSE,sep=",")
file.copy(from = featureFileSave, to = paste("subgroup/",featureFileSave))

print('## End  of features composition from protein sequence##')

print('####Generation of features composition from gene sequence###"')
library("rDNAse")
coding_list_essential <- coding_list_essential[(sapply(coding_list_essential, dnacheck))]
## DACC
b1 = t(sapply(coding_list_essential, extrDACC))
colnames(b1) <- paste("DACC", colnames(b1), sep = "_")
featureFileSave=paste(featureFile,"Feature_collection_DACC.csv",sep="")
write.table(b1,featureFileSave,quote=FALSE,sep=",")
file.copy(from = featureFileSave, to = paste("subgroup_gene/",featureFileSave))
#file.remove(paste("subgroup/",featureFileSave))

## DCC
b2 = t(sapply(coding_list_essential, extrDCC))
colnames(b2) <- paste("DCC", colnames(b2), sep = "_")
featureFileSave=paste(featureFile,"Feature_collection_DCC.csv",sep="")
write.table(b2,featureFileSave,quote=FALSE,sep=",")
file.copy(from = featureFileSave, to = paste("subgroup_gene/",featureFileSave))

## KMER_2
b3 = t(sapply(coding_list_essential, kmer, 2))
colnames(b3) <- paste("kmer_2", colnames(b3), sep = "_")
featureFileSave=paste(featureFile,"Feature_collection_KMER_2.csv",sep="")
write.table(b3,featureFileSave,quote=FALSE,sep=",")
file.copy(from = featureFileSave, to = paste("subgroup_gene/",featureFileSave))

## KMER_3
b4 = t(sapply(coding_list_essential, kmer, 3))
colnames(b4) <- paste("KMER_3", colnames(b4), sep = "_")
featureFileSave=paste(featureFile,"Feature_collection_KMER_3.csv",sep="")
write.table(b4,featureFileSave,quote=FALSE,sep=",")
file.copy(from = featureFileSave, to = paste("subgroup_gene/",featureFileSave))

## PseKNC_3
b5 = t(sapply(coding_list_essential, extrPseKNC, 3))
colnames(b5) <- paste("PseKNC_3", colnames(b5), sep = "_")
featureFileSave=paste(featureFile,"Feature_collection_PseKNC_3.csv",sep="")
write.table(b5,featureFileSave,quote=FALSE,sep=",")
file.copy(from = featureFileSave, to = paste("subgroup_gene/",featureFileSave))

## PseKNC_5
b6 = t(sapply(coding_list_essential, extrPseKNC, 5))
colnames(b6) <- paste("PseKNC_5", colnames(b6), sep = "_")
featureFileSave=paste(featureFile,"Feature_collection_PseKNC_5.csv",sep="")
write.table(b6,featureFileSave,quote=FALSE,sep=",")
file.copy(from = featureFileSave, to = paste("subgroup_gene/",featureFileSave))

## PseDNC
b7 = t(sapply(coding_list_essential, extrPseDNC))
colnames(b7) <- paste("PseDNC", colnames(b7), sep = "_")
featureFileSave=paste(featureFile,"Feature_collection_PseDNC.csv",sep="")
write.table(b7,featureFileSave,quote=FALSE,sep=",")
file.copy(from = featureFileSave, to = paste("subgroup_gene/",featureFileSave))

## TACC
b8 = t(sapply(coding_list_essential, extrTACC))
colnames(b8) <- paste("TACC", colnames(b8), sep = "_")
featureFileSave=paste(featureFile,"Feature_collection_TACC.csv",sep="")
write.table(b8,featureFileSave,quote=FALSE,sep=",")
file.copy(from = featureFileSave, to = paste("subgroup_gene/",featureFileSave))

## TCC
b9 = t(sapply(coding_list_essential, extrTCC))
colnames(b9) <- paste("TCC", colnames(b9), sep = "_")
featureFileSave=paste(featureFile,"Feature_collection_TCC.csv",sep="")
write.table(b9,featureFileSave,quote=FALSE,sep=",")
file.copy(from = featureFileSave, to = paste("subgroup_gene/",featureFileSave))
print('## End  of features composition from gene sequence##'