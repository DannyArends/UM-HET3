#
# time_range.R
#
# copyright (c) 2020-2030 - Danny Arends
#
# At Vita loci, determine the time range at which that loci act
#

setwd("/home/rqdt9/Github/UM-HET3")
source("adjustXprobs.R")
setwd("/home/rqdt9/OneDrive/Documents/HU-Berlin/UM-HET3/files")

lods.mM <- read.table("progressiveMapping_males.txt", sep = "\t", check.names=FALSE)
lods.fM <- read.table("progressiveMapping_females.txt", sep = "\t", check.names=FALSE)
lods.cM <- read.table("progressiveMapping_all.txt", sep = "\t", check.names=FALSE)

threshold <- 3.65

rownames(lods.cM)[which(lods.cM[,"1_3010272"] > threshold)]
rownames(lods.cM)[which(lods.cM[,"1_24042124"] > threshold)]
rownames(lods.cM)[which(lods.cM[,"2_112712327"] > threshold)]
rownames(lods.cM)[which(lods.cM[,"4_55012301"] > threshold)]
rownames(lods.cM)[which(lods.cM[,"6_107382038"] > threshold)]
rownames(lods.cM)[which(lods.cM[,"9_104091597"] > threshold)]
rownames(lods.cM)[which(lods.cM[,"10_72780332"] > threshold)]
rownames(lods.cM)[which(lods.cM[,"12_112855820"] > threshold)]
rownames(lods.cM)[which(lods.cM[,"13_89689878"] > threshold)]
rownames(lods.cM)[which(lods.cM[,"14_101437457"] > threshold)]
rownames(lods.cM)[which(lods.cM[,"15_74248242"] > threshold)]
rownames(lods.cM)[which(lods.cM[,"17_32883804"] > threshold)]
rownames(lods.cM)[which(lods.cM[,"18_60822951"] > threshold)]
rownames(lods.cM)[which(lods.cM[,"X_36008085"] > threshold)]
rownames(lods.cM)[which(lods.cM[,"X_156343080"] > threshold)]


rownames(lods.fM)[which(lods.fM[,"1_24042124"] > threshold)]
rownames(lods.fM)[which(lods.fM[,"2_139956785"] > threshold)]
rownames(lods.fM)[which(lods.fM[,"3_92135706"] > threshold)]
rownames(lods.fM)[which(lods.fM[,"9_29939029"] > threshold)]
rownames(lods.fM)[which(lods.fM[,"9_104091597"] > threshold)]
rownames(lods.fM)[which(lods.fM[,"11_82178599"] > threshold)]
rownames(lods.fM)[which(lods.fM[,"X_156343080"] > threshold)]


rownames(lods.mM)[which(lods.mM[,"1_3010272"] > threshold)]
rownames(lods.mM)[which(lods.mM[,"1_120474787"] > threshold)]
rownames(lods.mM)[which(lods.mM[,"2_112712327"] > threshold)]
rownames(lods.mM)[which(lods.mM[,"3_83838529"] > threshold)]
rownames(lods.mM)[which(lods.mM[,"4_52524395"] > threshold)]
rownames(lods.mM)[which(lods.mM[,"5_67573068"] > threshold)]
rownames(lods.mM)[which(lods.mM[,"6_134870385"] > threshold)]
rownames(lods.mM)[which(lods.mM[,"9_124056586"] > threshold)]
rownames(lods.mM)[which(lods.mM[,"10_72780332"] > threshold)]
rownames(lods.mM)[which(lods.mM[,"11_5628810"] > threshold)]
rownames(lods.mM)[which(lods.mM[,"15_62405371"] > threshold)]
rownames(lods.mM)[which(lods.mM[,"17_34460077"] > threshold)]


lods.mM <- read.table("progressiveMapping_pat_males.txt", sep = "\t", check.names=FALSE)
lods.fM <- read.table("progressiveMapping_pat_females.txt", sep = "\t", check.names=FALSE)
lods.cM <- read.table("progressiveMapping_pat_all.txt", sep = "\t", check.names=FALSE)


rownames(lods.cM)[which(lods.cM[,"1_3010274"] > threshold)]
rownames(lods.cM)[which(lods.cM[,"1_132295971"] > threshold)]
rownames(lods.cM)[which(lods.cM[,"2_139956785"] > threshold)]
rownames(lods.cM)[which(lods.cM[,"9_54904313"] > threshold)]
rownames(lods.cM)[which(lods.cM[,"9_124056586"] > threshold)]
rownames(lods.cM)[which(lods.cM[,"12_112429350"] > threshold)]

rownames(lods.fM)[which(lods.fM[,"2_139956785"] > threshold)]
rownames(lods.fM)[which(lods.fM[,"2_157112564"] > threshold)]

rownames(lods.mM)[which(lods.mM[,"2_119210799"] > threshold)]
rownames(lods.mM)[which(lods.mM[,"4_145301445"] > threshold)]
rownames(lods.mM)[which(lods.mM[,"6_108075853"] > threshold)]
rownames(lods.mM)[which(lods.mM[,"9_58081975"] > threshold)]

lods.mM <- read.table("progressiveMapping_mat_males.txt", sep = "\t", check.names=FALSE)
lods.fM <- read.table("progressiveMapping_mat_females.txt", sep = "\t", check.names=FALSE)
lods.cM <- read.table("progressiveMapping_mat_all.txt", sep = "\t", check.names=FALSE)

rownames(lods.cM)[which(lods.cM[,"4_52524395"] > threshold)]
rownames(lods.cM)[which(lods.cM[,"4_74811205"] > threshold)]
rownames(lods.cM)[which(lods.cM[,"15_79013396"] > threshold)]
rownames(lods.cM)[which(lods.cM[,"17_18001459"] > threshold)]
rownames(lods.cM)[which(lods.cM[,"17_68770703"] > threshold)]
rownames(lods.cM)[which(lods.cM[,"18_52488251"] > threshold)]

rownames(lods.fM)[which(lods.fM[,"8_36994142"] > threshold)]
rownames(lods.fM)[which(lods.fM[,"13_20905668"] > threshold)]
rownames(lods.fM)[which(lods.fM[,"18_52488251"] > threshold)]
rownames(lods.fM)[which(lods.fM[,"18_67884842"] > threshold)]

rownames(lods.mM)[which(lods.mM[,"4_52524395"] > threshold)]
rownames(lods.mM)[which(lods.mM[,"6_54992703"] > threshold)]
rownames(lods.mM)[which(lods.mM[,"17_18001459"] > threshold)]


