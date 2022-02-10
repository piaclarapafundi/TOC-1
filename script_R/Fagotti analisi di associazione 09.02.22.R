setwd("D:\\Progetti GSTeP\\Fagotti-Perrone TOC-1 - ongoing")
dati<-read_excel("TOC1.xlsx")
colnames(dati)

install.packages("survminer") #diagnostica su hazard proporzionali coeflo di Cox
install.packages("coxphw") #coeflo di Cox su hazard non proporzionali
install.packages("Hmisc") #matrice di correlazione (di Spearman)
install.packages("RVAideMemoire") #intervallo di confidenza su correlazione (di Sperman)
install.packages("spearmanCI") #intervallo di confidenza su correlazione (di Sperman)
install.packages("DescTools") #intervallo di confidenza su correlazione (di Sperman), mediante approssimazione Z
install.packages("My.stepwise")
install.packages("glmnet")
install.packages("reshape") 

library(reshape)
library(survival)
library(survminer)
library(coxphw)
library(Hmisc)
library(RVAideMemoire)
library(spearmanCI)
library(DescTools)
library(My.stepwise)
library(MASS)
library(glmnet)

#################################
##########VALID DATASET##########
#################################

dati<-data.frame(ID=dati$"Paziente, cod sanit", status_alive_dead=dati$"Status (0=alive, 1=DOD)", relapse=dati$"Relapse or progression", OS_diagn=dati$"OS dalla diagnosi", PFSultimoPLATINO=dati$"PFS dall'ultimo platino", DFS=dati$"DFS dalla citoriduzione",BRCA_dic=dati$"BRCA1\\2_mut=1\\ WT=0",BRCA=dati$"BRCA wt_1_2", age=dati$Età,ARnucleo=dati$ARnucleo_Score_num,PR=dati$PR_Score_num,ERalpha=dati$ERAlpha_Score_num,ERbeta1nucleo=dati$ERBeta1_NucleoScore_num,ERbeta1cito=dati$ERBeta1_CitoScore_num,ERbeta2nucleo=dati$ERBeta2_NucleoScore_num,ERbeta2cito=dati$ERBeta2_CitoScore_num,ERbeta5nucleo=dati$ERBeta5_NucleoScore_num,ERbeta5cito=dati$ERBeta5_CitoScore_num, a_b1n=dati$"alpha_beta1n", a_b2n=dati$"alpha_beta2n", a_b5n=dati$"alpha_beta5n",  somatico=dati$"somatico codificato", PLARES=dati$"Platinum-resistance (0:no; 1:si)",BMI=dati$bmi_at_diagnosis,menopausa=dati$menopausa_at_diagnosis,ascite=dati$chir1_ascite,chir_pi=dati$chir1_pi_lps_score,ct_Beva=dati$ct_beva,citorid_cm=dati$"Citoridotte =1 , Non citoridotte = 0",citorid_cm=dati$"Cytored_RT(in cm)")

dim(dati)
summary(dati)
BRCA=dati$BRCA; length(BRCA)
BRCA_dic=dati$BRCA_dic

View(dati)




hist(dati[,10])
hist(dati[,11])
hist(dati[,12])
hist(dati[,13])
hist(dati[,14])
hist(dati[,15])
hist(dati[,16])
hist(dati[,17])
hist(dati[,18])



apply(dati[,9:18],2,function(x) shapiro.test(x))

colnames(dati)


####################################
##ANALISI ASSOCIAZIONE OVERALL######
####################################


###OVERALL###

rcorr(as.matrix(dati[,10:18]),type="spearman")
rcorr(as.matrix(dati[,10:18]),type="pearson")

library(reshape2)
rdati<-data.frame(rcorr(as.matrix(dati[,10:18]),type="spearman")$r); rdati
melted_cormat <- melt(rdati); melted_cormat
melted_cormat2<-data.frame(rep(colnames(rdati),9),melted_cormat)
colnames(melted_cormat2)<-c("Var2","Var1","value")
head(melted_cormat2)
melted_cormat2[,1]<-factor(melted_cormat2[,1],levels=colnames(rdati))
head(melted_cormat2)

library(ggplot2)

#ggplot(data = melted_cormat2, aes(x=Var2, y=Var1, fill=value)) + geom_tile()

ggplot(data = melted_cormat2, aes(x=Var2, y=Var1, fill=value)) + 
 geom_tile(color = "white")+
 scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
   midpoint = 0, limit = c(-1,1), space = "Lab", 
   name="Spearman\nCorrelation") +
  theme_minimal()+ 
 theme(axis.text.x = element_text(angle = 45, vjust = 1, 
    size = 12, hjust = 1))+
 coord_fixed()

spearman.ci(dati[,10],dati[,11])
spearman.ci(dati[,10],dati[,12])
spearman.ci(dati[,10],dati[,13])
spearman.ci(dati[,10],dati[,14])
spearman.ci(dati[,10],dati[,15])
spearman.ci(dati[,10],dati[,16])
spearman.ci(dati[,10],dati[,17])
spearman.ci(dati[,10],dati[,18])

spearman.ci(dati[,11],dati[,12])
spearman.ci(dati[,11],dati[,13])
spearman.ci(dati[,11],dati[,14])
spearman.ci(dati[,11],dati[,15])
spearman.ci(dati[,11],dati[,16])
spearman.ci(dati[,11],dati[,17])
spearman.ci(dati[,11],dati[,18])

spearman.ci(dati[,12],dati[,13])
spearman.ci(dati[,12],dati[,14])
spearman.ci(dati[,12],dati[,15])
spearman.ci(dati[,12],dati[,16])
spearman.ci(dati[,12],dati[,17])
spearman.ci(dati[,12],dati[,18])

spearman.ci(dati[,13],dati[,14])
spearman.ci(dati[,13],dati[,15])
spearman.ci(dati[,13],dati[,16])
spearman.ci(dati[,13],dati[,17])
spearman.ci(dati[,13],dati[,18])

spearman.ci(dati[,14],dati[,15])
spearman.ci(dati[,14],dati[,16])
spearman.ci(dati[,14],dati[,17])
spearman.ci(dati[,14],dati[,18])

spearman.ci(dati[,15],dati[,16])
spearman.ci(dati[,15],dati[,17])
spearman.ci(dati[,15],dati[,18])

spearman.ci(dati[,16],dati[,17])
spearman.ci(dati[,16],dati[,18])

spearman.ci(dati[,17],dati[,18])


##################################
###ANALISI DI ASSOCIAZIONE BRCA###
##################################

###BRCA_0###

rcorr(as.matrix(dati[dati$BRCA==0,10:18]),type="spearman")
rcorr(as.matrix(dati[dati$BRCA==0,10:18]),type="pearson")

library(reshape2)
rdati<-data.frame(rcorr(as.matrix(dati[dati$BRCA==0,10:18]),type="spearman")$r); rdati
melted_cormat <- melt(rdati); melted_cormat
melted_cormat2<-data.frame(rep(colnames(rdati),9),melted_cormat)
colnames(melted_cormat2)<-c("Var2","Var1","value")
head(melted_cormat2)
melted_cormat2[,1]<-factor(melted_cormat2[,1],levels=colnames(rdati))
head(melted_cormat2)

library(ggplot2)

#ggplot(data = melted_cormat2, aes(x=Var2, y=Var1, fill=value)) + geom_tile()

ggplot(data = melted_cormat2, aes(x=Var2, y=Var1, fill=value)) + 
 geom_tile(color = "white")+
 scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
   midpoint = 0, limit = c(-1,1), space = "Lab", 
   name="Spearman\nCorrelation within BRCA 0") +
  theme_minimal()+ 
 theme(axis.text.x = element_text(angle = 45, vjust = 1, 
    size = 12, hjust = 1))+
 coord_fixed()


spearman.ci(dati[,10][dati$BRCA==0],dati[,11][dati$BRCA==0])
spearman.ci(dati[,10][dati$BRCA==0],dati[,12][dati$BRCA==0])
spearman.ci(dati[,10][dati$BRCA==0],dati[,13][dati$BRCA==0])
spearman.ci(dati[,10][dati$BRCA==0],dati[,14][dati$BRCA==0])
spearman.ci(dati[,10][dati$BRCA==0],dati[,15][dati$BRCA==0])
spearman.ci(dati[,10][dati$BRCA==0],dati[,16][dati$BRCA==0])
spearman.ci(dati[,10][dati$BRCA==0],dati[,17][dati$BRCA==0])
spearman.ci(dati[,10][dati$BRCA==0],dati[,18][dati$BRCA==0])

spearman.ci(dati[,11][dati$BRCA==0],dati[,12][dati$BRCA==0])
spearman.ci(dati[,11][dati$BRCA==0],dati[,13][dati$BRCA==0])
spearman.ci(dati[,11][dati$BRCA==0],dati[,14][dati$BRCA==0])
spearman.ci(dati[,11][dati$BRCA==0],dati[,15][dati$BRCA==0])
spearman.ci(dati[,11][dati$BRCA==0],dati[,16][dati$BRCA==0])
spearman.ci(dati[,11][dati$BRCA==0],dati[,17][dati$BRCA==0])
spearman.ci(dati[,11][dati$BRCA==0],dati[,18][dati$BRCA==0])

spearman.ci(dati[,12][dati$BRCA==0],dati[,13][dati$BRCA==0])
spearman.ci(dati[,12][dati$BRCA==0],dati[,14][dati$BRCA==0])
spearman.ci(dati[,12][dati$BRCA==0],dati[,15][dati$BRCA==0])
spearman.ci(dati[,12][dati$BRCA==0],dati[,16][dati$BRCA==0])
spearman.ci(dati[,12][dati$BRCA==0],dati[,17][dati$BRCA==0])
spearman.ci(dati[,12][dati$BRCA==0],dati[,18][dati$BRCA==0])

spearman.ci(dati[,13][dati$BRCA==0],dati[,14][dati$BRCA==0])
spearman.ci(dati[,13][dati$BRCA==0],dati[,15][dati$BRCA==0])
spearman.ci(dati[,13][dati$BRCA==0],dati[,16][dati$BRCA==0])
spearman.ci(dati[,13][dati$BRCA==0],dati[,17][dati$BRCA==0])
spearman.ci(dati[,13][dati$BRCA==0],dati[,18][dati$BRCA==0])

spearman.ci(dati[,14][dati$BRCA==0],dati[,15][dati$BRCA==0])
spearman.ci(dati[,14][dati$BRCA==0],dati[,16][dati$BRCA==0])
spearman.ci(dati[,14][dati$BRCA==0],dati[,17][dati$BRCA==0])
SpearmanRho(dati[,14][dati$BRCA==0],dati[,17][dati$BRCA==0],conf.level=0.95) 
#con DescTools
spearman.ci(dati[,14][dati$BRCA==0],dati[,18][dati$BRCA==0])

spearman.ci(dati[,15][dati$BRCA==0],dati[,16][dati$BRCA==0])
spearman.ci(dati[,15][dati$BRCA==0],dati[,17][dati$BRCA==0])
spearman.ci(dati[,15][dati$BRCA==0],dati[,18][dati$BRCA==0])

spearman.ci(dati[,16][dati$BRCA==0],dati[,17][dati$BRCA==0])
spearman.ci(dati[,16][dati$BRCA==0],dati[,18][dati$BRCA==0])

spearman.ci(dati[,17][dati$BRCA==0],dati[,18][dati$BRCA==0])


###BRCA 1###

rcorr(as.matrix(dati[dati$BRCA==1,10:18]),type="spearman")
rcorr(as.matrix(dati[dati$BRCA==1,10:18]),type="pearson")


library(reshape2)
rdati<-data.frame(rcorr(as.matrix(dati[dati$BRCA==1,10:18]),type="spearman")$r); rdati
melted_cormat <- melt(rdati); melted_cormat
melted_cormat2<-data.frame(rep(colnames(rdati),9),melted_cormat)
colnames(melted_cormat2)<-c("Var2","Var1","value")
head(melted_cormat2)
melted_cormat2[,1]<-factor(melted_cormat2[,1],levels=colnames(rdati))
head(melted_cormat2)

library(ggplot2)

#ggplot(data = melted_cormat2, aes(x=Var2, y=Var1, fill=value)) + geom_tile()

ggplot(data = melted_cormat2, aes(x=Var2, y=Var1, fill=value)) + 
 geom_tile(color = "white")+
 scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
   midpoint = 0, limit = c(-1,1), space = "Lab", 
   name="Spearman\nCorrelation within BRCA 1") +
  theme_minimal()+ 
 theme(axis.text.x = element_text(angle = 45, vjust = 1, 
    size = 12, hjust = 1))+
 coord_fixed()


spearman.ci(dati[,10][dati$BRCA==1],dati[,11][dati$BRCA==1])
spearman.ci(dati[,10][dati$BRCA==1],dati[,12][dati$BRCA==1])
spearman.ci(dati[,10][dati$BRCA==1],dati[,13][dati$BRCA==1])
spearman.ci(dati[,10][dati$BRCA==1],dati[,14][dati$BRCA==1])
spearman.ci(dati[,10][dati$BRCA==1],dati[,15][dati$BRCA==1])
spearman.ci(dati[,10][dati$BRCA==1],dati[,16][dati$BRCA==1])
spearman.ci(dati[,10][dati$BRCA==1],dati[,17][dati$BRCA==1])
spearman.ci(dati[,10][dati$BRCA==1],dati[,18][dati$BRCA==1])

spearman.ci(dati[,11][dati$BRCA==1],dati[,12][dati$BRCA==1])
spearman.ci(dati[,11][dati$BRCA==1],dati[,13][dati$BRCA==1])
spearman.ci(dati[,11][dati$BRCA==1],dati[,14][dati$BRCA==1])
spearman.ci(dati[,11][dati$BRCA==1],dati[,15][dati$BRCA==1])
spearman.ci(dati[,11][dati$BRCA==1],dati[,16][dati$BRCA==1])
spearman.ci(dati[,11][dati$BRCA==1],dati[,17][dati$BRCA==1])
spearman.ci(dati[,11][dati$BRCA==1],dati[,18][dati$BRCA==1])

spearman.ci(dati[,12][dati$BRCA==1],dati[,13][dati$BRCA==1])
spearman.ci(dati[,12][dati$BRCA==1],dati[,14][dati$BRCA==1])
spearman.ci(dati[,12][dati$BRCA==1],dati[,15][dati$BRCA==1])
spearman.ci(dati[,12][dati$BRCA==1],dati[,16][dati$BRCA==1])
spearman.ci(dati[,12][dati$BRCA==1],dati[,17][dati$BRCA==1])
spearman.ci(dati[,12][dati$BRCA==1],dati[,18][dati$BRCA==1])

spearman.ci(dati[,13][dati$BRCA==1],dati[,14][dati$BRCA==1])
spearman.ci(dati[,13][dati$BRCA==1],dati[,15][dati$BRCA==1])
spearman.ci(dati[,13][dati$BRCA==1],dati[,16][dati$BRCA==1])
spearman.ci(dati[,13][dati$BRCA==1],dati[,17][dati$BRCA==1])
spearman.ci(dati[,13][dati$BRCA==1],dati[,18][dati$BRCA==1])

spearman.ci(dati[,14][dati$BRCA==1],dati[,15][dati$BRCA==1])
spearman.ci(dati[,14][dati$BRCA==1],dati[,16][dati$BRCA==1])
spearman.ci(dati[,14][dati$BRCA==1],dati[,17][dati$BRCA==1])
SpearmanRho(dati[,14][dati$BRCA==1],dati[,17][dati$BRCA==1],conf.level=0.95) 
#con DescTools
spearman.ci(dati[,14][dati$BRCA==1],dati[,18][dati$BRCA==1])

spearman.ci(dati[,15][dati$BRCA==1],dati[,16][dati$BRCA==1])
spearman.ci(dati[,15][dati$BRCA==1],dati[,17][dati$BRCA==1])
spearman.ci(dati[,15][dati$BRCA==1],dati[,18][dati$BRCA==1])

spearman.ci(dati[,16][dati$BRCA==1],dati[,17][dati$BRCA==1])
spearman.ci(dati[,16][dati$BRCA==1],dati[,18][dati$BRCA==1])

spearman.ci(dati[,17][dati$BRCA==1],dati[,18][dati$BRCA==1])


###BRCA_2###

rcorr(as.matrix(dati[dati$BRCA==2,10:18]),type="spearman")
rcorr(as.matrix(dati[dati$BRCA==2,10:18]),type="pearson")

library(reshape2)
rdati<-data.frame(rcorr(as.matrix(dati[dati$BRCA==2,10:18]),type="spearman")$r); rdati
melted_cormat <- melt(rdati); melted_cormat
melted_cormat2<-data.frame(rep(colnames(rdati),9),melted_cormat)
colnames(melted_cormat2)<-c("Var2","Var1","value")
head(melted_cormat2)
melted_cormat2[,1]<-factor(melted_cormat2[,1],levels=colnames(rdati))
head(melted_cormat2)

library(ggplot2)

#ggplot(data = melted_cormat2, aes(x=Var2, y=Var1, fill=value)) + geom_tile()

ggplot(data = melted_cormat2, aes(x=Var2, y=Var1, fill=value)) + 
 geom_tile(color = "white")+
 scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
   midpoint = 0, limit = c(-1,1), space = "Lab", 
   name="Spearman\nCorrelation within BRCA 2") +
  theme_minimal()+ 
 theme(axis.text.x = element_text(angle = 45, vjust = 1, 
    size = 12, hjust = 1))+
 coord_fixed()


spearman.ci(dati[,10][dati$BRCA==2],dati[,11][dati$BRCA==2])
spearman.ci(dati[,10][dati$BRCA==2],dati[,12][dati$BRCA==2])
spearman.ci(dati[,10][dati$BRCA==2],dati[,13][dati$BRCA==2])
spearman.ci(dati[,10][dati$BRCA==2],dati[,14][dati$BRCA==2])
spearman.ci(dati[,10][dati$BRCA==2],dati[,15][dati$BRCA==2])
spearman.ci(dati[,10][dati$BRCA==2],dati[,16][dati$BRCA==2])
spearman.ci(dati[,10][dati$BRCA==2],dati[,17][dati$BRCA==2])
spearman.ci(dati[,10][dati$BRCA==2],dati[,18][dati$BRCA==2])

spearman.ci(dati[,11][dati$BRCA==2],dati[,12][dati$BRCA==2])
spearman.ci(dati[,11][dati$BRCA==2],dati[,13][dati$BRCA==2])
spearman.ci(dati[,11][dati$BRCA==2],dati[,14][dati$BRCA==2])
spearman.ci(dati[,11][dati$BRCA==2],dati[,15][dati$BRCA==2])
spearman.ci(dati[,11][dati$BRCA==2],dati[,16][dati$BRCA==2])
spearman.ci(dati[,11][dati$BRCA==2],dati[,17][dati$BRCA==2])
spearman.ci(dati[,11][dati$BRCA==2],dati[,18][dati$BRCA==2])

spearman.ci(dati[,12][dati$BRCA==2],dati[,13][dati$BRCA==2])
spearman.ci(dati[,12][dati$BRCA==2],dati[,14][dati$BRCA==2])
spearman.ci(dati[,12][dati$BRCA==2],dati[,15][dati$BRCA==2])
spearman.ci(dati[,12][dati$BRCA==2],dati[,16][dati$BRCA==2])
spearman.ci(dati[,12][dati$BRCA==2],dati[,17][dati$BRCA==2])
spearman.ci(dati[,12][dati$BRCA==2],dati[,18][dati$BRCA==2])

spearman.ci(dati[,13][dati$BRCA==2],dati[,14][dati$BRCA==2])
spearman.ci(dati[,13][dati$BRCA==2],dati[,15][dati$BRCA==2])
spearman.ci(dati[,13][dati$BRCA==2],dati[,16][dati$BRCA==2])
spearman.ci(dati[,13][dati$BRCA==2],dati[,17][dati$BRCA==2])
spearman.ci(dati[,13][dati$BRCA==2],dati[,18][dati$BRCA==2])

spearman.ci(dati[,14][dati$BRCA==2],dati[,15][dati$BRCA==2])
spearman.ci(dati[,14][dati$BRCA==2],dati[,16][dati$BRCA==2])
spearman.ci(dati[,14][dati$BRCA==2],dati[,17][dati$BRCA==2])
SpearmanRho(dati[,14][dati$BRCA==2],dati[,17][dati$BRCA==2],conf.level=0.95) 
#con DescTools
spearman.ci(dati[,14][dati$BRCA==2],dati[,18][dati$BRCA==2])

spearman.ci(dati[,15][dati$BRCA==2],dati[,16][dati$BRCA==2])
spearman.ci(dati[,15][dati$BRCA==2],dati[,17][dati$BRCA==2])
spearman.ci(dati[,15][dati$BRCA==2],dati[,18][dati$BRCA==2])

spearman.ci(dati[,16][dati$BRCA==2],dati[,17][dati$BRCA==2])
spearman.ci(dati[,16][dati$BRCA==2],dati[,18][dati$BRCA==2])

spearman.ci(dati[,17][dati$BRCA==2],dati[,18][dati$BRCA==2])






######################################
###ANALISI DI ASSOCIAZIONE BRCA_dic###
######################################

###BRCA_dic 0###

rcorr(as.matrix(dati[dati$BRCA_dic==0,10:18]),type="spearman")
rcorr(as.matrix(dati[dati$BRCA_dic==0,10:18]),type="pearson")

library(reshape2)
rdati<-data.frame(rcorr(as.matrix(dati[dati$BRCA_dic==0,10:18]),type="spearman")$r); rdati
melted_cormat <- melt(rdati); melted_cormat
melted_cormat2<-data.frame(rep(colnames(rdati),9),melted_cormat)
colnames(melted_cormat2)<-c("Var2","Var1","value")
head(melted_cormat2)
melted_cormat2[,1]<-factor(melted_cormat2[,1],levels=colnames(rdati))
head(melted_cormat2)

library(ggplot2)

#ggplot(data = melted_cormat2, aes(x=Var2, y=Var1, fill=value)) + geom_tile()

ggplot(data = melted_cormat2, aes(x=Var2, y=Var1, fill=value)) + 
 geom_tile(color = "white")+
 scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
   midpoint = 0, limit = c(-1,1), space = "Lab", 
   name="Spearman\nCorrelation within BRCA_dic 0") +
  theme_minimal()+ 
 theme(axis.text.x = element_text(angle = 45, vjust = 1, 
    size = 12, hjust = 1))+
 coord_fixed()


spearman.ci(dati[,10][dati$BRCA_dic==0],dati[,11][dati$BRCA_dic==0])
spearman.ci(dati[,10][dati$BRCA_dic==0],dati[,12][dati$BRCA_dic==0])
spearman.ci(dati[,10][dati$BRCA_dic==0],dati[,13][dati$BRCA_dic==0])
spearman.ci(dati[,10][dati$BRCA_dic==0],dati[,14][dati$BRCA_dic==0])
spearman.ci(dati[,10][dati$BRCA_dic==0],dati[,15][dati$BRCA_dic==0])
spearman.ci(dati[,10][dati$BRCA_dic==0],dati[,16][dati$BRCA_dic==0])
spearman.ci(dati[,10][dati$BRCA_dic==0],dati[,17][dati$BRCA_dic==0])
spearman.ci(dati[,10][dati$BRCA_dic==0],dati[,18][dati$BRCA_dic==0])

spearman.ci(dati[,11][dati$BRCA_dic==0],dati[,12][dati$BRCA_dic==0])
spearman.ci(dati[,11][dati$BRCA_dic==0],dati[,13][dati$BRCA_dic==0])
spearman.ci(dati[,11][dati$BRCA_dic==0],dati[,14][dati$BRCA_dic==0])
spearman.ci(dati[,11][dati$BRCA_dic==0],dati[,15][dati$BRCA_dic==0])
spearman.ci(dati[,11][dati$BRCA_dic==0],dati[,16][dati$BRCA_dic==0])
spearman.ci(dati[,11][dati$BRCA_dic==0],dati[,17][dati$BRCA_dic==0])
spearman.ci(dati[,11][dati$BRCA_dic==0],dati[,18][dati$BRCA_dic==0])

spearman.ci(dati[,12][dati$BRCA_dic==0],dati[,13][dati$BRCA_dic==0])
spearman.ci(dati[,12][dati$BRCA_dic==0],dati[,14][dati$BRCA_dic==0])
spearman.ci(dati[,12][dati$BRCA_dic==0],dati[,15][dati$BRCA_dic==0])
spearman.ci(dati[,12][dati$BRCA_dic==0],dati[,16][dati$BRCA_dic==0])
spearman.ci(dati[,12][dati$BRCA_dic==0],dati[,17][dati$BRCA_dic==0])
spearman.ci(dati[,12][dati$BRCA_dic==0],dati[,18][dati$BRCA_dic==0])

spearman.ci(dati[,13][dati$BRCA_dic==0],dati[,14][dati$BRCA_dic==0])
spearman.ci(dati[,13][dati$BRCA_dic==0],dati[,15][dati$BRCA_dic==0])
spearman.ci(dati[,13][dati$BRCA_dic==0],dati[,16][dati$BRCA_dic==0])
spearman.ci(dati[,13][dati$BRCA_dic==0],dati[,17][dati$BRCA_dic==0])
spearman.ci(dati[,13][dati$BRCA_dic==0],dati[,18][dati$BRCA_dic==0])

spearman.ci(dati[,14][dati$BRCA_dic==0],dati[,15][dati$BRCA_dic==0])
spearman.ci(dati[,14][dati$BRCA_dic==0],dati[,16][dati$BRCA_dic==0])
spearman.ci(dati[,14][dati$BRCA_dic==0],dati[,17][dati$BRCA_dic==0])
SpearmanRho(dati[,14][dati$BRCA_dic==0],dati[,17][dati$BRCA_dic==0],conf.level=0.95) 
#con DescTools
spearman.ci(dati[,14][dati$BRCA_dic==0],dati[,18][dati$BRCA_dic==0])

spearman.ci(dati[,15][dati$BRCA_dic==0],dati[,16][dati$BRCA_dic==0])
spearman.ci(dati[,15][dati$BRCA_dic==0],dati[,17][dati$BRCA_dic==0])
spearman.ci(dati[,15][dati$BRCA_dic==0],dati[,18][dati$BRCA_dic==0])

spearman.ci(dati[,16][dati$BRCA_dic==0],dati[,17][dati$BRCA_dic==0])
spearman.ci(dati[,16][dati$BRCA_dic==0],dati[,18][dati$BRCA_dic==0])

spearman.ci(dati[,17][dati$BRCA_dic==0],dati[,18][dati$BRCA_dic==0])


###BRCA_dic 1###

rcorr(as.matrix(dati[dati$BRCA_dic==1,10:18]),type="spearman")
rcorr(as.matrix(dati[dati$BRCA_dic==1,10:18]),type="pearson")


library(reshape2)
rdati<-data.frame(rcorr(as.matrix(dati[dati$BRCA_dic==1,10:18]),type="spearman")$r); rdati
melted_cormat <- melt(rdati); melted_cormat
melted_cormat2<-data.frame(rep(colnames(rdati),9),melted_cormat)
colnames(melted_cormat2)<-c("Var2","Var1","value")
head(melted_cormat2)
melted_cormat2[,1]<-factor(melted_cormat2[,1],levels=colnames(rdati))
head(melted_cormat2)

library(ggplot2)

#ggplot(data = melted_cormat2, aes(x=Var2, y=Var1, fill=value)) + geom_tile()

ggplot(data = melted_cormat2, aes(x=Var2, y=Var1, fill=value)) + 
 geom_tile(color = "white")+
 scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
   midpoint = 0, limit = c(-1,1), space = "Lab", 
   name="Spearman\nCorrelation within BRCA_dic 1") +
  theme_minimal()+ 
 theme(axis.text.x = element_text(angle = 45, vjust = 1, 
    size = 12, hjust = 1))+
 coord_fixed()


spearman.ci(dati[,10][dati$BRCA_dic==1],dati[,11][dati$BRCA_dic==1])
spearman.ci(dati[,10][dati$BRCA_dic==1],dati[,12][dati$BRCA_dic==1])
spearman.ci(dati[,10][dati$BRCA_dic==1],dati[,13][dati$BRCA_dic==1])
spearman.ci(dati[,10][dati$BRCA_dic==1],dati[,14][dati$BRCA_dic==1])
spearman.ci(dati[,10][dati$BRCA_dic==1],dati[,15][dati$BRCA_dic==1])
spearman.ci(dati[,10][dati$BRCA_dic==1],dati[,16][dati$BRCA_dic==1])
spearman.ci(dati[,10][dati$BRCA_dic==1],dati[,17][dati$BRCA_dic==1])
spearman.ci(dati[,10][dati$BRCA_dic==1],dati[,18][dati$BRCA_dic==1])

spearman.ci(dati[,11][dati$BRCA_dic==1],dati[,12][dati$BRCA_dic==1])
spearman.ci(dati[,11][dati$BRCA_dic==1],dati[,13][dati$BRCA_dic==1])
spearman.ci(dati[,11][dati$BRCA_dic==1],dati[,14][dati$BRCA_dic==1])
spearman.ci(dati[,11][dati$BRCA_dic==1],dati[,15][dati$BRCA_dic==1])
spearman.ci(dati[,11][dati$BRCA_dic==1],dati[,16][dati$BRCA_dic==1])
spearman.ci(dati[,11][dati$BRCA_dic==1],dati[,17][dati$BRCA_dic==1])
spearman.ci(dati[,11][dati$BRCA_dic==1],dati[,18][dati$BRCA_dic==1])

spearman.ci(dati[,12][dati$BRCA_dic==1],dati[,13][dati$BRCA_dic==1])
spearman.ci(dati[,12][dati$BRCA_dic==1],dati[,14][dati$BRCA_dic==1])
spearman.ci(dati[,12][dati$BRCA_dic==1],dati[,15][dati$BRCA_dic==1])
spearman.ci(dati[,12][dati$BRCA_dic==1],dati[,16][dati$BRCA_dic==1])
spearman.ci(dati[,12][dati$BRCA_dic==1],dati[,17][dati$BRCA_dic==1])
spearman.ci(dati[,12][dati$BRCA_dic==1],dati[,18][dati$BRCA_dic==1])

spearman.ci(dati[,13][dati$BRCA_dic==1],dati[,14][dati$BRCA_dic==1])
spearman.ci(dati[,13][dati$BRCA_dic==1],dati[,15][dati$BRCA_dic==1])
spearman.ci(dati[,13][dati$BRCA_dic==1],dati[,16][dati$BRCA_dic==1])
spearman.ci(dati[,13][dati$BRCA_dic==1],dati[,17][dati$BRCA_dic==1])
spearman.ci(dati[,13][dati$BRCA_dic==1],dati[,18][dati$BRCA_dic==1])

spearman.ci(dati[,14][dati$BRCA_dic==1],dati[,15][dati$BRCA_dic==1])
spearman.ci(dati[,14][dati$BRCA_dic==1],dati[,16][dati$BRCA_dic==1])
spearman.ci(dati[,14][dati$BRCA_dic==1],dati[,17][dati$BRCA_dic==1])
SpearmanRho(dati[,14][dati$BRCA_dic==1],dati[,17][dati$BRCA_dic==1],conf.level=0.95) 
#con DescTools
spearman.ci(dati[,14][dati$BRCA_dic==1],dati[,18][dati$BRCA_dic==1])

spearman.ci(dati[,15][dati$BRCA_dic==1],dati[,16][dati$BRCA_dic==1])
spearman.ci(dati[,15][dati$BRCA_dic==1],dati[,17][dati$BRCA_dic==1])
spearman.ci(dati[,15][dati$BRCA_dic==1],dati[,18][dati$BRCA_dic==1])

spearman.ci(dati[,16][dati$BRCA_dic==1],dati[,17][dati$BRCA_dic==1])
spearman.ci(dati[,16][dati$BRCA_dic==1],dati[,18][dati$BRCA_dic==1])

spearman.ci(dati[,17][dati$BRCA_dic==1],dati[,18][dati$BRCA_dic==1])



