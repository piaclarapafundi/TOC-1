setwd("D:\\Progetti GSTeP\\Fagotti-Perrone TOC-1")
dati<-read_excel("TOC1.1.xlsx")
dim(dati)
colnames(dati)
head(dati)
str(dati)
summary(dati)

#install.packages("survminer") #diagnostica su hazard proporzionali modello di Cox
#install.packages("coxphw") #modello di Cox su hazard non proporzionali
#install.packages("Hmisc") #matrice di correlazione (di Spearman)
#install.packages("RVAideMemoire") #intervallo di confidenza su correlazione (di Sperman)
#install.packages("spearmanCI") #intervallo di confidenza su correlazione (di Sperman)
#install.packages("DescTools") #intervallo di confidenza su correlazione (di Sperman), mediante approssimazione Z
#install.packages("My.stepwise")
#install.packages("glmnet")

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

citation(package = "survival", lib.loc = NULL) #Cox regression
citation(package = "survminer", lib.loc = NULL) #Diagnostics on Cox regression
citation(package = "coxphw", lib.loc = NULL) #Diagnostics on Cox regression
citation(package = "Hmisc", lib.loc = NULL) #Marker Association Analysis
citation(package = "RVAideMemoire", lib.loc = NULL) #Marker Association Analysis



###################
#####OUTCOMES######
###################
colnames(dati)


#outcome primario: 
#-OS con follow-up OS dalla diagnosi


#outcome secondari: 
#-relapse or progression con follow-up DFS dalla citoriduzione




###############################################
###FATTORI PROGNOSITICI -  MARKER MOLECOLARI###
###############################################

summary(dati$ARnucleo_Score_num)
summary(dati$PR_Score_num)
summary(dati$ERAlpha_Score_num)
summary(dati$ERBeta1_NucleoScore_num)
summary(dati$ERBeta1_CitoScore_num)
summary(dati$ERBeta2_NucleoScore_num)
summary(dati$ERBeta2_CitoScore_num)
summary(dati$ERBeta5_NucleoScore_num)
summary(dati$ERBeta5_CitoScore_num)
summary(dati$"PFS dall'ultimo platino")


#################################
##########VALID DATASET##########
#################################

dati<-data.frame(ID=dati$"Paziente, cod sanit", status_alive_dead=dati$"Status (0=alive, 1=DOD)", relapse=dati$"Relapse or progression", OS_diagn=dati$"OS dalla diagnosi", PFSultimoPLATINO=dati$"PFS dall'ultimo platino", DFS=dati$"DFS dalla citoriduzione",BRCA_dic=dati$"BRCA1\\2_mut=1\\ WT=0",BRCA=dati$"BRCA wt_1_2", age=dati$Età,ARnucleo=dati$ARnucleo_Score_num,PR=dati$PR_Score_num,ERalpha=dati$ERAlpha_Score_num,ERbeta1nucleo=dati$ERBeta1_NucleoScore_num,ERbeta1cito=dati$ERBeta1_CitoScore_num,ERbeta2nucleo=dati$ERBeta2_NucleoScore_num,ERbeta2cito=dati$ERBeta2_CitoScore_num,ERbeta5nucleo=dati$ERBeta5_NucleoScore_num,ERbeta5cito=dati$ERBeta5_CitoScore_num, a_b1n=dati$"alpha_beta1n", a_b2n=dati$"alpha_beta2n", a_b5n=dati$"alpha_beta5n",  somatico=dati$"somatico codificato", PLARES=dati$"Platinum-resistance (0:no; 1:si)",BMI=dati$bmi_at_diagnosis,menopausa=dati$menopausa_at_diagnosis,ascite=dati$chir1_ascite,chir_pi=dati$chir1_pi_lps_score,ct_Beva=dati$ct_beva,citorid_cm=dati$"Citoridotte =1 , Non citoridotte = 0",citorid_cm_num=dati$"Cytored_RT(in cm)")

BRCA=dati$BRCA; length(BRCA)
BRCA_dic=dati$BRCA_dic

dim(dati)
colnames(dati)
head(dati)
summary(dati)


#####################################
###MODELLO DI COX UNIVARIATO BRCA####
#####################################

######
#BRCA#
######

#MORTALITA'
addmargins(table(dati$status_alive_dead,BRCA))
by(dati$OS_diagn,BRCA,summary)
m0=coxph(Surv(dati$OS_diagn,dati$status_alive_dead) ~ factor(BRCA), data=dati)
summary(m0)
test.ph0 <- cox.zph(m0); test.ph0
library(survminer)
ggcoxzph(test.ph0)

m0<- coxphw(Surv(OS_diagn,status_alive_dead) ~ factor(BRCA), data = dati, template = "AHR")
summary(m0)

#RECIDIVA DFS
addmargins(table(dati$relapse,BRCA))
by(dati$DFS,BRCA,summary)
m1=coxph(Surv(dati$DFS,dati$relapse) ~ factor(BRCA), data=dati)
summary(m1)
test.ph1 <- cox.zph(m1); test.ph1
library(survminer)
ggcoxzph(test.ph1)

dati$timeDFS<-ifelse(dati$DFS==0,1,dati$DFS)
m1<-coxphw(Surv(timeDFS,relapse) ~ factor(BRCA), data = dati, template = "AHR")
summary(m1)


#################
#AR nucleo score#
#################

#MORTALITA'
by(dati$ARnucleo,dati$status_alive_dead,summary)
m2=coxph(Surv(dati$OS_diagn, dati$status_alive_dead) ~ dati$ARnucleo, data=dati)
summary(m2)
test.ph2 <- cox.zph(m2); test.ph2
library(survminer)
ggcoxzph(test.ph2)

m2=coxph(Surv(dati$OS_diagn, dati$status_alive_dead) ~ dati$ARnucleo*factor(BRCA), data=dati)
summary(m2)
test.ph2 <- cox.zph(m2); test.ph2
library(survminer)
ggcoxzph(test.ph2)

m2=coxph(Surv(dati$OS_diagn, dati$status_alive_dead) ~ dati$ARnucleo*factor(BRCA)+dati$age, data=dati)
summary(m2)
test.ph2 <- cox.zph(m2); test.ph2
library(survminer)
ggcoxzph(test.ph2)

m2 <- coxphw(Surv(OS_diagn, status_alive_dead) ~ ARnucleo, data=dati, template = "AHR")
summary(m2)

m2 <- coxphw(Surv(dati$OS_diagn, dati$status_alive_dead) ~ dati$ARnucleo*factor(BRCA), data = dati, template = "AHR")
summary(m2)


#RECIDIVA DFS
by(dati$ARnucleo,dati$relapse,summary)
m1=coxph(Surv(dati$DFS,dati$relapse) ~ dati$ARnucleo, data=dati)
summary(m1)
test.ph1 <- cox.zph(m1); test.ph1
library(survminer)
ggcoxzph(test.ph1)

m1=coxph(Surv(dati$DFS,dati$relapse) ~ dati$ARnucleo*factor(BRCA), data=dati)
summary(m1)
test.ph1 <- cox.zph(m1); test.ph1
library(survminer)
ggcoxzph(test.ph1)

m1=coxph(Surv(dati$DFS,dati$relapse) ~ dati$ARnucleo*factor(BRCA)+dati$age, data=dati)
summary(m1)
test.ph1 <- cox.zph(m1); test.ph1
library(survminer)
ggcoxzph(test.ph1)

dati$timeDFS<-ifelse(dati$DFS==0,1,dati$DFS)
m1<-coxphw(Surv(timeDFS,dati$relapse) ~ dati$ARnucleo, data = dati, template = "AHR")
summary(m1)

m1<-coxphw(Surv(timeDFS,dati$relapse) ~ dati$ARnucleo*factor(BRCA), data = dati, template = "AHR")
summary(m1)



##############
#PR_Score_num#
##############


#MORTALITA'
by(dati$PR,dati$status_alive_dead,summary)
m2=coxph(Surv(dati$OS_diagn, dati$status_alive_dead) ~ dati$PR, data=dati)
summary(m2)
test.ph2 <- cox.zph(m2); test.ph2
library(survminer)
ggcoxzph(test.ph2)

m2=coxph(Surv(dati$OS_diagn, dati$status_alive_dead) ~ dati$PR*factor(BRCA), data=dati)
summary(m2)
test.ph2 <- cox.zph(m2); test.ph2
library(survminer)
ggcoxzph(test.ph2)

m2=coxph(Surv(dati$OS_diagn, dati$status_alive_dead) ~ dati$PR*factor(BRCA)+dati$age, data=dati)
summary(m2)
test.ph2 <- cox.zph(m2); test.ph2
library(survminer)
ggcoxzph(test.ph2)

m2 <- coxphw(Surv(OS_diagn, status_alive_dead) ~ PR,data = dati, template = "AHR")
summary(m2)

m2 <- coxphw(Surv(OS_diagn, status_alive_dead) ~ PR*factor(BRCA),data = dati, template = "AHR")
summary(m2)


#RECIDIVA DFS
by(dati$PR,dati$relapse,summary)
m1=coxph(Surv(dati$DFS,dati$relapse) ~ dati$PR, data=dati)
summary(m1)
test.ph1 <- cox.zph(m1); test.ph1
library(survminer)
ggcoxzph(test.ph1)

m1=coxph(Surv(dati$DFS,dati$relapse) ~ dati$PR*factor(BRCA), data=dati)
summary(m1)
test.ph1 <- cox.zph(m1); test.ph1
library(survminer)
ggcoxzph(test.ph1)

m1=coxph(Surv(dati$DFS,dati$relapse) ~ dati$PR*factor(BRCA)+dati$age, data=dati)
summary(m1)
test.ph1 <- cox.zph(m1); test.ph1
library(survminer)
ggcoxzph(test.ph1)

dati$timeDFS<-ifelse(dati$DFS==0,1,dati$DFS)
m1<-coxphw(Surv(timeDFS,relapse) ~ PR, data = dati, template = "AHR")
summary(m1)

m1<-coxphw(Surv(timeDFS,relapse) ~ PR*factor(BRCA), data = dati, template = "AHR")
summary(m1)



###################
#ERAlpha_Score_num#
###################

#MORTALITA'
by(dati$ERalpha,dati$status_alive_dead,summary)
m4=coxph(Surv(dati$OS_diagn, dati$status_alive_dead) ~ dati$ERalpha, data=dati)
summary(m4)
test.ph4 <- cox.zph(m4); test.ph4
library(survminer)
ggcoxzph(test.ph4)

m4=coxph(Surv(dati$OS_diagn, dati$status_alive_dead) ~ dati$ERalpha*factor(BRCA), data=dati)
summary(m4)
test.ph4 <- cox.zph(m4); test.ph4
library(survminer)
ggcoxzph(test.ph4)

m4=coxph(Surv(dati$OS_diagn, dati$status_alive_dead) ~ dati$ERalpha*factor(BRCA)+dati$age, data=dati)
summary(m4)
test.ph4 <- cox.zph(m4); test.ph4
library(survminer)
ggcoxzph(test.ph4)

m4 <- coxphw(Surv(OS_diagn, status_alive_dead) ~ ERalpha,data = dati, template = "AHR")
summary(m4)

m4 <- coxphw(Surv(OS_diagn, status_alive_dead) ~ ERalpha*factor(BRCA),data = dati, template = "AHR")
summary(m4)


#RECIDIVA DFS
by(dati$ERalpha,dati$relapse,summary)
m4=coxph(Surv(dati$DFS,dati$relapse) ~ dati$ERalpha, data=dati)
summary(m4)
test.ph4 <- cox.zph(m4); test.ph4
library(survminer)
ggcoxzph(test.ph4)

m4=coxph(Surv(dati$DFS,dati$relapse) ~ dati$ERalpha*factor(BRCA), data=dati)
summary(m4)
test.ph4 <- cox.zph(m4); test.ph4
library(survminer)
ggcoxzph(test.ph4)

m4=coxph(Surv(dati$DFS,dati$relapse) ~ dati$ERalpha*factor(BRCA)+dati$age, data=dati)
summary(m4)
test.ph4 <- cox.zph(m4); test.ph4
library(survminer)
ggcoxzph(test.ph4)

dati$timeDFS<-ifelse(dati$DFS==0,1,dati$DFS)
m4<-coxphw(Surv(timeDFS,relapse) ~ ERalpha, data = dati, template = "AHR")
summary(m4)

m4<-coxphw(Surv(timeDFS,relapse) ~ ERalpha*factor(BRCA), data = dati, template = "AHR")
summary(m4)



#########################
#ERBeta1_NucleoScore_num#
#########################

#MORTALITA'
by(dati$ERbeta1nucleo,dati$status_alive_dead,summary)
m2=coxph(Surv(dati$OS_diagn, dati$status_alive_dead) ~ dati$ERbeta1nucleo, data=dati)
summary(m2)
test.ph2 <- cox.zph(m2); test.ph2
library(survminer)
ggcoxzph(test.ph2)

m2=coxph(Surv(dati$OS_diagn, dati$status_alive_dead) ~ dati$ERbeta1nucleo*factor(BRCA), data=dati)
summary(m2)
test.ph2 <- cox.zph(m2); test.ph2
library(survminer)
ggcoxzph(test.ph2)

m2=coxph(Surv(dati$OS_diagn, dati$status_alive_dead) ~ dati$ERbeta1nucleo*factor(BRCA)+dati$age, data=dati)
summary(m2)
test.ph2 <- cox.zph(m2); test.ph2
library(survminer)
ggcoxzph(test.ph2)

m2 <- coxphw(Surv(OS_diagn, status_alive_dead) ~ ERbeta1nucleo,data = dati, template = "AHR")
summary(m2)

m2 <- coxphw(Surv(OS_diagn, status_alive_dead) ~ ERbeta1nucleo*factor(BRCA),data = dati, template = "AHR")
summary(m2)


#RECIDIVA DFS
by(dati$ERbeta1nucleo,dati$relapse,summary)
m1=coxph(Surv(dati$DFS,dati$relapse) ~ dati$ERbeta1nucleo, data=dati)
summary(m1)
test.ph1 <- cox.zph(m1); test.ph1
library(survminer)
ggcoxzph(test.ph1)

m1=coxph(Surv(dati$DFS,dati$relapse) ~ dati$ERbeta1nucleo*factor(BRCA), data=dati)
summary(m1)
test.ph1 <- cox.zph(m1); test.ph1
library(survminer)
ggcoxzph(test.ph1)

m1=coxph(Surv(dati$DFS,dati$relapse) ~ dati$ERbeta1nucleo*factor(BRCA)+dati$age, data=dati)
summary(m1)
test.ph1 <- cox.zph(m1); test.ph1
library(survminer)
ggcoxzph(test.ph1)

dati$timeDFS<-ifelse(dati$DFS==0,1,dati$DFS)
m1<-coxphw(Surv(timeDFS,relapse) ~ ERbeta1nucleo, data = dati, template = "AHR")
summary(m1)

m1<-coxphw(Surv(timeDFS,relapse) ~ ERbeta1nucleo*factor(BRCA), data = dati, template = "AHR")
summary(m1)



#######################
#ERBeta1_CitoScore_num#
#######################

#MORTALITA'
by(dati$ERbeta1cito,dati$status_alive_dead,summary)
m2=coxph(Surv(dati$OS_diagn, dati$status_alive_dead) ~ dati$ERbeta1cito, data=dati)
summary(m2)
test.ph2 <- cox.zph(m2); test.ph2
library(survminer)
ggcoxzph(test.ph2)

m2=coxph(Surv(dati$OS_diagn, dati$status_alive_dead) ~ dati$ERbeta1cito*factor(BRCA), data=dati)
summary(m2)
test.ph2 <- cox.zph(m2); test.ph2
library(survminer)
ggcoxzph(test.ph2)

m2=coxph(Surv(dati$OS_diagn, dati$status_alive_dead) ~ dati$ERbeta1cito*factor(BRCA)+dati$age, data=dati)
summary(m2)
test.ph2 <- cox.zph(m2); test.ph2
library(survminer)
ggcoxzph(test.ph2)

m2 <- coxphw(Surv(OS_diagn, status_alive_dead) ~ ERbeta1cito,data = dati, template = "AHR")
summary(m2)

m2 <- coxphw(Surv(OS_diagn, status_alive_dead) ~ ERbeta1cito*factor(BRCA),data = dati, template = "AHR")
summary(m2)


#RECIDIVA DFS
by(dati$ERbeta1cito,dati$relapse,summary)
m1=coxph(Surv(dati$DFS,dati$relapse) ~ dati$ERbeta1cito, data=dati)
summary(m1)
test.ph1 <- cox.zph(m1); test.ph1
library(survminer)
ggcoxzph(test.ph1)

m1=coxph(Surv(dati$DFS,dati$relapse) ~ dati$ERbeta1cito*factor(BRCA), data=dati)
summary(m1)
test.ph1 <- cox.zph(m1); test.ph1
library(survminer)
ggcoxzph(test.ph1)

m1=coxph(Surv(dati$DFS,dati$relapse) ~ dati$ERbeta1cito*factor(BRCA)+dati$age, data=dati)
summary(m1)
test.ph1 <- cox.zph(m1); test.ph1
library(survminer)
ggcoxzph(test.ph1)

dati$timeDFS<-ifelse(dati$DFS==0,1,dati$DFS)
m1<-coxphw(Surv(timeDFS,relapse) ~ ERbeta1cito, data = dati, template = "AHR")
summary(m1)

m1<-coxphw(Surv(timeDFS,relapse) ~ ERbeta1cito*factor(BRCA), data = dati, template = "AHR")
summary(m1)



#########################
#ERBeta2_NucleoScore_num#
#########################

#MORTALITA'
by(dati$ERbeta2nucleo,dati$status_alive_dead,summary)
m2=coxph(Surv(dati$OS_diagn, dati$status_alive_dead) ~ dati$ERbeta2nucleo, data=dati)
summary(m2)
test.ph2 <- cox.zph(m2); test.ph2
library(survminer)
ggcoxzph(test.ph2)

m2=coxph(Surv(dati$OS_diagn, dati$status_alive_dead) ~ dati$ERbeta2nucleo*factor(BRCA), data=dati)
summary(m2)
test.ph2 <- cox.zph(m2); test.ph2
library(survminer)
ggcoxzph(test.ph2)

m2=coxph(Surv(dati$OS_diagn, dati$status_alive_dead) ~ dati$ERbeta2nucleo*factor(BRCA)+dati$age, data=dati)
summary(m2)
test.ph2 <- cox.zph(m2); test.ph2
library(survminer)
ggcoxzph(test.ph2)

m2 <- coxphw(Surv(OS_diagn, status_alive_dead) ~ ERbeta2nucleo,data = dati, template = "AHR")
summary(m2)

m2 <- coxphw(Surv(OS_diagn, status_alive_dead) ~ ERbeta2nucleo*factor(BRCA),data = dati, template = "AHR")
summary(m2)


#RECIDIVA DFS
by(dati$ERbeta2nucleo,dati$relapse,summary)
m1=coxph(Surv(dati$DFS,dati$relapse) ~ dati$ERbeta2nucleo, data=dati)
summary(m1)
test.ph1 <- cox.zph(m1); test.ph1
library(survminer)
ggcoxzph(test.ph1)

m1=coxph(Surv(dati$DFS,dati$relapse) ~ dati$ERbeta2nucleo*factor(BRCA), data=dati)
summary(m1)
test.ph1 <- cox.zph(m1); test.ph1
library(survminer)
ggcoxzph(test.ph1)

m1=coxph(Surv(dati$DFS,dati$relapse) ~ dati$ERbeta2nucleo*factor(BRCA)+dati$age, data=dati)
summary(m1)
test.ph1 <- cox.zph(m1); test.ph1
library(survminer)
ggcoxzph(test.ph1)

dati$timeDFS<-ifelse(dati$DFS==0,1,dati$DFS)
m1<-coxphw(Surv(timeDFS,relapse) ~ ERbeta2nucleo, data = dati, template = "AHR")
summary(m1)
m1<-coxphw(Surv(timeDFS,relapse) ~ ERbeta2nucleo*factor(BRCA), data = dati, template = "AHR")
summary(m1)


#######################
#ERBeta2_CitoScore_num#
#######################


#MORTALITA'
by(dati$ERbeta2cito,dati$status_alive_dead,summary)
m2=coxph(Surv(dati$OS_diagn, dati$status_alive_dead) ~ dati$ERbeta2cito, data=dati)
summary(m2)
test.ph2 <- cox.zph(m2); test.ph2
library(survminer)
ggcoxzph(test.ph2)

m2=coxph(Surv(dati$OS_diagn, dati$status_alive_dead) ~ dati$ERbeta2cito*factor(BRCA), data=dati)
summary(m2)
test.ph2 <- cox.zph(m2); test.ph2
library(survminer)
ggcoxzph(test.ph2)

m2=coxph(Surv(dati$OS_diagn, dati$status_alive_dead) ~ dati$ERbeta2cito*factor(BRCA)+dati$age, data=dati)
summary(m2)
test.ph2 <- cox.zph(m2); test.ph2
library(survminer)
ggcoxzph(test.ph2)

m2 <- coxphw(Surv(OS_diagn, status_alive_dead) ~ ERbeta2cito,data = dati, template = "AHR")
summary(m2)

m2 <- coxphw(Surv(OS_diagn, status_alive_dead) ~ ERbeta2cito*factor(BRCA),data = dati, template = "AHR")
summary(m2)


#RECIDIVA DFS
by(dati$ERbeta2cito,dati$relapse,summary)
m1=coxph(Surv(dati$DFS,dati$relapse) ~ dati$ERbeta2cito, data=dati)
summary(m1)
test.ph1 <- cox.zph(m1); test.ph1
library(survminer)
ggcoxzph(test.ph1)

m1=coxph(Surv(dati$DFS,dati$relapse) ~ dati$ERbeta2cito*factor(BRCA), data=dati)
summary(m1)
test.ph1 <- cox.zph(m1); test.ph1
library(survminer)
ggcoxzph(test.ph1)

m1=coxph(Surv(dati$DFS,dati$relapse) ~ dati$ERbeta2cito*factor(BRCA)+dati$age, data=dati)
summary(m1)
test.ph1 <- cox.zph(m1); test.ph1
library(survminer)
ggcoxzph(test.ph1)

dati$timeDFS<-ifelse(dati$DFS==0,1,dati$DFS)
m1<-coxphw(Surv(timeDFS,relapse) ~ ERbeta2cito, data = dati, template = "AHR")
summary(m1)

m1<-coxphw(Surv(timeDFS,relapse) ~ ERbeta2cito*factor(BRCA), data = dati, template = "AHR")
summary(m1)


#########################
#ERBeta5_NucleoScore_num#
#########################

#MORTALITA'
by(dati$ERbeta5nucleo,dati$status_alive_dead,summary)
m2=coxph(Surv(dati$OS_diagn, dati$status_alive_dead) ~ dati$ERbeta5nucleo, data=dati)
summary(m2)
test.ph2 <- cox.zph(m2); test.ph2
library(survminer)
ggcoxzph(test.ph2)

m2=coxph(Surv(dati$OS_diagn, dati$status_alive_dead) ~ dati$ERbeta5nucleo*factor(BRCA), data=dati)
summary(m2)
test.ph2 <- cox.zph(m2); test.ph2
library(survminer)
ggcoxzph(test.ph2)

m2=coxph(Surv(dati$OS_diagn, dati$status_alive_dead) ~ dati$ERbeta5nucleo*factor(BRCA)+dati$age, data=dati)
summary(m2)
test.ph2 <- cox.zph(m2); test.ph2
library(survminer)
ggcoxzph(test.ph2)

m2 <- coxphw(Surv(OS_diagn, status_alive_dead) ~ ERbeta5nucleo,data = dati, template = "AHR")
summary(m2)

m2 <- coxphw(Surv(OS_diagn, status_alive_dead) ~ ERbeta5nucleo*factor(BRCA),data = dati, template = "AHR")
summary(m2)


#RECIDIVA DFS
by(dati$ERbeta5nucleo,dati$relapse,summary)
m1=coxph(Surv(dati$DFS,dati$relapse) ~ dati$ERbeta5nucleo, data=dati)
summary(m1)
test.ph1 <- cox.zph(m1); test.ph1
library(survminer)
ggcoxzph(test.ph1)

m1=coxph(Surv(dati$DFS,dati$relapse) ~ dati$ERbeta5nucleo*factor(BRCA), data=dati)
summary(m1)
test.ph1 <- cox.zph(m1); test.ph1
library(survminer)
ggcoxzph(test.ph1)

m1=coxph(Surv(dati$DFS,dati$relapse) ~ dati$ERbeta5nucleo*factor(BRCA)+dati$age, data=dati)
summary(m1)
test.ph1 <- cox.zph(m1); test.ph1
library(survminer)
ggcoxzph(test.ph1)

dati$timeDFS<-ifelse(dati$DFS==0,1,dati$DFS)
m1<-coxphw(Surv(timeDFS,relapse) ~ ERbeta5nucleo, data = dati, template = "AHR")
summary(m1)

m1<-coxphw(Surv(timeDFS,relapse) ~ ERbeta5nucleo*factor(BRCA), data = dati, template = "AHR")
summary(m1)


#########################
#ERBeta5_CitoScore_num#
#########################

#MORTALITA'
by(dati$ERbeta5cito,dati$status_alive_dead,summary)
m2=coxph(Surv(dati$OS_diagn, dati$status_alive_dead) ~ dati$ERbeta5cito, data=dati)
summary(m2)
test.ph2 <- cox.zph(m2); test.ph2
library(survminer)
ggcoxzph(test.ph2)

m2=coxph(Surv(dati$OS_diagn, dati$status_alive_dead) ~ dati$ERbeta5cito*factor(BRCA), data=dati)
summary(m2)
test.ph2 <- cox.zph(m2); test.ph2
library(survminer)
ggcoxzph(test.ph2)

m2=coxph(Surv(dati$OS_diagn, dati$status_alive_dead) ~ dati$ERbeta5cito*factor(BRCA)+dati$age, data=dati)
summary(m2)
test.ph2 <- cox.zph(m2); test.ph2
library(survminer)
ggcoxzph(test.ph2)

m2 <- coxphw(Surv(OS_diagn, status_alive_dead) ~ ERbeta5cito,data = dati, template = "AHR")
summary(m2)

m2 <- coxphw(Surv(OS_diagn, status_alive_dead) ~ ERbeta5cito*factor(BRCA),data = dati, template = "AHR")
summary(m2)


#RECIDIVA DFS
by(dati$ERbeta5cito,dati$relapse,summary)
m1=coxph(Surv(dati$DFS,dati$relapse) ~ dati$ERbeta5cito, data=dati)
summary(m1)
test.ph1 <- cox.zph(m1); test.ph1
library(survminer)
ggcoxzph(test.ph1)

m1=coxph(Surv(dati$DFS,dati$relapse) ~ dati$ERbeta5cito*factor(BRCA), data=dati)
summary(m1)
test.ph1 <- cox.zph(m1); test.ph1
library(survminer)
ggcoxzph(test.ph1)

m1=coxph(Surv(dati$DFS,dati$relapse) ~ dati$ERbeta5cito*factor(BRCA)+dati$age, data=dati)
summary(m1)
test.ph1 <- cox.zph(m1); test.ph1
library(survminer)
ggcoxzph(test.ph1)

dati$timeDFS<-ifelse(dati$DFS==0,1,dati$DFS)
m1<-coxphw(Surv(timeDFS,relapse) ~ ERbeta5cito, data = dati, template = "AHR")
summary(m1)

m1<-coxphw(Surv(timeDFS,relapse) ~ ERbeta5cito*factor(BRCA), data = dati, template = "AHR")
summary(m1)

############BRCA##############
##############################
#ERalpha/ERbeta1nucleo########
#############################

#MORTALITA'
m2=coxph(Surv(dati$OS_diagn, dati$status_alive_dead) ~ dati$a_b1n, data=dati)
summary(m2)
test.ph2 <- cox.zph(m2); test.ph2
library(survminer)
ggcoxzph(test.ph2)

m3=coxph(Surv(dati$OS_diagn, dati$status_alive_dead) ~ dati$a_b1n*factor(BRCA), data=dati)
summary(m3)
test.ph2 <- cox.zph(m3); test.ph2
library(survminer)
ggcoxzph(test.ph2)

m4=coxph(Surv(dati$OS_diagn, dati$status_alive_dead) ~ dati$a_b1n*factor(BRCA)+dati$age, data=dati)
summary(m4)
test.ph2 <- cox.zph(m4); test.ph2
library(survminer)
ggcoxzph(test.ph2)



#RECIDIVA DFS
by(dati$a_b1n,dati$relapse,summary)
m1=coxph(Surv(dati$DFS,dati$relapse) ~ dati$a_b1n, data=dati)
summary(m1)
test.ph1 <- cox.zph(m1); test.ph1
library(survminer)
ggcoxzph(test.ph1)

m2=coxph(Surv(dati$DFS,dati$relapse) ~ dati$a_b1n*factor(BRCA), data=dati)
summary(m2)
test.ph1 <- cox.zph(m2); test.ph1
library(survminer)
ggcoxzph(test.ph1)

m3=coxph(Surv(dati$DFS,dati$relapse) ~ dati$a_b1n*factor(BRCA)+dati$age, data=dati)
summary(m3)
test.ph1 <- cox.zph(m3); test.ph1
library(survminer)
ggcoxzph(test.ph1)




###################
#alpha_beta2n######
###################


#MORTALITA'
m2=coxph(Surv(dati$OS_diagn, dati$status_alive_dead) ~ dati$a_b2n, data=dati)
summary(m2)
test.ph2 <- cox.zph(m2); test.ph2
library(survminer)
ggcoxzph(test.ph2)

m3=coxph(Surv(dati$OS_diagn, dati$status_alive_dead) ~ dati$a_b2n*factor(BRCA), data=dati)
summary(m3)
test.ph2 <- cox.zph(m3); test.ph2
library(survminer)
ggcoxzph(test.ph2)

m4=coxph(Surv(dati$OS_diagn, dati$status_alive_dead) ~ dati$a_b2n*factor(BRCA)+dati$age, data=dati)
summary(m4)
test.ph2 <- cox.zph(m4); test.ph2
library(survminer)
ggcoxzph(test.ph2)




#RECIDIVA DFS
by(dati$a_b2n,dati$relapse,summary)
m1=coxph(Surv(dati$DFS,dati$relapse) ~ dati$a_b2n, data=dati)
summary(m1)
test.ph1 <- cox.zph(m1); test.ph1
library(survminer)
ggcoxzph(test.ph1)

m1=coxph(Surv(dati$DFS,dati$relapse) ~ dati$a_b2n*factor(BRCA), data=dati)
summary(m1)
test.ph1 <- cox.zph(m1); test.ph1
library(survminer)
ggcoxzph(test.ph1)

m1=coxph(Surv(dati$DFS,dati$relapse) ~ dati$a_b2n*factor(BRCA)+dati$age, data=dati)
summary(m1)
test.ph1 <- cox.zph(m1); test.ph1
library(survminer)
ggcoxzph(test.ph1)




###################
#alpha_beta5n######
###################

#MORTALITA'
m4=coxph(Surv(dati$OS_diagn, dati$status_alive_dead) ~ dati$a_b5n, data=dati)
summary(m4)
test.ph4 <- cox.zph(m4); test.ph4
library(survminer)
ggcoxzph(test.ph4)

m4=coxph(Surv(dati$OS_diagn, dati$status_alive_dead) ~ dati$a_b5n*factor(BRCA), data=dati)
summary(m4)
test.ph4 <- cox.zph(m4); test.ph4
library(survminer)
ggcoxzph(test.ph4)

m4=coxph(Surv(dati$OS_diagn, dati$status_alive_dead) ~ dati$a_b5n*factor(BRCA)+dati$age, data=dati)
summary(m4)
test.ph4 <- cox.zph(m4); test.ph4
library(survminer)
ggcoxzph(test.ph4)


#RECIDIVA DFS
by(dati$a_b5n,dati$relapse,summary)
m4=coxph(Surv(dati$DFS,dati$relapse) ~ dati$a_b5n, data=dati)
summary(m4)
test.ph4 <- cox.zph(m4); test.ph4
library(survminer)
ggcoxzph(test.ph4)

m4=coxph(Surv(dati$DFS,dati$relapse) ~ dati$a_b5n*factor(BRCA), data=dati)
summary(m4)
test.ph4 <- cox.zph(m4); test.ph4
library(survminer)
ggcoxzph(test.ph4)

m4=coxph(Surv(dati$DFS,dati$relapse) ~ dati$a_b5n*factor(BRCA)+dati$age, data=dati)
summary(m4)
test.ph4 <- cox.zph(m4); test.ph4
library(survminer)
ggcoxzph(test.ph4)



#########################
#PFS dall'ultimo platino#
#########################

#MORTALITA'
by(dati$PFSultimoPLATINO,dati$status_alive_dead,summary)
m2=coxph(Surv(dati$OS_diagn, dati$status_alive_dead) ~ dati$PFSultimoPLATINO, data=dati)
summary(m2)
test.ph2 <- cox.zph(m2); test.ph2
library(survminer)
ggcoxzph(test.ph2)

m2=coxph(Surv(dati$OS_diagn, dati$status_alive_dead) ~ dati$PFSultimoPLATINO*factor(BRCA), data=dati)
summary(m2)
test.ph2 <- cox.zph(m2); test.ph2
library(survminer)
ggcoxzph(test.ph2)

m2=coxph(Surv(dati$OS_diagn, dati$status_alive_dead) ~ dati$PFSultimoPLATINO*factor(BRCA)+dati$age, data=dati)
summary(m2)
test.ph2 <- cox.zph(m2); test.ph2
library(survminer)
ggcoxzph(test.ph2)


#RECIDIVA DFS
by(dati$PFSultimoPLATINO,dati$relapse,summary)
m1=coxph(Surv(dati$DFS,dati$relapse) ~ dati$PFSultimoPLATINO, data=dati)
summary(m1)
test.ph1 <- cox.zph(m1); test.ph1
library(survminer)
ggcoxzph(test.ph1)

m1=coxph(Surv(dati$DFS,dati$relapse) ~ dati$PFSultimoPLATINO*factor(BRCA), data=dati)
summary(m1)
test.ph1 <- cox.zph(m1); test.ph1
library(survminer)
ggcoxzph(test.ph1)

m1=coxph(Surv(dati$DFS,dati$relapse) ~ dati$PFSultimoPLATINO*factor(BRCA)+dati$age, data=dati)
summary(m1)
test.ph1 <- cox.zph(m1); test.ph1
library(survminer)
ggcoxzph(test.ph1)



#########################
#platinum resistance#
#########################

#MORTALITA'
by(dati$PLARES,dati$status_alive_dead,summary)
m2=coxph(Surv(dati$OS_diagn, dati$status_alive_dead) ~ dati$PLARES, data=dati)
summary(m2)
test.ph2 <- cox.zph(m2); test.ph2
library(survminer)
ggcoxzph(test.ph2)

m2=coxph(Surv(dati$OS_diagn, dati$status_alive_dead) ~ dati$PLARES*factor(BRCA), data=dati)
summary(m2)
test.ph2 <- cox.zph(m2); test.ph2
library(survminer)
ggcoxzph(test.ph2)

m2=coxph(Surv(dati$OS_diagn, dati$status_alive_dead) ~ dati$PLARES*factor(BRCA)+dati$age, data=dati)
summary(m2)
test.ph2 <- cox.zph(m2); test.ph2
library(survminer)
ggcoxzph(test.ph2)


#RECIDIVA DFS
by(dati$PLARES,dati$relapse,summary)
m1=coxph(Surv(dati$DFS,dati$relapse) ~ dati$PLARES, data=dati)
summary(m1)
test.ph1 <- cox.zph(m1); test.ph1
library(survminer)
ggcoxzph(test.ph1)

m1=coxph(Surv(dati$DFS,dati$relapse) ~ dati$PLARES*factor(BRCA), data=dati)
summary(m1)
test.ph1 <- cox.zph(m1); test.ph1
library(survminer)
ggcoxzph(test.ph1)

m1=coxph(Surv(dati$DFS,dati$relapse) ~ dati$PLARES*factor(BRCA)+dati$age, data=dati)
summary(m1)
test.ph1 <- cox.zph(m1); test.ph1
library(survminer)
ggcoxzph(test.ph1)

#########################
#ascite#
#########################

#MORTALITA'
by(dati$ascite,dati$status_alive_dead,summary)
m2=coxph(Surv(dati$OS_diagn, dati$status_alive_dead) ~ dati$ascite, data=dati)
summary(m2)
test.ph2 <- cox.zph(m2); test.ph2
library(survminer)
ggcoxzph(test.ph2)

m2=coxph(Surv(dati$OS_diagn, dati$status_alive_dead) ~ dati$ascite*factor(BRCA), data=dati)
summary(m2)
test.ph2 <- cox.zph(m2); test.ph2
library(survminer)
ggcoxzph(test.ph2)

m2=coxph(Surv(dati$OS_diagn, dati$status_alive_dead) ~ dati$ascite*factor(BRCA)+dati$age, data=dati)
summary(m2)
test.ph2 <- cox.zph(m2); test.ph2
library(survminer)
ggcoxzph(test.ph2)


#RECIDIVA DFS
by(dati$ascite,dati$relapse,summary)
m1=coxph(Surv(dati$DFS,dati$relapse) ~ dati$ascite, data=dati)
summary(m1)
test.ph1 <- cox.zph(m1); test.ph1
library(survminer)
ggcoxzph(test.ph1)

m1=coxph(Surv(dati$DFS,dati$relapse) ~ dati$ascite*factor(BRCA), data=dati)
summary(m1)
test.ph1 <- cox.zph(m1); test.ph1
library(survminer)
ggcoxzph(test.ph1)

m1=coxph(Surv(dati$DFS,dati$relapse) ~ dati$ascite*factor(BRCA)+dati$age, data=dati)
summary(m1)
test.ph1 <- cox.zph(m1); test.ph1
library(survminer)
ggcoxzph(test.ph1)


#########################
#chir_pi_lps_score#
#########################

#MORTALITA'
by(dati$chir_pi,dati$status_alive_dead,summary)
m2=coxph(Surv(dati$OS_diagn, dati$status_alive_dead) ~ dati$chir_pi, data=dati)
summary(m2)
test.ph2 <- cox.zph(m2); test.ph2
library(survminer)
ggcoxzph(test.ph2)

m2=coxph(Surv(dati$OS_diagn, dati$status_alive_dead) ~ dati$chir_pi*factor(BRCA), data=dati)
summary(m2)
test.ph2 <- cox.zph(m2); test.ph2
library(survminer)
ggcoxzph(test.ph2)

m2=coxph(Surv(dati$OS_diagn, dati$status_alive_dead) ~ dati$chir_pi*factor(BRCA)+dati$age, data=dati)
summary(m2)
test.ph2 <- cox.zph(m2); test.ph2
library(survminer)
ggcoxzph(test.ph2)


#RECIDIVA DFS
by(dati$chir_pi,dati$relapse,summary)
m1=coxph(Surv(dati$DFS,dati$relapse) ~ dati$chir_pi, data=dati)
summary(m1)
test.ph1 <- cox.zph(m1); test.ph1
library(survminer)
ggcoxzph(test.ph1)

m1=coxph(Surv(dati$DFS,dati$relapse) ~ dati$chir_pi*factor(BRCA), data=dati)
summary(m1)
test.ph1 <- cox.zph(m1); test.ph1
library(survminer)
ggcoxzph(test.ph1)

m1=coxph(Surv(dati$DFS,dati$relapse) ~ dati$chir_pi*factor(BRCA)+dati$age, data=dati)
summary(m1)
test.ph1 <- cox.zph(m1); test.ph1
library(survminer)
ggcoxzph(test.ph1)


#########################
#BMI####################
#########################

#MORTALITA'
by(dati$BMI,dati$status_alive_dead,summary)
m2=coxph(Surv(dati$OS_diagn, dati$status_alive_dead) ~ dati$BMI, data=dati)
summary(m2)
test.ph2 <- cox.zph(m2); test.ph2
library(survminer)
ggcoxzph(test.ph2)

m2=coxph(Surv(dati$OS_diagn, dati$status_alive_dead) ~ dati$BMI*factor(BRCA), data=dati)
summary(m2)
test.ph2 <- cox.zph(m2); test.ph2
library(survminer)
ggcoxzph(test.ph2)

m2=coxph(Surv(dati$OS_diagn, dati$status_alive_dead) ~ dati$BMI*factor(BRCA)+dati$age, data=dati)
summary(m2)
test.ph2 <- cox.zph(m2); test.ph2
library(survminer)
ggcoxzph(test.ph2)


#RECIDIVA DFS
by(dati$BMI,dati$relapse,summary)
m1=coxph(Surv(dati$DFS,dati$relapse) ~ dati$BMI, data=dati)
summary(m1)
test.ph1 <- cox.zph(m1); test.ph1
library(survminer)
ggcoxzph(test.ph1)

m1=coxph(Surv(dati$DFS,dati$relapse) ~ dati$BMI*factor(BRCA), data=dati)
summary(m1)
test.ph1 <- cox.zph(m1); test.ph1
library(survminer)
ggcoxzph(test.ph1)

m1=coxph(Surv(dati$DFS,dati$relapse) ~ dati$BMI*factor(BRCA)+dati$age, data=dati)
summary(m1)
test.ph1 <- cox.zph(m1); test.ph1
library(survminer)
ggcoxzph(test.ph1)


#########################
#menopausa##############
#########################

#MORTALITA'
by(dati$menopausa,dati$status_alive_dead,summary)
m2=coxph(Surv(dati$OS_diagn, dati$status_alive_dead) ~ dati$menopausa, data=dati)
summary(m2)
test.ph2 <- cox.zph(m2); test.ph2
library(survminer)
ggcoxzph(test.ph2)

m2=coxph(Surv(dati$OS_diagn, dati$status_alive_dead) ~ dati$menopausa*factor(BRCA), data=dati)
summary(m2)
test.ph2 <- cox.zph(m2); test.ph2
library(survminer)
ggcoxzph(test.ph2)

m2=coxph(Surv(dati$OS_diagn, dati$status_alive_dead) ~ dati$menopausa*factor(BRCA)+dati$age, data=dati)
summary(m2)
test.ph2 <- cox.zph(m2); test.ph2
library(survminer)
ggcoxzph(test.ph2)


#RECIDIVA DFS
by(dati$menopausa,dati$relapse,summary)
m1=coxph(Surv(dati$DFS,dati$relapse) ~ dati$menopausa, data=dati)
summary(m1)
test.ph1 <- cox.zph(m1); test.ph1
library(survminer)
ggcoxzph(test.ph1)

m1=coxph(Surv(dati$DFS,dati$relapse) ~ dati$menopausa*factor(BRCA), data=dati)
summary(m1)
test.ph1 <- cox.zph(m1); test.ph1
library(survminer)
ggcoxzph(test.ph1)

m1=coxph(Surv(dati$DFS,dati$relapse) ~ dati$menopausa*factor(BRCA)+dati$age, data=dati)
summary(m1)
test.ph1 <- cox.zph(m1); test.ph1
library(survminer)
ggcoxzph(test.ph1)


#########################
#ct_beva#
#########################

#MORTALITA'
by(dati$ct_Beva,dati$status_alive_dead,summary)
m2=coxph(Surv(dati$OS_diagn, dati$status_alive_dead) ~ dati$ct_Beva, data=dati)
summary(m2)
test.ph2 <- cox.zph(m2); test.ph2
library(survminer)
ggcoxzph(test.ph2)

m2=coxph(Surv(dati$OS_diagn, dati$status_alive_dead) ~ dati$ct_Beva*factor(BRCA), data=dati)
summary(m2)
test.ph2 <- cox.zph(m2); test.ph2
library(survminer)
ggcoxzph(test.ph2)

m2=coxph(Surv(dati$OS_diagn, dati$status_alive_dead) ~ dati$ct_Beva*factor(BRCA)+dati$age, data=dati)
summary(m2)
test.ph2 <- cox.zph(m2); test.ph2
library(survminer)
ggcoxzph(test.ph2)


#RECIDIVA DFS
by(dati$ct_Beva,dati$relapse,summary)
m1=coxph(Surv(dati$DFS,dati$relapse) ~ dati$ct_Beva, data=dati)
summary(m1)
test.ph1 <- cox.zph(m1); test.ph1
library(survminer)
ggcoxzph(test.ph1)

m1=coxph(Surv(dati$DFS,dati$relapse) ~ dati$ct_Beva*factor(BRCA), data=dati)
summary(m1)
test.ph1 <- cox.zph(m1); test.ph1
library(survminer)
ggcoxzph(test.ph1)

m1=coxph(Surv(dati$DFS,dati$relapse) ~ dati$ct_Beva*factor(BRCA)+dati$age, data=dati)
summary(m1)
test.ph1 <- cox.zph(m1); test.ph1
library(survminer)
ggcoxzph(test.ph1)


#################################
#citorid_cm#CITORIDOTTE SI/NO####
#################################

#MORTALITA'
by(dati$citorid_cm,dati$status_alive_dead,summary)
m2=coxph(Surv(dati$OS_diagn, dati$status_alive_dead) ~ dati$citorid_cm, data=dati)
summary(m2)
test.ph2 <- cox.zph(m2); test.ph2
library(survminer)
ggcoxzph(test.ph2)

m2=coxph(Surv(dati$OS_diagn, dati$status_alive_dead) ~ dati$citorid_cm*factor(BRCA), data=dati)
summary(m2)
test.ph2 <- cox.zph(m2); test.ph2
library(survminer)
ggcoxzph(test.ph2)

m2=coxph(Surv(dati$OS_diagn, dati$status_alive_dead) ~ dati$citorid_cm*factor(BRCA)+dati$age, data=dati)
summary(m2)
test.ph2 <- cox.zph(m2); test.ph2
library(survminer)
ggcoxzph(test.ph2)

m2 <- coxphw(Surv(OS_diagn, status_alive_dead) ~ citorid_cm,data = dati, template = "AHR")
summary(m2)

m2 <- coxphw(Surv(OS_diagn, status_alive_dead) ~ citorid_cm*factor(BRCA),data = dati, template = "AHR")
summary(m2)


#RECIDIVA DFS
by(dati$citorid_cm,dati$relapse,summary)
m1=coxph(Surv(dati$DFS,dati$relapse) ~ dati$citorid_cm, data=dati)
summary(m1)
test.ph1 <- cox.zph(m1); test.ph1
library(survminer)
ggcoxzph(test.ph1)

m1=coxph(Surv(dati$DFS,dati$relapse) ~ dati$citorid_cm*factor(BRCA), data=dati)
summary(m1)
test.ph1 <- cox.zph(m1); test.ph1
library(survminer)
ggcoxzph(test.ph1)

m1=coxph(Surv(dati$DFS,dati$relapse) ~ dati$citorid_cm*factor(BRCA)+dati$age, data=dati)
summary(m1)
test.ph1 <- cox.zph(m1); test.ph1
library(survminer)
ggcoxzph(test.ph1)

dati$timeDFS<-ifelse(dati$DFS==0,1,dati$DFS)
m1<-coxphw(Surv(timeDFS,relapse) ~ citorid_cm, data = dati, template = "AHR")
summary(m1)

m1<-coxphw(Surv(timeDFS,relapse) ~ citorid_cm*factor(BRCA), data = dati, template = "AHR")
summary(m1)









#########################################
###MODELLO DI COX BRCA_dic####
#########################################

##########
#BRCA_dic#
##########

#MORTALITA'
addmargins(table(dati$status_alive_dead,BRCA_dic))
by(dati$OS_diagn,BRCA_dic,summary)
m0=coxph(Surv(dati$OS_diagn,dati$status_alive_dead) ~ factor(BRCA_dic), data=dati)
summary(m0)
test.ph0 <- cox.zph(m0); test.ph0
library(survminer)
ggcoxzph(test.ph0)

m0<- coxphw(Surv(OS_diagn,status_alive_dead) ~ factor(BRCA_dic), data = dati, template = "AHR")
summary(m0)

#RECIDIVA DFS
addmargins(table(dati$relapse,BRCA_dic))
by(dati$DFS,BRCA_dic,summary)
m1=coxph(Surv(dati$DFS,dati$relapse) ~ factor(BRCA_dic), data=dati)
summary(m1)
test.ph1 <- cox.zph(m1); test.ph1
library(survminer)
ggcoxzph(test.ph1)

dati$timeDFS<-ifelse(dati$DFS==0,1,dati$DFS)
m1<-coxphw(Surv(timeDFS,relapse) ~ factor(BRCA_dic), data = dati, template = "AHR")
summary(m1)


#################
#AR nucleo score#
#################

#MORTALITA'
by(dati$ARnucleo,dati$status_alive_dead,summary)
m2=coxph(Surv(dati$OS_diagn, dati$status_alive_dead) ~ dati$ARnucleo, data=dati)
summary(m2)
test.ph2 <- cox.zph(m2); test.ph2
library(survminer)
ggcoxzph(test.ph2)

m2=coxph(Surv(dati$OS_diagn, dati$status_alive_dead) ~ dati$ARnucleo*factor(BRCA_dic), data=dati)
summary(m2)
test.ph2 <- cox.zph(m2); test.ph2
library(survminer)
ggcoxzph(test.ph2)

m2=coxph(Surv(dati$OS_diagn, dati$status_alive_dead) ~ dati$ARnucleo*factor(BRCA_dic)+dati$age, data=dati)
summary(m2)
test.ph2 <- cox.zph(m2); test.ph2
library(survminer)
ggcoxzph(test.ph2)


#RECIDIVA DFS
by(dati$ARnucleo,dati$relapse,summary)
m1=coxph(Surv(dati$DFS,dati$relapse) ~ dati$ARnucleo, data=dati)
summary(m1)
test.ph1 <- cox.zph(m1); test.ph1
library(survminer)
ggcoxzph(test.ph1)

m1=coxph(Surv(dati$DFS,dati$relapse) ~ dati$ARnucleo*factor(BRCA_dic), data=dati)
summary(m1)
test.ph1 <- cox.zph(m1); test.ph1
library(survminer)
ggcoxzph(test.ph1)

m1=coxph(Surv(dati$DFS,dati$relapse) ~ dati$ARnucleo*factor(BRCA_dic)+dati$age, data=dati)
summary(m1)
test.ph1 <- cox.zph(m1); test.ph1
library(survminer)
ggcoxzph(test.ph1)



##############
#PR_Score_num#
##############


#MORTALITA'
by(dati$PR,dati$status_alive_dead,summary)
m2=coxph(Surv(dati$OS_diagn, dati$status_alive_dead) ~ dati$PR, data=dati)
summary(m2)
test.ph2 <- cox.zph(m2); test.ph2
library(survminer)
ggcoxzph(test.ph2)

m2=coxph(Surv(dati$OS_diagn, dati$status_alive_dead) ~ dati$PR*factor(BRCA_dic), data=dati)
summary(m2)
test.ph2 <- cox.zph(m2); test.ph2
library(survminer)
ggcoxzph(test.ph2)

m2=coxph(Surv(dati$OS_diagn, dati$status_alive_dead) ~ dati$PR*factor(BRCA_dic)+dati$age, data=dati)
summary(m2)
test.ph2 <- cox.zph(m2); test.ph2
library(survminer)
ggcoxzph(test.ph2)


#RECIDIVA DFS
by(dati$PR,dati$relapse,summary)
m1=coxph(Surv(dati$DFS,dati$relapse) ~ dati$PR, data=dati)
summary(m1)
test.ph1 <- cox.zph(m1); test.ph1
library(survminer)
ggcoxzph(test.ph1)

m1=coxph(Surv(dati$DFS,dati$relapse) ~ dati$PR*factor(BRCA_dic), data=dati)
summary(m1)
test.ph1 <- cox.zph(m1); test.ph1
library(survminer)
ggcoxzph(test.ph1)

m1=coxph(Surv(dati$DFS,dati$relapse) ~ dati$PR*factor(BRCA_dic)+dati$age, data=dati)
summary(m1)
test.ph1 <- cox.zph(m1); test.ph1
library(survminer)
ggcoxzph(test.ph1)



###################
#ERAlpha_Score_num#
###################

#MORTALITA'
by(dati$ERalpha,dati$status_alive_dead,summary)
m4=coxph(Surv(dati$OS_diagn, dati$status_alive_dead) ~ dati$ERalpha, data=dati)
summary(m4)
test.ph4 <- cox.zph(m4); test.ph4
library(survminer)
ggcoxzph(test.ph4)

m4=coxph(Surv(dati$OS_diagn, dati$status_alive_dead) ~ dati$ERalpha*factor(BRCA_dic), data=dati)
summary(m4)
test.ph4 <- cox.zph(m4); test.ph4
library(survminer)
ggcoxzph(test.ph4)

m4=coxph(Surv(dati$OS_diagn, dati$status_alive_dead) ~ dati$ERalpha*factor(BRCA_dic)+dati$age, data=dati)
summary(m4)
test.ph4 <- cox.zph(m4); test.ph4
library(survminer)
ggcoxzph(test.ph4)


#RECIDIVA DFS
by(dati$ERalpha,dati$relapse,summary)
m4=coxph(Surv(dati$DFS,dati$relapse) ~ dati$ERalpha, data=dati)
summary(m4)
test.ph4 <- cox.zph(m4); test.ph4
library(survminer)
ggcoxzph(test.ph4)

m4=coxph(Surv(dati$DFS,dati$relapse) ~ dati$ERalpha*factor(BRCA_dic), data=dati)
summary(m4)
test.ph4 <- cox.zph(m4); test.ph4
library(survminer)
ggcoxzph(test.ph4)

m4=coxph(Surv(dati$DFS,dati$relapse) ~ dati$ERalpha*factor(BRCA_dic)+dati$age, data=dati)
summary(m4)
test.ph4 <- cox.zph(m4); test.ph4
library(survminer)
ggcoxzph(test.ph4)


#########################
#ERBeta1_NucleoScore_num#
#########################

#MORTALITA'
by(dati$ERbeta1nucleo,dati$status_alive_dead,summary)
m2=coxph(Surv(dati$OS_diagn, dati$status_alive_dead) ~ dati$ERbeta1nucleo, data=dati)
summary(m2)
test.ph2 <- cox.zph(m2); test.ph2
library(survminer)
ggcoxzph(test.ph2)

m2=coxph(Surv(dati$OS_diagn, dati$status_alive_dead) ~ dati$ERbeta1nucleo*factor(BRCA_dic), data=dati)
summary(m2)
test.ph2 <- cox.zph(m2); test.ph2
library(survminer)
ggcoxzph(test.ph2)

m2=coxph(Surv(dati$OS_diagn, dati$status_alive_dead) ~ dati$ERbeta1nucleo*factor(BRCA_dic)+dati$age, data=dati)
summary(m2)
test.ph2 <- cox.zph(m2); test.ph2
library(survminer)
ggcoxzph(test.ph2)


#RECIDIVA DFS
by(dati$ERbeta1nucleo,dati$relapse,summary)
m1=coxph(Surv(dati$DFS,dati$relapse) ~ dati$ERbeta1nucleo, data=dati)
summary(m1)
test.ph1 <- cox.zph(m1); test.ph1
library(survminer)
ggcoxzph(test.ph1)

m1=coxph(Surv(dati$DFS,dati$relapse) ~ dati$ERbeta1nucleo*factor(BRCA_dic), data=dati)
summary(m1)
test.ph1 <- cox.zph(m1); test.ph1
library(survminer)
ggcoxzph(test.ph1)

m1=coxph(Surv(dati$DFS,dati$relapse) ~ dati$ERbeta1nucleo*factor(BRCA_dic)+dati$age, data=dati)
summary(m1)
test.ph1 <- cox.zph(m1); test.ph1
library(survminer)
ggcoxzph(test.ph1)


#######################
#ERBeta1_CitoScore_num#
#######################

#MORTALITA'
by(dati$ERbeta1cito,dati$status_alive_dead,summary)
m2=coxph(Surv(dati$OS_diagn, dati$status_alive_dead) ~ dati$ERbeta1cito, data=dati)
summary(m2)
test.ph2 <- cox.zph(m2); test.ph2
library(survminer)
ggcoxzph(test.ph2)

m2=coxph(Surv(dati$OS_diagn, dati$status_alive_dead) ~ dati$ERbeta1cito*factor(BRCA_dic), data=dati)
summary(m2)
test.ph2 <- cox.zph(m2); test.ph2
library(survminer)
ggcoxzph(test.ph2)

m2=coxph(Surv(dati$OS_diagn, dati$status_alive_dead) ~ dati$ERbeta1cito*factor(BRCA_dic)+dati$age, data=dati)
summary(m2)
test.ph2 <- cox.zph(m2); test.ph2
library(survminer)
ggcoxzph(test.ph2)


#RECIDIVA DFS
by(dati$ERbeta1cito,dati$relapse,summary)
m1=coxph(Surv(dati$DFS,dati$relapse) ~ dati$ERbeta1cito, data=dati)
summary(m1)
test.ph1 <- cox.zph(m1); test.ph1
library(survminer)
ggcoxzph(test.ph1)

m1=coxph(Surv(dati$DFS,dati$relapse) ~ dati$ERbeta1cito*factor(BRCA_dic), data=dati)
summary(m1)
test.ph1 <- cox.zph(m1); test.ph1
library(survminer)
ggcoxzph(test.ph1)

m1=coxph(Surv(dati$DFS,dati$relapse) ~ dati$ERbeta1cito*factor(BRCA_dic)+dati$age, data=dati)
summary(m1)
test.ph1 <- cox.zph(m1); test.ph1
library(survminer)
ggcoxzph(test.ph1)



#########################
#ERBeta2_NucleoScore_num#
#########################

#MORTALITA'
by(dati$ERbeta2nucleo,dati$status_alive_dead,summary)
m2=coxph(Surv(dati$OS_diagn, dati$status_alive_dead) ~ dati$ERbeta2nucleo, data=dati)
summary(m2)
test.ph2 <- cox.zph(m2); test.ph2
library(survminer)
ggcoxzph(test.ph2)

m2=coxph(Surv(dati$OS_diagn, dati$status_alive_dead) ~ dati$ERbeta2nucleo*factor(BRCA_dic), data=dati)
summary(m2)
test.ph2 <- cox.zph(m2); test.ph2
library(survminer)
ggcoxzph(test.ph2)

m2=coxph(Surv(dati$OS_diagn, dati$status_alive_dead) ~ dati$ERbeta2nucleo*factor(BRCA_dic)+dati$age, data=dati)
summary(m2)
test.ph2 <- cox.zph(m2); test.ph2
library(survminer)
ggcoxzph(test.ph2)


#RECIDIVA DFS
by(dati$ERbeta2nucleo,dati$relapse,summary)
m1=coxph(Surv(dati$DFS,dati$relapse) ~ dati$ERbeta2nucleo, data=dati)
summary(m1)
test.ph1 <- cox.zph(m1); test.ph1
library(survminer)
ggcoxzph(test.ph1)

m1=coxph(Surv(dati$DFS,dati$relapse) ~ dati$ERbeta2nucleo*factor(BRCA_dic), data=dati)
summary(m1)
test.ph1 <- cox.zph(m1); test.ph1
library(survminer)
ggcoxzph(test.ph1)

m1=coxph(Surv(dati$DFS,dati$relapse) ~ dati$ERbeta2nucleo*factor(BRCA_dic)+dati$age, data=dati)
summary(m1)
test.ph1 <- cox.zph(m1); test.ph1
library(survminer)
ggcoxzph(test.ph1)

dati$timeDFS<-ifelse(dati$DFS==0,1,dati$DFS)
m1<-coxphw(Surv(timeDFS,relapse) ~ ERbeta2nucleo, data = dati, template = "AHR")
summary(m1)
m1<-coxphw(Surv(timeDFS,relapse) ~ ERbeta2nucleo*factor(BRCA_dic), data = dati, template = "AHR")
summary(m1)


#######################
#ERBeta2_CitoScore_num#
#######################


#MORTALITA'
by(dati$ERbeta2cito,dati$status_alive_dead,summary)
m2=coxph(Surv(dati$OS_diagn, dati$status_alive_dead) ~ dati$ERbeta2cito, data=dati)
summary(m2)
test.ph2 <- cox.zph(m2); test.ph2
library(survminer)
ggcoxzph(test.ph2)

m2=coxph(Surv(dati$OS_diagn, dati$status_alive_dead) ~ dati$ERbeta2cito*factor(BRCA_dic), data=dati)
summary(m2)
test.ph2 <- cox.zph(m2); test.ph2
library(survminer)
ggcoxzph(test.ph2)

m2=coxph(Surv(dati$OS_diagn, dati$status_alive_dead) ~ dati$ERbeta2cito*factor(BRCA_dic)+dati$age, data=dati)
summary(m2)
test.ph2 <- cox.zph(m2); test.ph2
library(survminer)
ggcoxzph(test.ph2)



#RECIDIVA DFS
by(dati$ERbeta2cito,dati$relapse,summary)
m1=coxph(Surv(dati$DFS,dati$relapse) ~ dati$ERbeta2cito, data=dati)
summary(m1)
test.ph1 <- cox.zph(m1); test.ph1
library(survminer)
ggcoxzph(test.ph1)

m1=coxph(Surv(dati$DFS,dati$relapse) ~ dati$ERbeta2cito*factor(BRCA_dic), data=dati)
summary(m1)
test.ph1 <- cox.zph(m1); test.ph1
library(survminer)
ggcoxzph(test.ph1)

m1=coxph(Surv(dati$DFS,dati$relapse) ~ dati$ERbeta2cito*factor(BRCA_dic)+dati$age, data=dati)
summary(m1)
test.ph1 <- cox.zph(m1); test.ph1
library(survminer)
ggcoxzph(test.ph1)

dati$timeDFS<-ifelse(dati$DFS==0,1,dati$DFS)
m1<-coxphw(Surv(timeDFS,relapse) ~ ERbeta2cito, data = dati, template = "AHR")
summary(m1)

m1<-coxphw(Surv(timeDFS,relapse) ~ ERbeta2cito*factor(BRCA_dic), data = dati, template = "AHR")
summary(m1)


#########################
#ERBeta5_NucleoScore_num#
#########################

#MORTALITA'
by(dati$ERbeta5nucleo,dati$status_alive_dead,summary)
m2=coxph(Surv(dati$OS_diagn, dati$status_alive_dead) ~ dati$ERbeta5nucleo, data=dati)
summary(m2)
test.ph2 <- cox.zph(m2); test.ph2
library(survminer)
ggcoxzph(test.ph2)

m2=coxph(Surv(dati$OS_diagn, dati$status_alive_dead) ~ dati$ERbeta5nucleo*factor(BRCA_dic), data=dati)
summary(m2)
test.ph2 <- cox.zph(m2); test.ph2
library(survminer)
ggcoxzph(test.ph2)

m2=coxph(Surv(dati$OS_diagn, dati$status_alive_dead) ~ dati$ERbeta5nucleo*factor(BRCA_dic)+dati$age, data=dati)
summary(m2)
test.ph2 <- cox.zph(m2); test.ph2
library(survminer)
ggcoxzph(test.ph2)

m2 <- coxphw(Surv(OS_diagn, status_alive_dead) ~ ERbeta5nucleo,data = dati, template = "AHR")
summary(m2)

m2 <- coxphw(Surv(OS_diagn, status_alive_dead) ~ ERbeta5nucleo*factor(BRCA_dic),data = dati, template = "AHR")
summary(m2)


#RECIDIVA DFS
by(dati$ERbeta5nucleo,dati$relapse,summary)
m1=coxph(Surv(dati$DFS,dati$relapse) ~ dati$ERbeta5nucleo, data=dati)
summary(m1)
test.ph1 <- cox.zph(m1); test.ph1
library(survminer)
ggcoxzph(test.ph1)

m1=coxph(Surv(dati$DFS,dati$relapse) ~ dati$ERbeta5nucleo*factor(BRCA_dic), data=dati)
summary(m1)
test.ph1 <- cox.zph(m1); test.ph1
library(survminer)
ggcoxzph(test.ph1)

m1=coxph(Surv(dati$DFS,dati$relapse) ~ dati$ERbeta5nucleo*factor(BRCA_dic)+dati$age, data=dati)
summary(m1)
test.ph1 <- cox.zph(m1); test.ph1
library(survminer)
ggcoxzph(test.ph1)

dati$timeDFS<-ifelse(dati$DFS==0,1,dati$DFS)
m1<-coxphw(Surv(timeDFS,relapse) ~ ERbeta5nucleo, data = dati, template = "AHR")
summary(m1)

m1<-coxphw(Surv(timeDFS,relapse) ~ ERbeta5nucleo*factor(BRCA_dic), data = dati, template = "AHR")
summary(m1)


#########################
#ERBeta5_CitoScore_num#
#########################

#MORTALITA'
by(dati$ERbeta5cito,dati$status_alive_dead,summary)
m2=coxph(Surv(dati$OS_diagn, dati$status_alive_dead) ~ dati$ERbeta5cito, data=dati)
summary(m2)
test.ph2 <- cox.zph(m2); test.ph2
library(survminer)
ggcoxzph(test.ph2)

m2=coxph(Surv(dati$OS_diagn, dati$status_alive_dead) ~ dati$ERbeta5cito*factor(BRCA_dic), data=dati)
summary(m2)
test.ph2 <- cox.zph(m2); test.ph2
library(survminer)
ggcoxzph(test.ph2)

m2=coxph(Surv(dati$OS_diagn, dati$status_alive_dead) ~ dati$ERbeta5cito*factor(BRCA_dic)+dati$age, data=dati)
summary(m2)
test.ph2 <- cox.zph(m2); test.ph2
library(survminer)
ggcoxzph(test.ph2)

m2 <- coxphw(Surv(OS_diagn, status_alive_dead) ~ ERbeta5cito,data = dati, template = "AHR")
summary(m2)

m2 <- coxphw(Surv(OS_diagn, status_alive_dead) ~ ERbeta5cito*factor(BRCA_dic),data = dati, template = "AHR")
summary(m2)


#RECIDIVA DFS
by(dati$ERbeta5cito,dati$relapse,summary)
m1=coxph(Surv(dati$DFS,dati$relapse) ~ dati$ERbeta5cito, data=dati)
summary(m1)
test.ph1 <- cox.zph(m1); test.ph1
library(survminer)
ggcoxzph(test.ph1)

m1=coxph(Surv(dati$DFS,dati$relapse) ~ dati$ERbeta5cito*factor(BRCA_dic), data=dati)
summary(m1)
test.ph1 <- cox.zph(m1); test.ph1
library(survminer)
ggcoxzph(test.ph1)

m1=coxph(Surv(dati$DFS,dati$relapse) ~ dati$ERbeta5cito*factor(BRCA_dic)+dati$age, data=dati)
summary(m1)
test.ph1 <- cox.zph(m1); test.ph1
library(survminer)
ggcoxzph(test.ph1)

dati$timeDFS<-ifelse(dati$DFS==0,1,dati$DFS)
m1<-coxphw(Surv(timeDFS,relapse) ~ ERbeta5cito, data = dati, template = "AHR")
summary(m1)

m1<-coxphw(Surv(timeDFS,relapse) ~ ERbeta5cito*factor(BRCA_dic), data = dati, template = "AHR")
summary(m1)

############BRCA_dic##############
##############################
#ERalpha/ERbeta1nucleo########
#############################

#MORTALITA'
m2=coxph(Surv(dati$OS_diagn, dati$status_alive_dead) ~ dati$a_b1n, data=dati)
summary(m2)
test.ph2 <- cox.zph(m2); test.ph2
library(survminer)
ggcoxzph(test.ph2)

m3=coxph(Surv(dati$OS_diagn, dati$status_alive_dead) ~ dati$a_b1n*factor(BRCA_dic), data=dati)
summary(m3)
test.ph2 <- cox.zph(m3); test.ph2
library(survminer)
ggcoxzph(test.ph2)

m4=coxph(Surv(dati$OS_diagn, dati$status_alive_dead) ~ dati$a_b1n*factor(BRCA_dic)+dati$age, data=dati)
summary(m4)
test.ph2 <- cox.zph(m4); test.ph2
library(survminer)
ggcoxzph(test.ph2)



#RECIDIVA DFS
by(dati$a_b1n,dati$relapse,summary)
m1=coxph(Surv(dati$DFS,dati$relapse) ~ dati$a_b1n, data=dati)
summary(m1)
test.ph1 <- cox.zph(m1); test.ph1
library(survminer)
ggcoxzph(test.ph1)

m2=coxph(Surv(dati$DFS,dati$relapse) ~ dati$a_b1n*factor(BRCA_dic), data=dati)
summary(m2)
test.ph1 <- cox.zph(m2); test.ph1
library(survminer)
ggcoxzph(test.ph1)

m3=coxph(Surv(dati$DFS,dati$relapse) ~ dati$a_b1n*factor(BRCA_dic)+dati$age, data=dati)
summary(m3)
test.ph1 <- cox.zph(m3); test.ph1
library(survminer)
ggcoxzph(test.ph1)




###################
#alpha_beta2n######
###################


#MORTALITA'
m2=coxph(Surv(dati$OS_diagn, dati$status_alive_dead) ~ dati$a_b2n, data=dati)
summary(m2)
test.ph2 <- cox.zph(m2); test.ph2
library(survminer)
ggcoxzph(test.ph2)

m3=coxph(Surv(dati$OS_diagn, dati$status_alive_dead) ~ dati$a_b2n*factor(BRCA_dic), data=dati)
summary(m3)
test.ph2 <- cox.zph(m3); test.ph2
library(survminer)
ggcoxzph(test.ph2)

m4=coxph(Surv(dati$OS_diagn, dati$status_alive_dead) ~ dati$a_b2n*factor(BRCA_dic)+dati$age, data=dati)
summary(m4)
test.ph2 <- cox.zph(m4); test.ph2
library(survminer)
ggcoxzph(test.ph2)




#RECIDIVA DFS
by(dati$a_b2n,dati$relapse,summary)
m1=coxph(Surv(dati$DFS,dati$relapse) ~ dati$a_b2n, data=dati)
summary(m1)
test.ph1 <- cox.zph(m1); test.ph1
library(survminer)
ggcoxzph(test.ph1)

m1=coxph(Surv(dati$DFS,dati$relapse) ~ dati$a_b2n*factor(BRCA_dic), data=dati)
summary(m1)
test.ph1 <- cox.zph(m1); test.ph1
library(survminer)
ggcoxzph(test.ph1)

m1=coxph(Surv(dati$DFS,dati$relapse) ~ dati$a_b2n*factor(BRCA_dic)+dati$age, data=dati)
summary(m1)
test.ph1 <- cox.zph(m1); test.ph1
library(survminer)
ggcoxzph(test.ph1)




###################
#alpha_beta5n######
###################

#MORTALITA'
m4=coxph(Surv(dati$OS_diagn, dati$status_alive_dead) ~ dati$a_b5n, data=dati)
summary(m4)
test.ph4 <- cox.zph(m4); test.ph4
library(survminer)
ggcoxzph(test.ph4)

m4=coxph(Surv(dati$OS_diagn, dati$status_alive_dead) ~ dati$a_b5n*factor(BRCA_dic), data=dati)
summary(m4)
test.ph4 <- cox.zph(m4); test.ph4
library(survminer)
ggcoxzph(test.ph4)

m4=coxph(Surv(dati$OS_diagn, dati$status_alive_dead) ~ dati$a_b5n*factor(BRCA_dic)+dati$age, data=dati)
summary(m4)
test.ph4 <- cox.zph(m4); test.ph4
library(survminer)
ggcoxzph(test.ph4)


#RECIDIVA DFS
by(dati$a_b5n,dati$relapse,summary)
m4=coxph(Surv(dati$DFS,dati$relapse) ~ dati$a_b5n, data=dati)
summary(m4)
test.ph4 <- cox.zph(m4); test.ph4
library(survminer)
ggcoxzph(test.ph4)

m4=coxph(Surv(dati$DFS,dati$relapse) ~ dati$a_b5n*factor(BRCA_dic), data=dati)
summary(m4)
test.ph4 <- cox.zph(m4); test.ph4
library(survminer)
ggcoxzph(test.ph4)

m4=coxph(Surv(dati$DFS,dati$relapse) ~ dati$a_b5n*factor(BRCA_dic)+dati$age, data=dati)
summary(m4)
test.ph4 <- cox.zph(m4); test.ph4
library(survminer)
ggcoxzph(test.ph4)



#########################
#PFS dall'ultimo platino#
#########################

#MORTALITA'
by(dati$PFSultimoPLATINO,dati$status_alive_dead,summary)
m2=coxph(Surv(dati$OS_diagn, dati$status_alive_dead) ~ dati$PFSultimoPLATINO, data=dati)
summary(m2)
test.ph2 <- cox.zph(m2); test.ph2
library(survminer)
ggcoxzph(test.ph2)

m2=coxph(Surv(dati$OS_diagn, dati$status_alive_dead) ~ dati$PFSultimoPLATINO*factor(BRCA_dic), data=dati)
summary(m2)
test.ph2 <- cox.zph(m2); test.ph2
library(survminer)
ggcoxzph(test.ph2)

m2=coxph(Surv(dati$OS_diagn, dati$status_alive_dead) ~ dati$PFSultimoPLATINO*factor(BRCA_dic)+dati$age, data=dati)
summary(m2)
test.ph2 <- cox.zph(m2); test.ph2
library(survminer)
ggcoxzph(test.ph2)


#RECIDIVA DFS
by(dati$PFSultimoPLATINO,dati$relapse,summary)
m1=coxph(Surv(dati$DFS,dati$relapse) ~ dati$PFSultimoPLATINO, data=dati)
summary(m1)
test.ph1 <- cox.zph(m1); test.ph1
library(survminer)
ggcoxzph(test.ph1)

m1=coxph(Surv(dati$DFS,dati$relapse) ~ dati$PFSultimoPLATINO*factor(BRCA_dic), data=dati)
summary(m1)
test.ph1 <- cox.zph(m1); test.ph1
library(survminer)
ggcoxzph(test.ph1)

m1=coxph(Surv(dati$DFS,dati$relapse) ~ dati$PFSultimoPLATINO*factor(BRCA_dic)+dati$age, data=dati)
summary(m1)
test.ph1 <- cox.zph(m1); test.ph1
library(survminer)
ggcoxzph(test.ph1)



#########################
#platinum resistance#
#########################

#MORTALITA'
by(dati$PLARES,dati$status_alive_dead,summary)
m2=coxph(Surv(dati$OS_diagn, dati$status_alive_dead) ~ dati$PLARES, data=dati)
summary(m2)
test.ph2 <- cox.zph(m2); test.ph2
library(survminer)
ggcoxzph(test.ph2)

m2=coxph(Surv(dati$OS_diagn, dati$status_alive_dead) ~ dati$PLARES*factor(BRCA_dic), data=dati)
summary(m2)
test.ph2 <- cox.zph(m2); test.ph2
library(survminer)
ggcoxzph(test.ph2)

m2=coxph(Surv(dati$OS_diagn, dati$status_alive_dead) ~ dati$PLARES*factor(BRCA_dic)+dati$age, data=dati)
summary(m2)
test.ph2 <- cox.zph(m2); test.ph2
library(survminer)
ggcoxzph(test.ph2)


#RECIDIVA DFS
by(dati$PLARES,dati$relapse,summary)
m1=coxph(Surv(dati$DFS,dati$relapse) ~ dati$PLARES, data=dati)
summary(m1)
test.ph1 <- cox.zph(m1); test.ph1
library(survminer)
ggcoxzph(test.ph1)

m1=coxph(Surv(dati$DFS,dati$relapse) ~ dati$PLARES*factor(BRCA_dic), data=dati)
summary(m1)
test.ph1 <- cox.zph(m1); test.ph1
library(survminer)
ggcoxzph(test.ph1)

m1=coxph(Surv(dati$DFS,dati$relapse) ~ dati$PLARES*factor(BRCA_dic)+dati$age, data=dati)
summary(m1)
test.ph1 <- cox.zph(m1); test.ph1
library(survminer)
ggcoxzph(test.ph1)

#########################
#ascite#
#########################

#MORTALITA'
by(dati$ascite,dati$status_alive_dead,summary)
m2=coxph(Surv(dati$OS_diagn, dati$status_alive_dead) ~ dati$ascite, data=dati)
summary(m2)
test.ph2 <- cox.zph(m2); test.ph2
library(survminer)
ggcoxzph(test.ph2)

m2=coxph(Surv(dati$OS_diagn, dati$status_alive_dead) ~ dati$ascite*factor(BRCA_dic), data=dati)
summary(m2)
test.ph2 <- cox.zph(m2); test.ph2
library(survminer)
ggcoxzph(test.ph2)

m2=coxph(Surv(dati$OS_diagn, dati$status_alive_dead) ~ dati$ascite*factor(BRCA_dic)+dati$age, data=dati)
summary(m2)
test.ph2 <- cox.zph(m2); test.ph2
library(survminer)
ggcoxzph(test.ph2)


#RECIDIVA DFS
by(dati$ascite,dati$relapse,summary)
m1=coxph(Surv(dati$DFS,dati$relapse) ~ dati$ascite, data=dati)
summary(m1)
test.ph1 <- cox.zph(m1); test.ph1
library(survminer)
ggcoxzph(test.ph1)

m1=coxph(Surv(dati$DFS,dati$relapse) ~ dati$ascite*factor(BRCA_dic), data=dati)
summary(m1)
test.ph1 <- cox.zph(m1); test.ph1
library(survminer)
ggcoxzph(test.ph1)

m1=coxph(Surv(dati$DFS,dati$relapse) ~ dati$ascite*factor(BRCA_dic)+dati$age, data=dati)
summary(m1)
test.ph1 <- cox.zph(m1); test.ph1
library(survminer)
ggcoxzph(test.ph1)


#########################
#chir_pi_lps_score#
#########################

#MORTALITA'
by(dati$chir_pi,dati$status_alive_dead,summary)
m2=coxph(Surv(dati$OS_diagn, dati$status_alive_dead) ~ dati$chir_pi, data=dati)
summary(m2)
test.ph2 <- cox.zph(m2); test.ph2
library(survminer)
ggcoxzph(test.ph2)

m2=coxph(Surv(dati$OS_diagn, dati$status_alive_dead) ~ dati$chir_pi*factor(BRCA_dic), data=dati)
summary(m2)
test.ph2 <- cox.zph(m2); test.ph2
library(survminer)
ggcoxzph(test.ph2)

m2=coxph(Surv(dati$OS_diagn, dati$status_alive_dead) ~ dati$chir_pi*factor(BRCA_dic)+dati$age, data=dati)
summary(m2)
test.ph2 <- cox.zph(m2); test.ph2
library(survminer)
ggcoxzph(test.ph2)


#RECIDIVA DFS
by(dati$chir_pi,dati$relapse,summary)
m1=coxph(Surv(dati$DFS,dati$relapse) ~ dati$chir_pi, data=dati)
summary(m1)
test.ph1 <- cox.zph(m1); test.ph1
library(survminer)
ggcoxzph(test.ph1)

m1=coxph(Surv(dati$DFS,dati$relapse) ~ dati$chir_pi*factor(BRCA_dic), data=dati)
summary(m1)
test.ph1 <- cox.zph(m1); test.ph1
library(survminer)
ggcoxzph(test.ph1)

m1=coxph(Surv(dati$DFS,dati$relapse) ~ dati$chir_pi*factor(BRCA_dic)+dati$age, data=dati)
summary(m1)
test.ph1 <- cox.zph(m1); test.ph1
library(survminer)
ggcoxzph(test.ph1)


#########################
#BMI####################
#########################

#MORTALITA'
by(dati$BMI,dati$status_alive_dead,summary)
m2=coxph(Surv(dati$OS_diagn, dati$status_alive_dead) ~ dati$BMI, data=dati)
summary(m2)
test.ph2 <- cox.zph(m2); test.ph2
library(survminer)
ggcoxzph(test.ph2)

m2=coxph(Surv(dati$OS_diagn, dati$status_alive_dead) ~ dati$BMI*factor(BRCA_dic), data=dati)
summary(m2)
test.ph2 <- cox.zph(m2); test.ph2
library(survminer)
ggcoxzph(test.ph2)

m2=coxph(Surv(dati$OS_diagn, dati$status_alive_dead) ~ dati$BMI*factor(BRCA_dic)+dati$age, data=dati)
summary(m2)
test.ph2 <- cox.zph(m2); test.ph2
library(survminer)
ggcoxzph(test.ph2)


#RECIDIVA DFS
by(dati$BMI,dati$relapse,summary)
m1=coxph(Surv(dati$DFS,dati$relapse) ~ dati$BMI, data=dati)
summary(m1)
test.ph1 <- cox.zph(m1); test.ph1
library(survminer)
ggcoxzph(test.ph1)

m1=coxph(Surv(dati$DFS,dati$relapse) ~ dati$BMI*factor(BRCA_dic), data=dati)
summary(m1)
test.ph1 <- cox.zph(m1); test.ph1
library(survminer)
ggcoxzph(test.ph1)

m1=coxph(Surv(dati$DFS,dati$relapse) ~ dati$BMI*factor(BRCA_dic)+dati$age, data=dati)
summary(m1)
test.ph1 <- cox.zph(m1); test.ph1
library(survminer)
ggcoxzph(test.ph1)


#########################
#menopausa##############
#########################

#MORTALITA'
by(dati$menopausa,dati$status_alive_dead,summary)
m2=coxph(Surv(dati$OS_diagn, dati$status_alive_dead) ~ dati$menopausa, data=dati)
summary(m2)
test.ph2 <- cox.zph(m2); test.ph2
library(survminer)
ggcoxzph(test.ph2)

m2=coxph(Surv(dati$OS_diagn, dati$status_alive_dead) ~ dati$menopausa*factor(BRCA_dic), data=dati)
summary(m2)
test.ph2 <- cox.zph(m2); test.ph2
library(survminer)
ggcoxzph(test.ph2)

m2=coxph(Surv(dati$OS_diagn, dati$status_alive_dead) ~ dati$menopausa*factor(BRCA_dic)+dati$age, data=dati)
summary(m2)
test.ph2 <- cox.zph(m2); test.ph2
library(survminer)
ggcoxzph(test.ph2)


#RECIDIVA DFS
by(dati$menopausa,dati$relapse,summary)
m1=coxph(Surv(dati$DFS,dati$relapse) ~ dati$menopausa, data=dati)
summary(m1)
test.ph1 <- cox.zph(m1); test.ph1
library(survminer)
ggcoxzph(test.ph1)

m1=coxph(Surv(dati$DFS,dati$relapse) ~ dati$menopausa*factor(BRCA_dic), data=dati)
summary(m1)
test.ph1 <- cox.zph(m1); test.ph1
library(survminer)
ggcoxzph(test.ph1)

m1=coxph(Surv(dati$DFS,dati$relapse) ~ dati$menopausa*factor(BRCA_dic)+dati$age, data=dati)
summary(m1)
test.ph1 <- cox.zph(m1); test.ph1
library(survminer)
ggcoxzph(test.ph1)


#########################
#ct_beva#
#########################

#MORTALITA'
by(dati$ct_Beva,dati$status_alive_dead,summary)
m2=coxph(Surv(dati$OS_diagn, dati$status_alive_dead) ~ dati$ct_Beva, data=dati)
summary(m2)
test.ph2 <- cox.zph(m2); test.ph2
library(survminer)
ggcoxzph(test.ph2)

m2=coxph(Surv(dati$OS_diagn, dati$status_alive_dead) ~ dati$ct_Beva*factor(BRCA_dic), data=dati)
summary(m2)
test.ph2 <- cox.zph(m2); test.ph2
library(survminer)
ggcoxzph(test.ph2)

m2=coxph(Surv(dati$OS_diagn, dati$status_alive_dead) ~ dati$ct_Beva*factor(BRCA_dic)+dati$age, data=dati)
summary(m2)
test.ph2 <- cox.zph(m2); test.ph2
library(survminer)
ggcoxzph(test.ph2)


#RECIDIVA DFS
by(dati$ct_Beva,dati$relapse,summary)
m1=coxph(Surv(dati$DFS,dati$relapse) ~ dati$ct_Beva, data=dati)
summary(m1)
test.ph1 <- cox.zph(m1); test.ph1
library(survminer)
ggcoxzph(test.ph1)

m1=coxph(Surv(dati$DFS,dati$relapse) ~ dati$ct_Beva*factor(BRCA_dic), data=dati)
summary(m1)
test.ph1 <- cox.zph(m1); test.ph1
library(survminer)
ggcoxzph(test.ph1)

m1=coxph(Surv(dati$DFS,dati$relapse) ~ dati$ct_Beva*factor(BRCA_dic)+dati$age, data=dati)
summary(m1)
test.ph1 <- cox.zph(m1); test.ph1
library(survminer)
ggcoxzph(test.ph1)


#################################
#citorid_cm#CITORIDOTTE SI/NO####
#################################

#MORTALITA'
by(dati$citorid_cm,dati$status_alive_dead,summary)
m2=coxph(Surv(dati$OS_diagn, dati$status_alive_dead) ~ dati$citorid_cm, data=dati)
summary(m2)
test.ph2 <- cox.zph(m2); test.ph2
library(survminer)
ggcoxzph(test.ph2)

m2=coxph(Surv(dati$OS_diagn, dati$status_alive_dead) ~ dati$citorid_cm*factor(BRCA_dic), data=dati)
summary(m2)
test.ph2 <- cox.zph(m2); test.ph2
library(survminer)
ggcoxzph(test.ph2)

m2=coxph(Surv(dati$OS_diagn, dati$status_alive_dead) ~ dati$citorid_cm*factor(BRCA_dic)+dati$age, data=dati)
summary(m2)
test.ph2 <- cox.zph(m2); test.ph2
library(survminer)
ggcoxzph(test.ph2)



#RECIDIVA DFS
by(dati$citorid_cm,dati$relapse,summary)
m1=coxph(Surv(dati$DFS,dati$relapse) ~ dati$citorid_cm, data=dati)
summary(m1)
test.ph1 <- cox.zph(m1); test.ph1
library(survminer)
ggcoxzph(test.ph1)

m1=coxph(Surv(dati$DFS,dati$relapse) ~ dati$citorid_cm*factor(BRCA_dic), data=dati)
summary(m1)
test.ph1 <- cox.zph(m1); test.ph1
library(survminer)
ggcoxzph(test.ph1)

m1=coxph(Surv(dati$DFS,dati$relapse) ~ dati$citorid_cm*factor(BRCA_dic)+dati$age, data=dati)
summary(m1)
test.ph1 <- cox.zph(m1); test.ph1
library(survminer)
ggcoxzph(test.ph1)




####################################
#######analisi stratificata menopausa##########
######################################

m00=coxph(Surv(PFSultimoPLATINO,relapse) ~ menopausa, data=subset(dati,BRCA_dic==0))
summary(m00)

m01=coxph(Surv(PFSultimoPLATINO,relapse) ~ menopausa, data=subset(dati,BRCA_dic==1))
summary(m01)
