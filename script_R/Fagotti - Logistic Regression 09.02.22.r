setwd("D:\\Progetti GSTeP\\Fagotti-Perrone TOC-1 - ongoing")
dati<-read_excel("TOC1.xlsx")
dim(dati)
colnames(dati)
head(dati)
str(dati)
summary(dati)

install.packages("survminer") #diagnostica su hazard proporzionali coeflo di Cox
install.packages("coxphw") #coeflo di Cox su hazard non proporzionali
install.packages("Hmisc") #matrice di correlazione (di Spearman)
install.packages("RVAideMemoire") #intervallo di confidenza su correlazione (di Sperman)
install.packages("spearmanCI") #intervallo di confidenza su correlazione (di Sperman)
install.packages("DescTools") #intervallo di confidenza su correlazione (di Sperman), mediante approssimazione Z
install.packages("My.stepwise")
install.packages("glmnet")

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


#################################
##########VALID DATASET##########
#################################

dati<-data.frame(ID=dati$"Paziente, cod sanit", status_alive_dead=dati$"Status (0=alive, 1=DOD)", relapse=dati$"Relapse or progression", OS_diagn=dati$"OS dalla diagnosi", PFSultimoPLATINO=dati$"PFS dall'ultimo platino", DFS=dati$"DFS dalla citoriduzione",BRCA_dic=dati$"BRCA1\\2_mut=1\\ WT=0",BRCA=dati$"BRCA wt_1_2", age=dati$Età,ARnucleo=dati$ARnucleo_Score_num,PR=dati$PR_Score_num,ERalpha=dati$ERAlpha_Score_num,ERbeta1nucleo=dati$ERBeta1_NucleoScore_num,ERbeta1cito=dati$ERBeta1_CitoScore_num,ERbeta2nucleo=dati$ERBeta2_NucleoScore_num,ERbeta2cito=dati$ERBeta2_CitoScore_num,ERbeta5nucleo=dati$ERBeta5_NucleoScore_num,ERbeta5cito=dati$ERBeta5_CitoScore_num, a_b1n=dati$"alpha_beta1n", a_b2n=dati$"alpha_beta2n", a_b5n=dati$"alpha_beta5n",  somatico=dati$"somatico codificato", PLARES=dati$"Platinum-resistance (0:no; 1:si)",BMI=dati$bmi_at_diagnosis,menopausa=dati$menopausa_at_diagnosis,ascite=dati$chir1_ascite,chir_pi=dati$chir1_pi_lps_score,ct_Beva=dati$ct_beva,citorid_cm=dati$"Citoridotte =1 , Non citoridotte = 0",citorid_cm_num=dati$"Cytored_RT(in cm)")

dim(dati)
summary(dati)
BRCA=dati$BRCA; length(BRCA)
BRCA_dic=dati$BRCA_dic




#####################################
###MODELLO LOGISTICO BRCA####
#####################################


################################
#########univariata####
################################

m0<-glm(formula = dati$PLARES ~ dati$age, data = dati, family = binomial(link = logit))
summary(m0)
exp(confint.default(m0))
exp(coef(m0))

m1<-glm(formula = dati$PLARES ~ dati$BRCA, data = dati, family = binomial(link = logit))
summary(m1)
exp(confint.default(m1))
exp(coef(m1))

m2<-glm(formula = dati$PLARES ~ dati$BMI, data = dati, family = binomial(link = logit))
summary(m2)
exp(confint.default(m2))
exp(coef(m2))

m3<-glm(formula = dati$PLARES ~ dati$menopausa, data = dati, family = binomial(link = logit))
summary(m3)
exp(confint.default(m3))
exp(coef(m3))

m4<-glm(formula = dati$PLARES ~ dati$ARnucleo, data = dati, family = binomial(link = logit))
summary(m4)
exp(confint.default(m4))
exp(coef(m4))

m5<-glm(formula = dati$PLARES ~ dati$PR, data = dati, family = binomial(link = logit))
summary(m5)
exp(confint.default(m5))
exp(coef(m5))

m6<-glm(formula = dati$PLARES ~ dati$ERalpha, data = dati, family = binomial(link = logit))
summary(m6)
exp(confint.default(m6))
exp(coef(m6))

m7<-glm(formula = dati$PLARES ~ dati$ERbeta1nucleo, data = dati, family = binomial(link = logit))
summary(m7)
exp(confint.default(m7))
exp(coef(m7))

m8<-glm(formula = dati$PLARES ~ dati$ERbeta1cito, data = dati, family = binomial(link = logit))
summary(m8)
exp(confint.default(m8))
exp(coef(m8))

m9<-glm(formula = dati$PLARES ~ dati$ERbeta2nucleo, data = dati, family = binomial(link = logit))
summary(m9)
exp(confint.default(m9))
exp(coef(m9))

m10<-glm(formula = dati$PLARES ~ dati$ERbeta2cito, data = dati, family = binomial(link = logit))
summary(m10)
exp(confint.default(m10))
exp(coef(m10))

m11<-glm(formula = dati$PLARES ~ dati$ERbeta5nucleo, data = dati, family = binomial(link = logit))
summary(m11)
exp(confint.default(m11))
exp(coef(m11))

m12<-glm(formula = dati$PLARES ~ dati$ERbeta5cito, data = dati, family = binomial(link = logit))
summary(m12)
exp(confint.default(m12))
exp(coef(m12))

m13<-glm(formula = dati$PLARES ~ dati$ascite, data = dati, family = binomial(link = logit))
summary(m13)
exp(confint.default(m13))
exp(coef(m13))

m14<-glm(formula = dati$PLARES ~ dati$chir_pi, data = dati, family = binomial(link = logit))
summary(m14)
exp(confint.default(m14))
exp(coef(m14))

m15<-glm(formula = dati$PLARES ~ dati$ct_Beva, data = dati, family = binomial(link = logit))
summary(m15)
exp(confint.default(m15))
exp(coef(m15))

m16<-glm(formula = dati$PLARES ~ dati$citorid, data = dati, family = binomial(link = logit))
summary(m16)
exp(confint.default(m16))
exp(coef(m16))

m17<-glm(formula = dati$PLARES ~ dati$a_b1n, data = dati, family = binomial(link = logit))
summary(m17)
exp(confint.default(m17))
exp(coef(m17))

m18<-glm(formula = dati$PLARES ~ dati$a_b2n, data = dati, family = binomial(link = logit))
summary(m18)
exp(confint.default(m18))
exp(coef(m18))

m19<-glm(formula = dati$PLARES ~ dati$a_b5n, data = dati, family = binomial(link = logit))
summary(m19)
exp(confint.default(m19))
exp(coef(m19))



##########################
###interaction BRCA dic###
##########################

m0<-glm(formula = dati$PLARES ~ dati$age*factor(BRCA_dic), data = dati, family = binomial(link = logit))
summary(m0)
exp(confint.default(m0)))
exp(coef(m0))

m2<-glm(formula = dati$PLARES ~ dati$BMI*factor(BRCA_dic), data = dati, family = binomial(link = logit))
summary(m2)
exp(coef(m2))
exp(confint.default(m2))

m3<-glm(formula = dati$PLARES ~ dati$menopausa*factor(BRCA_dic), data = dati, family = binomial(link = logit))
summary(m3)
exp(coef(m3))
exp(confint.default(m3))

m4<-glm(formula = dati$PLARES ~ dati$ARnucleo*factor(BRCA_dic), data = dati, family = binomial(link = logit))
summary(m4)
exp(coef(m4))
exp(confint.default(m4))

m5<-glm(formula = dati$PLARES ~ dati$PR*factor(BRCA_dic), data = dati, family = binomial(link = logit))
summary(m5)
exp(coef(m5))
exp(confint.default(m5))

m6<-glm(formula = dati$PLARES ~ dati$ERalpha*factor(BRCA_dic), data = dati, family = binomial(link = logit))
summary(m6)
exp(coef(m6))
exp(confint.default(m6))

m7<-glm(formula = dati$PLARES ~ dati$ERbeta1nucleo*factor(BRCA_dic), data = dati, family = binomial(link = logit))
summary(m7)
exp(coef(m7))
exp(confint.default(m7))

m8<-glm(formula = dati$PLARES ~ dati$ERbeta1cito*factor(BRCA_dic), data = dati, family = binomial(link = logit))
summary(m8)
exp(coef(m8))
exp(confint.default(m8))

m9<-glm(formula = dati$PLARES ~ dati$ERbeta2nucleo*factor(BRCA_dic), data = dati, family = binomial(link = logit))
summary(m9)
exp(coef(m9))
exp(confint.default(m9))

m10<-glm(formula = dati$PLARES ~ dati$ERbeta2cito*factor(BRCA_dic), data = dati, family = binomial(link = logit))
summary(m10)
exp(confint.default(m10))
exp(coef(m10))

m11<-glm(formula = dati$PLARES ~ dati$ERbeta5nucleo*factor(BRCA_dic), data = dati, family = binomial(link = logit))
summary(m11)
exp(confint.default(m11))
exp(coef(m11))

m12<-glm(formula = dati$PLARES ~ dati$ERbeta5cito*factor(BRCA_dic), data = dati, family = binomial(link = logit))
summary(m12)
exp(confint.default(m12))
exp(coef(m12))

m13<-glm(formula = dati$PLARES ~ dati$ascite*factor(BRCA_dic), data = dati, family = binomial(link = logit))
summary(m13)
exp(confint.default(m13))
exp(coef(m13))

m14<-glm(formula = dati$PLARES ~ dati$chir_pi*factor(BRCA_dic), data = dati, family = binomial(link = logit))
summary(m14)
exp(confint.default(m14))
exp(coef(m14))

m15<-glm(formula = dati$PLARES ~ dati$ct_Beva*factor(BRCA_dic), data = dati, family = binomial(link = logit))
summary(m15)
exp(confint.default(m15))
exp(coef(m15))

m16<-glm(formula = dati$PLARES ~ dati$citorid_cm*factor(BRCA_dic), data = dati, family = binomial(link = logit))
summary(m16)
exp(confint.default(m16))
exp(coef(m16))

m17<-glm(formula = dati$PLARES ~ dati$a_b1n*factor(BRCA_dic), data = dati, family = binomial(link = logit))
summary(m17)
exp(confint.default(m17))
exp(coef(m17))

m18<-glm(formula = dati$PLARES ~ dati$a_b2n*factor(BRCA_dic), data = dati, family = binomial(link = logit))
summary(m18)
exp(confint.default(m18))
exp(coef(m18))

m19<-glm(formula = dati$PLARES ~ dati$a_b5n*factor(BRCA_dic), data = dati, family = binomial(link = logit))
summary(m19)
exp(confint.default(m19))
exp(coef(m19))

##########################
###interaction BRCA###
##########################

m0<-glm(formula = dati$PLARES ~ dati$age*factor(BRCA), data = dati, family = binomial(link = logit))
summary(m0)
exp(confint.default(m0))
exp(coef(m0))

m3<-glm(formula = dati$PLARES ~ dati$menopausa*factor(BRCA), data = dati, family = binomial(link = logit))
summary(m3)
exp(confint.default(m3))
exp(coef(m3))

m4<-glm(formula = dati$PLARES ~ dati$ARnucleo*factor(BRCA), data = dati, family = binomial(link = logit))
summary(m4)
exp(confint.default(m4))
exp(coef(m4))

m5<-glm(formula = dati$PLARES ~ dati$PR*factor(BRCA), data = dati, family = binomial(link = logit))
summary(m5)
exp(confint.default(m5))
exp(coef(m5))

m6<-glm(formula = dati$PLARES ~ dati$ERalpha*factor(BRCA), data = dati, family = binomial(link = logit))
summary(m6)
exp(confint.default(m6))
exp(coef(m6))

m7<-glm(formula = dati$PLARES ~ dati$ERbeta1nucleo*factor(BRCA), data = dati, family = binomial(link = logit))
summary(m7)
exp(confint.default(m7))
exp(coef(m7))

m8<-glm(formula = dati$PLARES ~ dati$ERbeta1cito*factor(BRCA), data = dati, family = binomial(link = logit))
summary(m8)
exp(confint.default(m8))
exp(coef(m8))

m9<-glm(formula = dati$PLARES ~ dati$ERbeta2nucleo*factor(BRCA), data = dati, family = binomial(link = logit))
summary(m9)
exp(confint.default(m9))
exp(coef(m9))

m10<-glm(formula = dati$PLARES ~ dati$ERbeta2cito*factor(BRCA), data = dati, family = binomial(link = logit))
summary(m10)
exp(confint.default(m10))
exp(coef(m10))

m11<-glm(formula = dati$PLARES ~ dati$ERbeta5nucleo*factor(BRCA), data = dati, family = binomial(link = logit))
summary(m11)
exp(confint.default(m11))
exp(coef(m11))

m12<-glm(formula = dati$PLARES ~ dati$ERbeta5cito*factor(BRCA), data = dati, family = binomial(link = logit))
summary(m12)
exp(confint.default(m12))
exp(coef(m12))

m13<-glm(formula = dati$PLARES ~ dati$ascite*factor(BRCA), data = dati, family = binomial(link = logit))
summary(m13)
exp(confint.default(m13))
exp(coef(m13))

m14<-glm(formula = dati$PLARES ~ dati$chir_pi*factor(BRCA), data = dati, family = binomial(link = logit))
summary(m14)
exp(confint.default(m14))
exp(coef(m14))

m15<-glm(formula = dati$PLARES ~ dati$ct_Beva*factor(BRCA), data = dati, family = binomial(link = logit))
summary(m15)
exp(confint.default(m15))
exp(coef(m15))

m16<-glm(formula = dati$PLARES ~ dati$citorid_cm*factor(BRCA), data = dati, family = binomial(link = logit))
summary(m16)
exp(confint.default(m16))
exp(coef(m16))

m17<-glm(formula = dati$PLARES ~ dati$a_b1n*factor(BRCA), data = dati, family = binomial(link = logit))
summary(m17)
exp(confint.default(m17))
exp(coef(m17))

m18<-glm(formula = dati$PLARES ~ dati$a_b2n*factor(BRCA), data = dati, family = binomial(link = logit))
summary(m18)
exp(confint.default(m18))
exp(coef(m18))

m19<-glm(formula = dati$PLARES ~ dati$a_b5n*factor(BRCA), data = dati, family = binomial(link = logit))
summary(m19)
exp(confint.default(m19))
exp(coef(m19))


########################################
###interaction BRCA dic age adjusted###
#######################################

m2<-glm(formula = dati$PLARES ~ dati$BMI*factor(BRCA_dic) + dati$age, data = dati, family = binomial(link = logit))
summary(m2)
exp(coef(m2))
exp(confint.default(m2))

m3<-glm(formula = dati$PLARES ~ dati$menopausa*factor(BRCA_dic) + dati$age, data = dati, family = binomial(link = logit))
summary(m3)
exp(coef(m3))
exp(confint.default(m3))

m4<-glm(formula = dati$PLARES ~ dati$ARnucleo*factor(BRCA_dic) + dati$age, data = dati, family = binomial(link = logit))
summary(m4)
exp(coef(m4))
exp(confint.default(m4))

m5<-glm(formula = dati$PLARES ~ dati$PR*factor(BRCA_dic) + dati$age, data = dati, family = binomial(link = logit))
summary(m5)
exp(coef(m5))
exp(confint.default(m5))

m6<-glm(formula = dati$PLARES ~ dati$ERalpha*factor(BRCA_dic) + dati$age, data = dati, family = binomial(link = logit))
summary(m6)
exp(coef(m6))
exp(confint.default(m6))

m7<-glm(formula = dati$PLARES ~ dati$ERbeta1nucleo*factor(BRCA_dic) + dati$age, data = dati, family = binomial(link = logit))
summary(m7)
exp(coef(m7))
exp(confint.default(m7))

m8<-glm(formula = dati$PLARES ~ dati$ERbeta1cito*factor(BRCA_dic) + dati$age, data = dati, family = binomial(link = logit))
summary(m8)
exp(coef(m8))
exp(confint.default(m8))

m9<-glm(formula = dati$PLARES ~ dati$ERbeta2nucleo*factor(BRCA_dic) + dati$age, data = dati, family = binomial(link = logit))
summary(m9)
exp(coef(m9))
exp(confint.default(m9))

m10<-glm(formula = dati$PLARES ~ dati$ERbeta2cito*factor(BRCA_dic) + dati$age, data = dati, family = binomial(link = logit))
summary(m10)
exp(confint.default(m10))
exp(coef(m10))

m11<-glm(formula = dati$PLARES ~ dati$ERbeta5nucleo*factor(BRCA_dic) + dati$age, data = dati, family = binomial(link = logit))
summary(m11)
exp(confint.default(m11))
exp(coef(m11))

m12<-glm(formula = dati$PLARES ~ dati$ERbeta5cito*factor(BRCA_dic) + dati$age, data = dati, family = binomial(link = logit))
summary(m12)
exp(confint.default(m12))
exp(coef(m12))

m13<-glm(formula = dati$PLARES ~ dati$ascite*factor(BRCA_dic) + dati$age, data = dati, family = binomial(link = logit))
summary(m13)
exp(confint.default(m13))
exp(coef(m13))

m14<-glm(formula = dati$PLARES ~ dati$chir_pi*factor(BRCA_dic) + dati$age, data = dati, family = binomial(link = logit))
summary(m14)
exp(confint.default(m14))
exp(coef(m14))

m15<-glm(formula = dati$PLARES ~ dati$ct_Beva*factor(BRCA_dic) + dati$age, data = dati, family = binomial(link = logit))
summary(m15)
exp(confint.default(m15))
exp(coef(m15))

m16<-glm(formula = dati$PLARES ~ dati$citorid_cm*factor(BRCA_dic) + dati$age, data = dati, family = binomial(link = logit))
summary(m16)
exp(confint.default(m16))
exp(coef(m16))

m17<-glm(formula = dati$PLARES ~ dati$a_b1n*factor(BRCA_dic) + dati$age, data = dati, family = binomial(link = logit))
summary(m17)
exp(confint.default(m17))
exp(coef(m17))

m18<-glm(formula = dati$PLARES ~ dati$a_b2n*factor(BRCA_dic) + dati$age, data = dati, family = binomial(link = logit))
summary(m18)
exp(confint.default(m18))
exp(coef(m18))

m19<-glm(formula = dati$PLARES ~ dati$a_b5n*factor(BRCA_dic) + dati$age, data = dati, family = binomial(link = logit))
summary(m19)
exp(confint.default(m19))
exp(coef(m19))

###################################
###interaction BRCA age adjusted###
##################################


m2<-glm(formula = dati$PLARES ~ dati$BMI*factor(BRCA) + dati$age, data = dati, family = binomial(link = logit))
summary(m2)
exp(confint.default(m2))
exp(coef(m2))

m3<-glm(formula = dati$PLARES ~ dati$menopausa*factor(BRCA) + dati$age, data = dati, family = binomial(link = logit))
summary(m3)
exp(confint.default(m3))
exp(coef(m3))

m4<-glm(formula = dati$PLARES ~ dati$ARnucleo*factor(BRCA) + dati$age, data = dati, family = binomial(link = logit))
summary(m4)
exp(confint.default(m4))
exp(coef(m4))

m5<-glm(formula = dati$PLARES ~ dati$PR*factor(BRCA) + dati$age, data = dati, family = binomial(link = logit))
summary(m5)
exp(confint.default(m5))
exp(coef(m5))

m6<-glm(formula = dati$PLARES ~ dati$ERalpha*factor(BRCA) + dati$age, data = dati, family = binomial(link = logit))
summary(m6)
exp(confint.default(m6))
exp(coef(m6))

m7<-glm(formula = dati$PLARES ~ dati$ERbeta1nucleo*factor(BRCA) + dati$age, data = dati, family = binomial(link = logit))
summary(m7)
exp(confint.default(m7))
exp(coef(m7))

m8<-glm(formula = dati$PLARES ~ dati$ERbeta1cito*factor(BRCA) + dati$age, data = dati, family = binomial(link = logit))
summary(m8)
exp(confint.default(m8))
exp(coef(m8))

m9<-glm(formula = dati$PLARES ~ dati$ERbeta2nucleo*factor(BRCA) + dati$age, data = dati, family = binomial(link = logit))
summary(m9)
exp(confint.default(m9))
exp(coef(m9))

m10<-glm(formula = dati$PLARES ~ dati$ERbeta2cito*factor(BRCA) + dati$age, data = dati, family = binomial(link = logit))
summary(m10)
exp(confint.default(m10))
exp(coef(m10))

m11<-glm(formula = dati$PLARES ~ dati$ERbeta5nucleo*factor(BRCA) + dati$age, data = dati, family = binomial(link = logit))
summary(m11)
exp(confint.default(m11))
exp(coef(m11))

m12<-glm(formula = dati$PLARES ~ dati$ERbeta5cito*factor(BRCA) + dati$age, data = dati, family = binomial(link = logit))
summary(m12)
exp(confint.default(m12))
exp(coef(m12))

m13<-glm(formula = dati$PLARES ~ dati$ascite*factor(BRCA) + dati$age, data = dati, family = binomial(link = logit))
summary(m13)
exp(confint.default(m13))
exp(coef(m13))

m14<-glm(formula = dati$PLARES ~ dati$chir_pi*factor(BRCA) + dati$age, data = dati, family = binomial(link = logit))
summary(m14)
exp(confint.default(m14))
exp(coef(m14))

m15<-glm(formula = dati$PLARES ~ dati$ct_Beva*factor(BRCA) + dati$age, data = dati, family = binomial(link = logit))
summary(m15)
exp(confint.default(m15))
exp(coef(m15))

m16<-glm(formula = dati$PLARES ~ dati$citorid_cm*factor(BRCA) + dati$age, data = dati, family = binomial(link = logit))
summary(m16)
exp(confint.default(m16))
exp(coef(m16))

m17<-glm(formula = dati$PLARES ~ dati$a_b1n*factor(BRCA) + dati$age, data = dati, family = binomial(link = logit))
summary(m17)
exp(confint.default(m17))
exp(coef(m17))

m18<-glm(formula = dati$PLARES ~ dati$a_b2n*factor(BRCA) + dati$age, data = dati, family = binomial(link = logit))
summary(m18)
exp(confint.default(m18))
exp(coef(m18))

m19<-glm(formula = dati$PLARES ~ dati$a_b5n*factor(BRCA) + dati$age, data = dati, family = binomial(link = logit))
summary(m19)
exp(confint.default(m19))
exp(coef(m19))





