##########processing the data from Friedrich2021########
rm(list=ls())
library(readr)

###import μ
μ_Friedrich2021<- read_csv("Data/growth_rate/growth_rate.csv")#h-1
colnames(μ_Friedrich2021)<-c("condition","μ")


###import phi0
phi0<-0.08 


###import the calculation function of effective mR
source("cal_eff_mR.R")


###import ki
ki <- read_csv("Data/elongation_speed/ki.csv")
colnames(ki)<-c("geneid","ki(da/s)")

###import xi
ribo_wt_glu <- read_csv("Data/riboseq/wt-glu.csv")
colnames(ribo_wt_glu)<-c("geneid","Xi")

ribo_da_glu <- read_csv("Data/riboseq/da-glu.csv")  
colnames(ribo_da_glu)<-c("geneid","Xi")

ribo_wt_gly <- read_csv("Data/riboseq/wt-gly.csv")
colnames(ribo_wt_gly)<-c("geneid","Xi")

ribo_da_gly <- read_csv("Data/riboseq/da-gly.csv") 
colnames(ribo_da_gly)<-c("geneid","Xi")

###import proteomics data
pro_wt_glu <- read_csv("Data/proteomics/data_L_calib_057/wt_glu.csv")
colnames(pro_wt_glu)<-c("geneid","phi")

pro_da_glu <- read_csv("Data/proteomics/data_L_calib_057/da_glu.csv")  
colnames(pro_da_glu)<-c("geneid","phi")

pro_wt_gly <- read_csv("Data/proteomics/data_L_calib_057/wt_gly.csv")
colnames(pro_wt_gly)<-c("geneid","phi")

pro_da_gly <- read_csv("Data/proteomics/data_L_calib_057/da_gly.csv") 
colnames(pro_da_gly)<-c("geneid","phi")

###import proteomics data without calibration of phi for Figure S2b
#pro_wt_glu <- read_csv("Data/proteomics/data_no_calib/wt_glu.csv")
#colnames(pro_wt_glu)<-c("geneid","phi")

#pro_da_glu <- read_csv("Data/proteomics/data_no_calib/da_glu.csv")  
#colnames(pro_da_glu)<-c("geneid","phi")

#pro_wt_gly <- read_csv("Data/proteomics/data_no_calib/wt_gly.csv")
#colnames(pro_wt_gly)<-c("geneid","phi")

#pro_da_gly <- read_csv("Data/proteomics/data_no_calib/da_gly.csv") 
#colnames(pro_da_gly)<-c("geneid","phi")

###import proteomics data in which phi is calibrated with L−1 for Figure S2c
#pro_wt_glu <- read_csv("Data/proteomics/data_L_calib_1/wt_glu.csv")
#colnames(pro_wt_glu)<-c("geneid","phi")

#pro_da_glu <- read_csv("Data/proteomics/data_L_calib_1/da_glu.csv")  
#colnames(pro_da_glu)<-c("geneid","phi")

#pro_wt_gly <- read_csv("Data/proteomics/data_L_calib_1/wt_gly.csv")
#colnames(pro_wt_gly)<-c("geneid","phi")

#pro_da_gly <- read_csv("Data/proteomics/data_L_calib_1/da_gly.csv") 
#colnames(pro_da_gly)<-c("geneid","phi")

###import alpha
a_1 <- read_csv("Data/degradation_rate/a_1.csv")


###import ribosome protein list and mass
ribopro <- read_csv("Data/ribosomal_protein/ribopro_eu.csv")
ribomass <- read_csv("Data/ribosomal_protein/ribo_eu_sub_mass.csv")



##########process the data
sum_data_Friedrich2021<-data.frame()
cor_phipre_chi<-data.frame()

####wt-glu
temp1<-merge(pro_wt_glu,ribopro,by="geneid")
phR<-sum(temp1$phi)
mR<-cal_eff_mR(as.vector(temp1$geneid),ribomass)

temp1<-merge(pro_wt_glu,ribo_wt_glu,by="geneid")
temp1<-merge(temp1,ki,by="geneid")
temp1<-merge(temp1,a_1,by="geneid")

temp1$phi<-temp1$phi/sum(temp1$phi)
temp1$Xi<-temp1$Xi/sum(temp1$Xi)


temp2<-merge(temp1,ribopro,by="geneid")
kR<-sum(temp2$`ki(da/s)`*temp2$Xi)/sum(temp2$Xi)
aR<-sum(temp2$a_h*temp2$phi)/sum(temp2$phi)


temp3<-setdiff(temp1$geneid,ribopro$geneid)
temp3<-data.frame(geneid=temp3)
temp4<-merge(temp3,temp1,by="geneid")
mean_a_oth<-mean(temp4$a_h)
mean_k_oth<-mean(temp4$`ki(da/s)`)


Ikx<-cor(temp4$`ki(da/s)`,temp4$Xi)*(sd(temp4$`ki(da/s)`)/mean(temp4$`ki(da/s)`))*(sd(temp4$Xi/mean(temp4$Xi)))
Iax_phi<-cor(temp4$a_h,temp4$phi)*(sd(temp4$a_h)/mean(temp4$a_h))*(sd(temp4$phi/mean(temp4$phi)))

temp1<-c(Ikx,Iax_phi,mR,phR,as.numeric(μ_Friedrich2021[1,2]),"wt-glu")
sum_data_Friedrich2021<-rbind(sum_data_Friedrich2021,temp1,stringsAsFactors = F)



##calculate the correlation between phi_pre and other parameters
#phi_pre
temp2<-temp4$`ki(da/s)`*temp4$Xi/((as.numeric(μ_Friedrich2021[1,2])+temp4$a_h)/3600)
temp3<-temp2/sum(temp2)
temp5<-cor(temp4$phi,temp4$Xi) #phi and chi
temp6<-cor(temp3,temp4$Xi)     #phi_pre and chi
temp7<-cor(temp3,temp4$phi)    #phi_pre and phi
#the case in which alpha=0
temp8<-temp4$`ki(da/s)`*temp4$Xi/((as.numeric(μ_Friedrich2021[1,2]))/3600)
temp9<-temp8/sum(temp8)
temp10<-cor(temp9,temp4$Xi)    #phi_pre and chi
temp11<-cor(temp9,temp4$phi)   #phi_pre and phi
#the case in which k=k_ave
temp8<-mean_k_oth*temp4$Xi/((as.numeric(μ_Friedrich2021[1,2])+temp4$a_h)/3600)
temp9<-temp8/sum(temp8)
temp12<-cor(temp9,temp4$Xi)    #phi_pre and chi
temp13<-cor(temp9,temp4$phi)   #phi_pre and phi
#the case in which k=k_ave and alpha=0
temp8<-mean_k_oth*temp4$Xi/((as.numeric(μ_Friedrich2021[1,2]))/3600)
temp9<-temp8/sum(temp8)
temp14<-cor(temp9,temp4$Xi)    #phi_pre and chi
temp15<-cor(temp9,temp4$phi)   #phi_pre and phi
cor_phipre_chi<-rbind(cor_phipre_chi,c(temp5,temp6,temp7,temp10,temp11,temp12,temp13,temp14,temp15))


####wt-gly
temp1<-merge(pro_wt_gly,ribopro,by="geneid")
phR<-sum(temp1$phi)
mR<-cal_eff_mR(as.vector(temp1$geneid),ribomass)

temp1<-merge(pro_wt_gly,ribo_wt_gly,by="geneid")
temp1<-merge(temp1,ki,by="geneid")
temp1<-merge(temp1,a_1,by="geneid")

temp1$phi<-temp1$phi/sum(temp1$phi)
temp1$Xi<-temp1$Xi/sum(temp1$Xi)

temp2<-merge(temp1,ribopro,by="geneid")
kR<-sum(temp2$`ki(da/s)`*temp2$Xi)/sum(temp2$Xi)
aR<-sum(temp2$a_h*temp2$phi)/sum(temp2$phi)


temp3<-setdiff(temp1$geneid,ribopro$geneid)
temp3<-data.frame(geneid=temp3)
temp4<-merge(temp3,temp1,by="geneid")
mean_a_oth<-mean(temp4$a_h)
mean_k_oth<-mean(temp4$`ki(da/s)`)

Ikx<-cor(temp4$`ki(da/s)`,temp4$Xi)*(sd(temp4$`ki(da/s)`)/mean(temp4$`ki(da/s)`))*(sd(temp4$Xi/mean(temp4$Xi)))
Iax_phi<-cor(temp4$a_h,temp4$phi)*(sd(temp4$a_h)/mean(temp4$a_h))*(sd(temp4$phi/mean(temp4$phi)))

temp1<-c(Ikx,Iax_phi,mR,phR,as.numeric(μ_Friedrich2021[5,2]),"wt-gly")
sum_data_Friedrich2021<-rbind(sum_data_Friedrich2021,temp1,stringsAsFactors = F)

##calculate the correlation between phi_pre and other parameters
#phi_pre
temp2<-temp4$`ki(da/s)`*temp4$Xi/((as.numeric(μ_Friedrich2021[5,2])+temp4$a_h)/3600)
temp3<-temp2/sum(temp2)
temp5<-cor(temp4$phi,temp4$Xi) #phi and chi
temp6<-cor(temp3,temp4$Xi)     #phi_pre and chi
temp7<-cor(temp3,temp4$phi)    #phi_pre and phi
#the case in which alpha=0
temp8<-temp4$`ki(da/s)`*temp4$Xi/((as.numeric(μ_Friedrich2021[5,2]))/3600)
temp9<-temp8/sum(temp8)
temp10<-cor(temp9,temp4$Xi)    #phi_pre and chi
temp11<-cor(temp9,temp4$phi)   #phi_pre and phi
#the case in which k=k_ave
temp8<-mean_k_oth*temp4$Xi/((as.numeric(μ_Friedrich2021[5,2])+temp4$a_h)/3600)
temp9<-temp8/sum(temp8)
temp12<-cor(temp9,temp4$Xi)    #phi_pre and chi
temp13<-cor(temp9,temp4$phi)   #phi_pre and phi
#the case in which k=k_ave and alpha=0
temp8<-mean_k_oth*temp4$Xi/((as.numeric(μ_Friedrich2021[5,2]))/3600)
temp9<-temp8/sum(temp8)
temp14<-cor(temp9,temp4$Xi)    #phi_pre and chi
temp15<-cor(temp9,temp4$phi)   #phi_pre and phi
cor_phipre_chi<-rbind(cor_phipre_chi,c(temp5,temp6,temp7,temp10,temp11,temp12,temp13,temp14,temp15))

####da-glu
temp1<-merge(pro_da_glu,ribopro,by="geneid")
phR<-sum(temp1$phi)
mR<-cal_eff_mR(as.vector(temp1$geneid),ribomass)

temp1<-merge(pro_da_glu,ribo_da_glu,by="geneid")
temp1<-merge(temp1,ki,by="geneid")
temp1<-merge(temp1,a_1,by="geneid")

temp1$phi<-temp1$phi/sum(temp1$phi)
temp1$Xi<-temp1$Xi/sum(temp1$Xi)

temp2<-merge(temp1,ribopro,by="geneid")
kR<-sum(temp2$`ki(da/s)`*temp2$Xi)/sum(temp2$Xi)
aR<-sum(temp2$a_h*temp2$phi)/sum(temp2$phi)


temp3<-setdiff(temp1$geneid,ribopro$geneid)
temp3<-data.frame(geneid=temp3)
temp4<-merge(temp3,temp1,by="geneid")
mean_a_oth<-mean(temp4$a_h)
mean_k_oth<-mean(temp4$`ki(da/s)`)

Ikx<-cor(temp4$`ki(da/s)`,temp4$Xi)*(sd(temp4$`ki(da/s)`)/mean(temp4$`ki(da/s)`))*(sd(temp4$Xi/mean(temp4$Xi)))
Iax_phi<-cor(temp4$a_h,temp4$phi)*(sd(temp4$a_h)/mean(temp4$a_h))*(sd(temp4$phi/mean(temp4$phi)))

temp1<-c(Ikx,Iax_phi,mR,phR,as.numeric(μ_Friedrich2021[2,2]),"da-glu")
sum_data_Friedrich2021<-rbind(sum_data_Friedrich2021,temp1,stringsAsFactors = F)

##calculate the correlation between phi_pre and other parameters
#phi_pre
temp2<-temp4$`ki(da/s)`*temp4$Xi/((as.numeric(μ_Friedrich2021[2,2])+temp4$a_h)/3600)
temp3<-temp2/sum(temp2)
temp5<-cor(temp4$phi,temp4$Xi) #phi and chi
temp6<-cor(temp3,temp4$Xi)     #phi_pre and chi
temp7<-cor(temp3,temp4$phi)    #phi_pre and phi
#the case in which alpha=0
temp8<-temp4$`ki(da/s)`*temp4$Xi/((as.numeric(μ_Friedrich2021[2,2]))/3600)
temp9<-temp8/sum(temp8)
temp10<-cor(temp9,temp4$Xi)    #phi_pre and chi
temp11<-cor(temp9,temp4$phi)   #phi_pre and phi
#the case in which k=k_ave
temp8<-mean_k_oth*temp4$Xi/((as.numeric(μ_Friedrich2021[2,2])+temp4$a_h)/3600)
temp9<-temp8/sum(temp8)
temp12<-cor(temp9,temp4$Xi)    #phi_pre and chi
temp13<-cor(temp9,temp4$phi)   #phi_pre and phi
#the case in which k=k_ave and alpha=0
temp8<-mean_k_oth*temp4$Xi/((as.numeric(μ_Friedrich2021[2,2]))/3600)
temp9<-temp8/sum(temp8)
temp14<-cor(temp9,temp4$Xi)    #phi_pre and chi
temp15<-cor(temp9,temp4$phi)   #phi_pre and phi
cor_phipre_chi<-rbind(cor_phipre_chi,c(temp5,temp6,temp7,temp10,temp11,temp12,temp13,temp14,temp15))

####da-gly
temp1<-merge(pro_da_gly,ribopro,by="geneid")
phR<-sum(temp1$phi)
mR<-cal_eff_mR(as.vector(temp1$geneid),ribomass)

temp1<-merge(pro_da_gly,ribo_da_gly,by="geneid")
temp1<-merge(temp1,ki,by="geneid")
temp1<-merge(temp1,a_1,by="geneid")

temp1$phi<-temp1$phi/sum(temp1$phi)
temp1$Xi<-temp1$Xi/sum(temp1$Xi)

temp2<-merge(temp1,ribopro,by="geneid")
kR<-sum(temp2$`ki(da/s)`*temp2$Xi)/sum(temp2$Xi)
aR<-sum(temp2$a_h*temp2$phi)/sum(temp2$phi)


temp3<-setdiff(temp1$geneid,ribopro$geneid)
temp3<-data.frame(geneid=temp3)
temp4<-merge(temp3,temp1,by="geneid")
mean_a_oth<-mean(temp4$a_h)
mean_k_oth<-mean(temp4$`ki(da/s)`)

Ikx<-cor(temp4$`ki(da/s)`,temp4$Xi)*(sd(temp4$`ki(da/s)`)/mean(temp4$`ki(da/s)`))*(sd(temp4$Xi/mean(temp4$Xi)))
Iax_phi<-cor(temp4$a_h,temp4$phi)*(sd(temp4$a_h)/mean(temp4$a_h))*(sd(temp4$phi/mean(temp4$phi)))

temp1<-c(Ikx,Iax_phi,mR,phR,as.numeric(μ_Friedrich2021[6,2]),"da-gly")
sum_data_Friedrich2021<-rbind(sum_data_Friedrich2021,temp1,stringsAsFactors = F)

##calculate the correlation between phi_pre and other parameters
#phi_pre
temp2<-temp4$`ki(da/s)`*temp4$Xi/((as.numeric(μ_Friedrich2021[6,2])+temp4$a_h)/3600)
temp3<-temp2/sum(temp2)
temp5<-cor(temp4$phi,temp4$Xi) #phi and chi
temp6<-cor(temp3,temp4$Xi)     #phi_pre and chi
temp7<-cor(temp3,temp4$phi)    #phi_pre and phi
#the case in which alpha=0
temp8<-temp4$`ki(da/s)`*temp4$Xi/((as.numeric(μ_Friedrich2021[6,2]))/3600)
temp9<-temp8/sum(temp8)
temp10<-cor(temp9,temp4$Xi)    #phi_pre and chi
temp11<-cor(temp9,temp4$phi)   #phi_pre and phi
#the case in which k=k_ave
temp8<-mean_k_oth*temp4$Xi/((as.numeric(μ_Friedrich2021[6,2])+temp4$a_h)/3600)
temp9<-temp8/sum(temp8)
temp12<-cor(temp9,temp4$Xi)    #phi_pre and chi
temp13<-cor(temp9,temp4$phi)   #phi_pre and phi
#the case in which k=k_ave and alpha=0
temp8<-mean_k_oth*temp4$Xi/((as.numeric(μ_Friedrich2021[6,2]))/3600)
temp9<-temp8/sum(temp8)
temp14<-cor(temp9,temp4$Xi)    #phi_pre and chi
temp15<-cor(temp9,temp4$phi)   #phi_pre and phi
cor_phipre_chi<-rbind(cor_phipre_chi,c(temp5,temp6,temp7,temp10,temp11,temp12,temp13,temp14,temp15))

#############################
colnames(sum_data_Friedrich2021)<-c("Ixk","Ipha","mR","phR","mu","comb")
sum_data_Friedrich2021$mu<-as.numeric(sum_data_Friedrich2021$mu)/60
colnames(sum_data_Friedrich2021)<-c("Ixk","Ipha","mR(Da)","phR","mu(1/min)","comb")

colnames(cor_phipre_chi)<-c("rho_phi_chi","rho_phipre_chi","rho_phipre_phi","rho_phipre_chi_a0","rho_phipre_phi_a0","rho_phipre_chi_keq","rho_phipre_phi_keq","rho_phipre_chi_a0keq","rho_phipre_phi_a0keq")


###The results for validation in matlab (Simu_theory_validation.mlx)
write.csv(sum_data_Friedrich2021,"sum_data_Friedrich2021.csv")
#write.csv(sum_data_Friedrich2021,"sum_data_Friedrich2021_nocalib.csv")
#write.csv(sum_data_Friedrich2021,"sum_data_Friedrich2021_Lcalib.csv")

###The results for Table S2
write.csv(cor_phipre_chi,"cor_phipre_chi.csv")

