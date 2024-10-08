library(readxl)

sheet_SCD<-excel_sheets("./SCD.xlsx")
sheet_SCG<-excel_sheets("./SCG.xlsx")


i=1
temp7<-vector()
for (i in 1:8) 
{
  if (i<5) 
  {
    temp1<- read_excel("./SCD.xlsx",sheet = i)    
  }else{
    temp1<- read_excel("./SCG.xlsx",sheet = i-4)  
  }
  
  temp1<-temp1[,-3]
  temp1$ln<-unlist(log(temp1[,2])[,1])
  colnames(temp1)<-c("time","OD","ln")
  
  j=1
  temp2<-vector()
  for (j in 1:(nrow(temp1)-4)) 
  {
    temp3<-lm(ln~time,temp1[j:(j+4),])
    temp2<-c(temp2,temp3[["coefficients"]][["time"]])
  }
  temp4<-max(temp2)
  temp5<-which(temp2>0.9*temp4)
  temp6<-lm(ln~time,temp1[min(temp5):(max(temp5)+4),])[["coefficients"]][["time"]]
  temp7<-c(temp7,temp6)
  
}
temp7<-as.data.frame(temp7)
rownames(temp7)<-c(paste0(sheet_SCD,"-glu"),paste0(sheet_SCG,"-gly"))
colnames(temp7)<-"growth_rate/h-1"
write.csv(temp7,file = "growth_rate.csv")
