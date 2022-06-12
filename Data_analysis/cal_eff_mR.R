#######calculate effective ribosome mass


#input: 
#1 single-col dataframe cotaining ribosomal ORFs 
#2 ribosomal protein mass 
cal_eff_mR<-function(ribo_gene_df,ribo_mass)
{

  #ribo_gene_df<-as.vector(temp1$geneid)
  temp1<-data.frame(geneid=ribo_gene_df,des=rep("ribo",length(ribo_gene_df)))
  temp2<-merge(temp1,ribo_mass,by="geneid")
  temp3<-levels(as.factor(temp2$sub))
  mR<-0
  i=15
  for (i in 1:length(temp3)) 
  {
    temp4<-mean(temp2$mass[temp2$sub==temp3[i]])
    mR<-mR+temp4 
  }
  return(mR)
  
}

