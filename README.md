# Codes for *"*Environment-specificity and universality of the microbial growth law*"*
Qirun Wang & Jie Lin
20220611
#


## Softwares and packages
> MATLAB R2020b
* Curve Fitting Toolbox 3.5.12
* Parallel Computing Toolbox 7.3
> R studio 1.2.1335

> R 3.6.1
* readr 1.3.1
* readxl 1.3.1
* DESeq2 1.24.0
* clusterProfiler 3.12.0
* org.Sc.sgd.db 3.8.2
* ggplot2 3.3.0
* ggthemes 4.2.4
#

## Overview
> ### Simulations
  All MATLAB codes as well as example datasets are available in the directory *Simulations*. *Simu_kave.mlx*, *Simu_alpha0.mlx*, *Simu_all_noise_chi.mlx*, and *Simu_all_noise_k.mlx* are codes for figure 1, 2, 3 and S1, while *Simu_theory_validation.mlx* is used to validate the model with experimental data in figure 4a, b, and S2.
> ### Fit model to data
  All MATLAB codes, $\phi_R$ - $\mu$ data, and the fitting results are available in the directory *Data_fit*. *fit_Dai.mlx* and *fit_Metzl* are codes for figure 5.
> ### Data analysis
  All R codes, experimental data, and results are available in the directory *Data_analysis*. *sum_data_Friedrich2021.R* is the code for calculating $I_{\chi,k}$, $I_{\phi,\alpha}$, effective $m_R$, and $\phi_R$ from experimental data. It also exports the correlation coefficients between predicted $\phi_i$ and other parameters for table S2. *DDS_gsea.R* is the code for differential expression analysis and GSEA analysis in figure 4c and S3. The *Data* directory contains processed data, including protein degradation rates, elongation speeds, growth rates $\mu$, $\phi_i$ from proteomics, $\chi_i$ from ribo-seq, and ribosomal protein list.
#

## The structure of the directory
+ Read_me.md
+ Simulations
     + multi_cor_new.m
     + muphi_large.mat
     + muphi_noise_chi.mat
     + muphi_noise_k.mat
     + muphi_stm.mat
     + sampling_2norm.m
     + Simu_all.mlx
     + Simu_all_noise_chi.mlx
     + Simu_all_noise_k.mlx
     + Simu_alpha0.mlx
     + Simu_kave.mlx
     + Simu_theory_validation.mlx
+ Data_fit
     + Dai2016.xls
     + Dai_para.csv
     + fit_Dai.mlx
     + fit_Metzl.mlx
     + Metzl2017.xls
     + Metzl_para.csv     
+ Data_analysis
   + .RData
   + .Rhistory
   + cal_eff_mR.R
   + cor_phipre_chi.csv
   + Data_analysis.Rproj
   + DDS_gsea.R
   + sum_data_Friedrich2021.R
   + sum_data_Friedrich2021.csv
   + Data
     + degradation_rate
         +  a_1.csv  
     + elongation_speed
         +  Data_Riba.xlsx
         +  ki.csv   
     + growth_rate
         +  .RData
         +   .Rhistory
         +   growthcurve_new.R
         +  growth_rate.csv
         +  growth_rate.Rproj
         +  SCD.xlsx
         +  SCG.xlsx
     + proteomics
         +  data_calib_1
             + da_glu.csv
             + da_gly.csv
             + wt_glu.csv
             +  wt_gly.csv    
         +  data_L_calib_057
             + da_glu.csv
             + da_gly.csv
             + wt_glu.csv
             + wt_gly.csv  
         +  data_no_calib
             + da_glu.csv
             + da_gly.csv
             + wt_glu.csv
             + wt_gly.csv        
      + riboseq
         +  da-glu.csv
         +  da-gly.csv
         +  riboseq_res_minus.xlsx
         +  wt-glu.csv
         +  wt-gly.csv
      + ribosomal_protein
         +  ribopro_eu.csv
         +  ribo_eu_sub_mass.csv
#

## Descriptions of all files
+  ### Read_me.md
   
   The readme file. 
+  ### Simulations
   Contains all files for simualtions.
     + ### Simu_kave.mlx     
       The numeric simulations for figure 1a where all the elongation speed $k_i= \langle k \rangle$.
     + ### Simu_alpha0.mlx
       The numeric simulations for figure 1b where all the protein degradation rate $\alpha_i= 0$.
     + ### Simu_all.mlx
       The numeric simulations for figure 2 where $k_i$ and $\alpha_i$ are heterogeneous.
     + ### muphi_large.mat
       Example dataset for *Simu_all.mlx*.  
     + ### Simu_all_noise_k.mlx
       The numeric simulations for figure S1a with noised $k_i$.     
     + ### muphi_noise_k.mat
       Example dataset for *Simu_all_noise_k.mlx*.  
     + ### Simu_all_noise_chi.mlx
       The numeric simulations for figure S1b with noised $\chi_i$.  
     + ### muphi_noise_chi.mat
       Example dataset for *Simu_all_noise_chi.mlx*.  
     + ### Simu_theory_validation.mlx
       Validation of our theory with experimental data (figure 4a).
     + ### muphi_stm.mat
       Dataset of STM result.
     + ### multi_cor_new.m 
       The function is used to generate a set of numbers whose correlation coefficients with other two vectors are known.
     + ### sampling_2norm.m
       Sampling data points whose distribution will be two-dimensional Gaussian distribution.
+  ### Data_fit
   Contains codes for data fitting.
     + ### fit_Dai.mlx
       Fit our model to the data from Dai, et al., 2016 (Figure 5b and c).
     + ### Dai2016.xls
       $\phi_R$ - $\mu$ data from Dai, et al., 2016.
     + ### Dai_para.csv
       Detailed fitting results of Dai, et al., 2016 (Figure 5c).
     + ### fit_Metzl.mlx
       Fit our model to the data from Metzl, et al., 2017 (Figure 5a and c).
     + ### Metzl2017.xls
       $\phi_R$ - $\mu$ data from Metzl, et al., 2017.
     + ### Metzl_para.csv 
       Detailed fitting results of Metzl, et al., 2017 (Figure 5c).
+  ### Data_analysis
     + ### Data_analysis.Rproj & .RData & .Rhistory
       R studio project files.
     + ### sum_data_Friedrich2021.R
       Codes for processing the data from Friedrich, et al., 2021 (figure 4a, b and d).
     + ### sum_data_Friedrich2021.csv
       Results including $I_{\chi,k}$, $I_{\phi,\alpha}$, $\phi_R$, $\mu$, and effective $m_R$ of sum_data_Friedrich2021.R used in *Simu_theory_validation.mlx* and figure 4d.
     + ### cor_phipre_chi.csv
       Results of correlation coefficients between predicted $\phi_i$ and other parameters of *sum_data_Friedrich2021.R* (table S2).
     + ### cal_eff_mR.R
       Codes for calculating the effective molecular mass of a ribosome.
     + ### DDS_gsea.R
       Codes for differential expression analysis and GSEA functional analysis of genes (figure 4c and figure S3).


     + ### Data
         + ### degradation_rate
           Contains the degradation rates of proteins from Lahtvee et al., 2017.  
         + ### elongation_speed
            + ### Data_Riba.xlsx
               Calculation of elongation speeds from Riba et al., 2019. 
            + ### ki.csv 
               Elongation speeds from *Data_Riba*.xlsx.               
         + ### growth_rate
            + ### growth_rate.Rproj & .RData & .Rhistory
               R studio project files.
            + ### growthcurve_new.R
               The codes for calculating growth rates from OD600-time curves.            
            + ### growth_rate.csv
               Calculated growth rates from *growthcurve_new.R*.             
            + ### SCD.xlsx
              OD600-time curves of cell in SC+2% glucose medium (SCD) from Friedrich et al., 2021.
            + ### SCG.xlsx
              OD600-time curves of cell in SC+2% glycerol medium (SCG) from Friedrich et al., 2021.
         + ### proteomics
            + ### data_calib_1
              Contains xTop-corrected proteomics data from Friedrich, et al., 2021, in which the experimental $\phi_i$ is calibrated with $L^{−1}$. 
            + ### data_L_calib_057
              Contains xTop-corrected proteomics data from Friedrich, et al., 2021, in which the experimental $\phi_i$ is calibrated with $L^{−0.57}$. 
            + ### no_calib
              Contains xTop-corrected proteomics data from Friedrich, et al., 2021, in which the experimental $\phi_i$ is not calibrated.        
         + ### riboseq
            + ### riboseq_res_minus.xlsx
               Ribo-seq data read counts from Friedrich, et al., 2021. 
            + ### da-glu.csv & da-gly.csv & wt-glu.csv & wt-gly.csv
               Ribosome allocation fraction $\chi$ calculated from ribo-seq data.                
         + ### ribosomal_protein
            + ### ribopro_eu.csv
               All ribosomal protein genes except ribosomes in mitochondria from SGD.            
            + ### ribo_eu_sub_mass.csv
               The molecular mass of all ribosomal proteins except ribosomes in mitochondria.  

