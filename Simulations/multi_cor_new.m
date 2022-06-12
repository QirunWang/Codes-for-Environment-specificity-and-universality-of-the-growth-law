function chi=multi_cor_new(sigma_knorm,mean_knorm,k_norm,sigma_anorm,mean_anorm,a_norm,rho_xk,rho_xa,N,chiR)


%Generate a vector 
%whose correlation coefficients with other two vectors are known.



mu = [0 mean_knorm mean_anorm];
sigma = [1                   rho_xk*sigma_knorm  rho_xa*sigma_anorm;
         rho_xk*sigma_knorm  sigma_knorm^2       0;
         rho_xa*sigma_anorm  0                   sigma_anorm^2];

if det(sigma)>0
R = chol(sigma)';

temp1=randn(N,1);
tempk=k_norm-mu(:,2);
tempa=a_norm-mu(:,3);
temp3=tempa/R(3,3);
temp2=(tempk-R(3,2)*temp3)/R(2,2);

z = repmat(mu,N,1) + [temp1,temp2,temp3]*R;


temp_i=1;
while temp_i ~=0
    CV_chi=rand(1)+4.5;
    sigmay=sqrt(log(1+CV_chi^2));
    meany=log(1)-sigmay^2/2;       
    chi_oth=exp(meany+sigmay*randn(N,1));
    chi_oth=(1-chiR)*chi_oth/sum(chi_oth);
    temp3= var(chi_oth)^0.5/mean(chi_oth);
    if temp3<5.5&&temp3>4.5
        temp_i=0;
    end
end

[~,i1] = sort(z(:,1));
chi_oth(i1) = sort(chi_oth);
chi=chi_oth;




else
chi=-1*ones(N,1);    
end   
    
end
