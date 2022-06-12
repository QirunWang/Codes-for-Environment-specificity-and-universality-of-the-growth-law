%Draw numbers from a known dataset 
%The drawn numbers will have a normal distribution


%temp1 is a n*3 matrix
%temp1(:,1) is the xs
%temp1(:,2) is the ys
%temp1(:,3) is the orders

%mean is [mu_x mu_y]
%sigma is [s_x^2  , cov_xy;
%          cov_xy , s_x^2 ]

%seq defines bin range in the sampling

%times = 1

function res=sampling_2norm(temp1,mean,sigma,sep,times)



temp2=mvncdf(mean-sep,mean+sep,mean,sigma);
temp3=temp2;


temp4=length(temp1(temp1(:,1)<(mean(1)+sep)& ...
                   temp1(:,2)<(mean(2)+sep)& ...
                   temp1(:,1)>(mean(1)-sep)& ...
                   temp1(:,2)>(mean(2)-sep), ...
                   1));
temp5=temp4/(temp3*times);
%calculate the total number



temp11=[];
for i=min(temp1(:,1)):sep:max(temp1(:,1))
    for j=min(temp1(:,2)):sep:max(temp1(:,2))
        
        temp6=mvncdf([i,j],[i+sep,j+sep],mean,sigma);
        temp7=temp6;
        temp8=floor(temp5*temp7);
         %find out how many numbers should fall in [i,j;i+sep,j+sep]
        
        if temp8==0
        continue
        else
        temp9=temp1(temp1(:,1)<(i+sep)& ...
                   temp1(:,2)<(j+sep)& ...
                   temp1(:,1)>i& ...
                   temp1(:,2)>j, ...
                   :);
        if length(temp9(:,1))<temp8
            temp10=temp9;
        else
            temp10=temp9(1:temp8,:);  
        end

        temp11=[temp11;temp10];
        end 
         
        res=temp11;
    end

end
%test the normal distribution
%hist3([temp11(:,1),temp11(:,2)],[50 50])
%histogram(temp11(temp11(:,2)<(mean(2)+sep)&temp11(:,2)>(mean(2)-sep),1))
%jbtest(temp11(temp11(:,2)<(mean(2)+sep)&temp11(:,2)>(mean(2)-sep),1))














end