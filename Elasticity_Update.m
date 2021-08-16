function [ New_PR, New_YM] = Elasticity_Update( Sigma_xx,Sigma_yy,Total_Strain )
%ELASTICITY_UPDATE 
%   
New_PR = (Total_Strain(1,1)*Sigma_yy-Total_Strain(2,2)*Sigma_xx)/(Sigma_xx+Sigma_yy)/(Total_Strain(1,1)-Total_Strain(2,2));
New_YM = (1+New_PR)*(1-2*New_PR)*(Sigma_xx+Sigma_yy)/(Total_Strain(1,1)+Total_Strain(2,2));
end

