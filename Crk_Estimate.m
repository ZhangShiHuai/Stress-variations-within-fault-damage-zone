function [Crk_inf,Critical_inf,NonCri_inf] = Crk_Estimate(Crk_num,Crk_len,Crk_mu,Crk_Orien,Shear_G,PoissonRatio,Sigma_xx,Sigma_yy,Sigma_xy,k_s,k_n,Pp,Dila_coe)
% Input stress conditions and fracture constitutive parameters
% Output fracture information in a matrix "Crk_inf"
Crk_inf = zeros(Crk_num,14); % Crk_No,Crk_len, Crk_mu, Crk_Orien,Crack km (local matrix-fracture stiffness),Stress drop,Relative displacement vector, Normal_stress, Shear stress, fracture shear stress, mean displacement, tao_p <---> All cracks
Crk_inf(:,1) = linspace(1,Crk_num,Crk_num)';
Crk_inf(:,2) = Crk_len;
Crk_inf(:,3) = Crk_mu;
Crk_inf(:,4) = Crk_Orien;
Crk_inf(:,5) = 4*Shear_G./(Crk_len*pi*(1-PoissonRatio)); % km

for Crk_No = 1:Crk_num
    [Stress_drop,Vector_S,Normal_stress,Shear_stress,ShearStress_crk, Mean_disp_Shear,Mean_disp_Normal, tao_p] = Crack_Criticality(Sigma_xx,Sigma_yy,Sigma_xy,Crk_Orien(Crk_No),Crk_mu(Crk_No),Crk_inf(Crk_No,5), k_s, k_n, Pp,Dila_coe);

    Crk_inf(Crk_No,6) = Stress_drop;
    Crk_inf(Crk_No,7:8) = -Vector_S';
    Crk_inf(Crk_No,9) = Normal_stress;
    Crk_inf(Crk_No,10) = Shear_stress;
    Crk_inf(Crk_No,11) = ShearStress_crk;
    Crk_inf(Crk_No,12) = Mean_disp_Shear;
    Crk_inf(Crk_No,13) = Mean_disp_Normal;
    Crk_inf(Crk_No,14) = tao_p;
end

Critical_inf = Crk_inf(find(Crk_inf(:,6)),:);
NonCri_inf = Crk_inf(Crk_inf(:,6)==0,:);

end

