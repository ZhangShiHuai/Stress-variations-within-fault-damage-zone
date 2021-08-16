function [Stress_drop,Vector_S,Normal_stress,Shear_stress,ShearStress_crk, Mean_disp,Mean_disp_Normal, tao_p] = Crack_Criticality( Sigma_xx,Sigma_yy,Sigma_xy, alpha, mu, k_m, k_s,k_n,Pp,Dila_coe)
%STRESS_CRITICALITY 
% 
Far_StressTensor=[Sigma_xx Sigma_xy;
                  Sigma_xy Sigma_yy;];                       
%
Unit_normal = [cos(alpha);sin(alpha);];
Vector_Sg = Unit_normal'*Far_StressTensor*(eye(2)-kron(Unit_normal,Unit_normal'));
Vector_S = Vector_Sg'/norm(Vector_Sg);
%
% Normal_stress =  (Sigma_xx+Sigma_yy)/2+(Sigma_xx-Sigma_yy)*cos(2*alpha)/2-Sigma_xy*sin(2*alpha)-Pp;
% Shear_stress = abs((Sigma_xx-Sigma_yy)*sin(2*alpha)/2+Sigma_xy*cos(2*alpha));
%
Normal_stress =  Unit_normal'*Far_StressTensor*Unit_normal-Pp;
Shear_stress = Unit_normal'*Far_StressTensor*Vector_S;
tao_p = mu*Normal_stress;

if Shear_stress < tao_p*(k_m+k_s)/k_s
    Stress_drop = 0;
    Mean_disp = Shear_stress/(k_m+k_s);
    ShearStress_crk = k_s*Mean_disp;
    Mean_disp_Normal = -Normal_stress/k_n;
else
    Stress_drop = Shear_stress-tao_p;
    Mean_disp = (Shear_stress-tao_p)/k_m;
    ShearStress_crk = tao_p;
    Mean_disp_Normal =-Normal_stress/k_n + Dila_coe*(Mean_disp-tao_p/k_s);
end

end

