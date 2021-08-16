function [ Crk_inf_temp, Delta_Epsilon] = Unloading(Length_update,Width_update,Sigma_xx,Sigma_yy,Sigma_xy,New_Crk_inf,Crk_num,k_s,k_n,Pp)
%UNLOADING 此处显示有关此函数的摘要
%   此处显示详细说明
Far_StressTensor=[Sigma_xx Sigma_xy;
                  Sigma_xy Sigma_yy;];  
Crk_inf_temp = New_Crk_inf;
Delta_Epsilon = zeros(2,2);
for i=1:Crk_num
    alpha = New_Crk_inf(i,4);
    Unit_normal = [cos(alpha);sin(alpha);];
    Vector_Sg = Unit_normal'*Far_StressTensor*(eye(2)-kron(Unit_normal,Unit_normal'));
    Vector_S = Vector_Sg'/norm(Vector_Sg);
    Normal_stress =  Unit_normal'*Far_StressTensor*Unit_normal-Pp;
    Shear_stress = Unit_normal'*Far_StressTensor*Vector_S;
    % update stress/strain after unloading
    if New_Crk_inf(i,6) == 0 % Elastic fractures
        Crk_inf_temp(i,9) = Normal_stress;% Normal stress
        Crk_inf_temp(i,10) = Shear_stress;% Shear stress
        Crk_inf_temp(i,11) = Shear_stress*k_s/(New_Crk_inf(i,5)+k_s);% Shear stress on fracture
        Crk_inf_temp(i,12) = New_Crk_inf(i,12)-(New_Crk_inf(i,10)-Shear_stress)/(New_Crk_inf(i,5)+k_s);% Shear disp
    else % Plastic fractures
        Crk_inf_temp(i,9) = Normal_stress;% Normal stress
        Crk_inf_temp(i,10) = Shear_stress;% Shear stress
        if Shear_stress > New_Crk_inf(i,3)* Normal_stress%still plastic
            Crk_inf_temp(i,11) = New_Crk_inf(i,3)*Normal_stress;% Shear stress on fracture
            Crk_inf_temp(i,12) = New_Crk_inf(i,12)-(New_Crk_inf(i,11)-Crk_inf_temp(i,11))/k_s;% Shear disp
            Crk_inf_temp(i,6) = Shear_stress -  Crk_inf_temp(i,11); % stress drop
        else % elastic
            Crk_inf_temp(i,11) = Shear_stress*k_s/(New_Crk_inf(i,5)+k_s);% Shear stress on fracture
            Crk_inf_temp(i,12) = New_Crk_inf(i,12) - (New_Crk_inf(i,11)-Crk_inf_temp(i,11))/k_s;% Shear disp
            Crk_inf_temp(i,6) = 0;% stress drop
        end
    end
    %
    Crk_inf_temp(i,13) = New_Crk_inf(i,13)-(-Normal_stress/k_n);% normal disp
    %
    Normal_Vector = [cos(alpha);sin(alpha);];
    Disp_UnitVector = New_Crk_inf(i,7:8)';
    [D_Epsilon_Crk ] = Disp_Strain(Length_update,Width_update,Crk_inf_temp(i,2),Crk_inf_temp(i,12),Crk_inf_temp(i,13),Normal_Vector,Disp_UnitVector);
    %
    Delta_Epsilon = Delta_Epsilon + D_Epsilon_Crk;
end
%
end

