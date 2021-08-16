function [Del_Stress_yy,New_Stress_yy,Mat_Strain_yy,Delta_Epsilon,k,a,b,fp,Crk_inf_temp] = Bisection(a,Delta_Sigmayy_ini,Stress_xx_input,Stress_yy_input,Strain_yy_input,Stress_xy_input,PR_input,YM_input,Crk_num,k_s,k_n,Pp,Length_update,Width_update,New_Crk_inf,Strain_Crk_input,Theta,err)
format long

Del_Stress_yy = (a+Delta_Sigmayy_ini)/2;% minius sign
if Theta>=45
    if (Stress_yy_input-Del_Stress_yy) >= Stress_xx_input
        b = -(Stress_xx_input-Stress_yy_input);
    else
        b = Delta_Sigmayy_ini;
    end
else
    if (Stress_yy_input-Del_Stress_yy) <= Stress_xx_input
        b = -(Stress_xx_input-Stress_yy_input);
    else
        b = Delta_Sigmayy_ini;
    end
end
%
for k = 1:10000000
    %
    Del_Stress_yy = (a+b)/2;
    New_Stress_yy = Stress_yy_input-Del_Stress_yy; % Updated stress in y-direction
    %
    Del_Mat_Strain_yy = (1-PR_input^2)*Del_Stress_yy/YM_input;
    Mat_Strain_yy = Strain_yy_input-Del_Mat_Strain_yy;
    %
    % Unloading
    [Crk_inf_temp, Delta_Epsilon] = Unloading(Length_update,Width_update,Stress_xx_input,New_Stress_yy,Stress_xy_input,New_Crk_inf,Crk_num,k_s,k_n,Pp);
    %
    %
    New_Strain_yy = Mat_Strain_yy+Delta_Epsilon(2,2);
    
    fp = New_Strain_yy-Strain_yy_input;
    fa = Strain_Crk_input(2,2);
    %
    if (fp==0 ||abs((b - a))/2 < err)
        break
    end
    %
    if(fa * fp < 0)
        b=Del_Stress_yy;
    else
        a = Del_Stress_yy;
    end
    %
    if k>10000000
        disp('Non Convergence in Bisection')
    end
end

