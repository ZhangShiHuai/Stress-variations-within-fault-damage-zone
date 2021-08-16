function [ Strain_tensor ] = Strain_Solver(Sigma_xx, Sigma_yy, Sigma_xy, PoissonRatio, Young_Modulus)
%CONSTITUTIVERELATIONS 
%   Plane strain
Comp_mat = ((1-PoissonRatio^2)/Young_Modulus)*[1 -PoissonRatio/(1-PoissonRatio) 0;
                                               -PoissonRatio/(1-PoissonRatio) 1 0;
                                               0  0   1/(1-PoissonRatio);         ];
%
Strain_tensor_temp = Comp_mat*[Sigma_xx;Sigma_yy;Sigma_xy;];

Strain_tensor = [Strain_tensor_temp(1) Strain_tensor_temp(3);
                 Strain_tensor_temp(3) Strain_tensor_temp(2) ;];

end

