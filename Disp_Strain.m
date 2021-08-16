function [D_Epsilon ] = Disp_Strain(Length,Width,Crk_len,Mean_disp_Shear,Mean_disp_Normal,Normal_Vector,Disp_UnitVector)
% Project fracture displacements onto the element boundary and calculate the model strain.
%   Pollard & Segall, 1987
%   Kachanov, 1992

Disp_Vector =Mean_disp_Shear*Disp_UnitVector+Mean_disp_Normal*Normal_Vector;

D_Epsilon_xx = -Crk_len/(2*Length*Width)*(2*Disp_Vector(1)*Normal_Vector(1)); % compression is positive
D_Epsilon_xy = -Crk_len/(2*Length*Width)*(Disp_Vector(1)*Normal_Vector(2)+Disp_Vector(2)*Normal_Vector(1)); % compression is positive
D_Epsilon_yy = -Crk_len/(2*Length*Width)*(2*Disp_Vector(2)*Normal_Vector(2));

D_Epsilon = [D_Epsilon_xx D_Epsilon_xy;
             D_Epsilon_xy D_Epsilon_yy;];
         
end

