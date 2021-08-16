function [c] = Mohr_Circle_line( sigma1, sigma3)
%Plot Mohr circle in line.
%
Mohr_angle = linspace(0,pi,10000)';
Mohr_Center = (sigma1+sigma3)/2;
Mohr_Radius = (sigma1-sigma3)/2;
Mohr_x = Mohr_Center+Mohr_Radius.*cos(Mohr_angle);
Mohr_y = Mohr_Radius.*sin(Mohr_angle);
%
c = plot(Mohr_x,Mohr_y,'-','Color',[0.8 0.8 0.8],'LineWidth',2); 
end

