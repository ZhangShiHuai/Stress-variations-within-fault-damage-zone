function [c] = Mohr_Circle_dot( sigma1, sigma3)
%MOHR_CIRCLE �˴���ʾ�йش˺�����ժҪ
%   �˴���ʾ��ϸ˵��
Mohr_angle = linspace(0,pi,10000)';
Mohr_Center = (sigma1+sigma3)/2;
Mohr_Radius = (sigma1-sigma3)/2;
Mohr_x = Mohr_Center+Mohr_Radius.*cos(Mohr_angle);
Mohr_y = Mohr_Radius.*sin(Mohr_angle);
%
hold on
c = plot(Mohr_x,Mohr_y,'k.'); 
hold on
end

