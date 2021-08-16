function [ C_his,fc ] = Histogram_Statistics( statistic_objt,low_bnd,up_bnd,bin_num )
%HISTOGRAM_STATISTICS �˴���ʾ�йش˺�����ժҪ
%   �˴���ʾ��ϸ˵��
xbins = linspace(low_bnd,up_bnd,bin_num);
[num,fc]=hist(statistic_objt,xbins);
C_his = linspecer(bin_num);
ind=bin_num:-1:1;
C_his_use=C_his(ind,:);
%
figure
bar_width=fc(2)-fc(1);
for i=1:bin_num
    h=bar(fc(i),num(i),bar_width);
    set(h,'FaceColor',C_his_use(i,:),'EdgeColor','w');
    hold on
end
set(gcf,'Colormap',C_his);
colormap(flipud(C_his));
caxis([low_bnd up_bnd]);
cb = colorbar;
get(cb,'Position');
% cb.Label.String = 'frictional coefficient';
end

