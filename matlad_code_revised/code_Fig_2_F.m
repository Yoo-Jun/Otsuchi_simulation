clc; clear all; close all; 

load('dataset_obs_TS')
load('dataset_TS')

% period=10;


a=[20140626;
    20140728;
    20141030;
    20150216;
    20150317;
    20150423;
    20150527;
    20150612;
    20150624;
    20150707;
    20150722;
    20150819;
    20150903;
    20151023;
    20151106];
[c,h]=contour(si,thetai,dens,[21 22 23 24 25 26 27 28],':','linecolor',[0.6 0.6 0.6]);
hold on
clabel(c,h,'LabelSpacing',1000,'color',[0.6 0.6 0.6],'fontsize',30);
xlabel('Salinity (psu)','FontSize',30)
ylabel('Temperature (^oC)','FontSize',30)
hold on;
axis([29.5,34.5,0,23])

p1=plot(obs_sal(:,1),obs_temp(:,1),'o','color',[0.5 0.5 0.5],'linewidth',2);
p2=plot(obs_sal(:,2),obs_temp(:,2),'x','color',[0.5 0.5 0.5],'linewidth',2);
p3=plot(obs_sal(:,3),obs_temp(:,3),'s','color',[0.5 0.5 0.5],'linewidth',2);
% p4=plot(obs_sal(:,4),obs_temp(:,4),'bs','linewidth',10);
p5=plot(obs_sal(:,5),obs_temp(:,5),'+','color',[0.5 0.5 0.5],'linewidth',2);
p6=plot(obs_sal(:,6),obs_temp(:,6),'*','color',[0.5 0.5 0.5],'linewidth',2);
p7=plot(obs_sal(:,7),obs_temp(:,7),'v','color',[0.5 0.5 0.5],'linewidth',2);
p8=plot(obs_sal(:,8),obs_temp(:,8),'d','color',[0.5 0.5 0.5],'linewidth',2);
p9=plot(obs_sal(:,9),obs_temp(:,9),'>','color',[0.5 0.5 0.5],'linewidth',2);
% p11=plot(obs_sal(:,11),obs_temp(:,11),'ks','linewidth',10);
p12=plot(obs_sal(:,12),obs_temp(:,12),'<','color',[0.5 0.5 0.5],'linewidth',2);
p13=plot(obs_sal(:,13),obs_temp(:,13),'^','color',[0.5 0.5 0.5],'linewidth',2);
p14=plot(obs_sal(:,14),obs_temp(:,14),'p','color',[0.5 0.5 0.5],'linewidth',2);
p15=plot(obs_sal(:,15),obs_temp(:,15),'h','color',[0.5 0.5 0.5],'linewidth',2);
p4=plot(obs_sal(:,4),obs_temp(:,4),'_','color',[0.5 0.5 0.5],'linewidth',2);
p11=plot(obs_sal(:,11),obs_temp(:,11),'|','color',[0.5 0.5 0.5],'linewidth',2);
p10=plot(obs_sal(:,10),obs_temp(:,10),'rs','linewidth',10);
% p11=plot(obs_sal(:,11),obs_temp(:,11),'ks','linewidth',10);
legend([p1 p2 p3 p4 p5 p6 p7 p8 p9 p10 p11 p12 p13 p14 p15],...
    {'20140626','20140728','20141030','20150216','20150317','20150423','20150527','20150612','20150624','20150707',...
    '20150722','20150819','20150903','20151023','20151106'},'Location','northeastoutside')


pbaspect([1 1 1])
ax = gca;
ax.FontSize = 30; 
grid on

