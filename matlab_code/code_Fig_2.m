clc; clear all; close all; 

load('dataset_obs_TS')


smin=29;
smax=35;
thetamin=0;
thetamax=23;
xdim=round((smax-smin)*10+1);
ydim=round((thetamax-thetamin)*10+1);
dens=zeros(ydim,xdim);
thetai=((1:ydim)-1)*1+thetamin;
si=((1:xdim)-1)*0.1+smin;

for j=1:ydim
    for i=1:xdim
        dens(j,i)=sw_dens0(si(i),thetai(j));
    end
end
dens=dens-1000;

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

p1=plot(obs_sal(:,1),obs_temp(:,1),'ko','linewidth',2);
p2=plot(obs_sal(:,2),obs_temp(:,2),'bo','linewidth',2);
p3=plot(obs_sal(:,3),obs_temp(:,3),'o','color',[0.6350 0.0780 0.1840],'linewidth',2);
p4=plot(obs_sal(:,4),obs_temp(:,4),'mo','linewidth',2);
p5=plot(obs_sal(:,5),obs_temp(:,5),'o','color',[0.9290 0.6940 0.1250],'linewidth',2);
p6=plot(obs_sal(:,6),obs_temp(:,6),'o','color',[0.3010 0.7450 0.3],'linewidth',2);
p7=plot(obs_sal(:,7),obs_temp(:,7),'go','linewidth',2);
p8=plot(obs_sal(:,8),obs_temp(:,8),'o','color',[1 0.4470 0.7410],'linewidth',2);
p9=plot(obs_sal(:,9),obs_temp(:,9),'o','color',[0.8500 0.3250 0.0980],'linewidth',2);
p11=plot(obs_sal(:,11),obs_temp(:,11),'o','color',[1 0.30 0.9330],'linewidth',2);
p12=plot(obs_sal(:,12),obs_temp(:,12),'co','linewidth',2);
p13=plot(obs_sal(:,13),obs_temp(:,13),'o','color',[0.4660 0.6740 0.1880],'linewidth',2);
p14=plot(obs_sal(:,14),obs_temp(:,14),'o','color',[0.4940 0.5 0.5560],'linewidth',2);
p15=plot(obs_sal(:,15),obs_temp(:,15),'o','color',[0.4940 0.1840 0.5560],'linewidth',2);
p10=plot(obs_sal(:,10),obs_temp(:,10),'rs','linewidth',10);
legend([p1 p2 p3 p4 p5 p6 p7 p8 p9 p10 p11 p12 p13 p14 p15],...
    {'20140626','20140728','20141030','20150216','20150317','20150423','20150527','20150612','20150624','20150707',...
    '20150722','20150819','20150903','20151023','20151106'},'Location','northeastoutside')

pbaspect([1 1 1])
ax = gca;
ax.FontSize = 30; 
grid on

