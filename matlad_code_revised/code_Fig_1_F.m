clc; clear all; close all;
%need m_map %need dataset 
%e.g.,
% addpath('H:\Otsuchi\code_9\Figure_codes_revised\dataset_revised')
% addpath('H:\Otsuchi\m_map1.4\m_map')
%

load('dataset_coastline_river')
load('dataset_geometry')
load('dataset_geometry_2')
load('dataset_obs_station')
mycmap=[0    0.6000    0.8000
    0.1149    0.6907    0.8666
    0.2181    0.7698    0.9215
    0.3096    0.8372    0.9648
    0.3895    0.8930    0.9964
    0.4544    0.9354    1.0000
    0.5033    0.9640    1.0000
    0.5446    0.9830    1.0000
    0.5866    0.9965    1.0000
    0.6357    1.0000    1.0000
    0.6878    1.0000    1.0000
    0.7384    1.0000    1.0000
    0.7835    1.0000    1.0000
    0.8222    0.9930    0.9930
    0.8631    0.9783    0.9783
    0.8958    0.9570    0.9570
    0.9075    0.9249    0.9249
    0.8851    0.8768    0.8768
    0.8223    0.8075    0.8075
    0.7371    0.7271    0.7271
    0.6507    0.6483    0.6483
    0.5832    0.5832    0.5832
    0.5306    0.5306    0.5306
    0.4841    0.4841    0.4841
    0.4397    0.4397    0.4397
    0.3931    0.3931    0.3931
    0.3448    0.3448    0.3448
    0.2966    0.2966    0.2966
    0.2483    0.2483    0.2483
    0.2000    0.2000    0.2000];
figure(1)
m_proj('Mercator','long',[115 180],'lat',[15 55]);
m_gshhs_l('patch',[0.5 0.5 0.5])
m_grid('box','fancy','tickdir','in','tickstyle','dd','xtick',([120, 130, 140, 150, 160, 170, 180]),'ytick',([20, 30, 40, 50]),...
    'fontsize',40);
hold on
m_plot([136.818 147.234],[35.595 35.595],'color',[0.4940 0.1840 0.5560],'linewidth',3)
m_plot([136.818 147.234],[43.058 43.058],'color',[0.4940 0.1840 0.5560],'linewidth',3)
m_plot([136.818 136.818],[35.595 43.058],'color',[0.4940 0.1840 0.5560],'linewidth',3)
m_plot([147.234 147.234],[35.595 43.058],'color',[0.4940 0.1840 0.5560],'linewidth',3)

m_plot([139.639 144.355],[37.855 37.855],'color',[0.4940 0.1840 0.5560],'linewidth',3)
m_plot([139.639 144.355],[40.875 40.875],'color',[0.4940 0.1840 0.5560],'linewidth',3)
m_plot([139.639 139.639],[37.855 40.875],'color',[0.4940 0.1840 0.5560],'linewidth',3)
m_plot([144.355 144.355],[37.855 40.875],'color',[0.4940 0.1840 0.5560],'linewidth',3)
daspect([1 1 1])
% lon0 = 142.50;
% lat0 = 39.1;
% dlon_100km = 100 / (111.32 * cosd(lat0));
% m_plot([lon0, lon0 + dlon_100km],[lat0, lat0],'m','linewidth',5);

figure(2)
m_proj('Mercator','long',[141.045, 142.775],'lat',[38.957,39.781]);
m_gshhs_h('patch',[0.5 0.5 0.5],'facealpha',0.5,'linewidth',4);
hold on
% TOPO3(TOPO3>1001)=1001;
TOPO3(TOPO3<-1501)=-1501;
m_contourf(gx,gy,TOPO2,[-1600 -1500 -1400 -1300 -1200 -1100 -1000 -900 -800 -700 -600 -500 -400 -300 -200 -100 ... 
    100 200 300 400 500 600 700 800 900 1000 1100 1200 1300 1400 1500 1600],'linecolor','k','linewidth',1);
[c,h]=m_contour(gx,gy,TOPO2,[-1500 -1000 -500 500 1000],'linecolor','k','linewidth',4);
clabel(c,h,'FontSize',25,'LabelSpacing',5000)
colorbar('fontsize',20)
colormap(mycmap)
caxis([-1600 1600])
m_grid('box','fancy','tickdir','in','tickstyle','dd','xtick',[141.5 142 142.5],'ytick',[39 39.5],'fontsize',30);
hold on
m_plot([141.818 142.222],[39.283 39.283],'Color', [0.4940 0.1840 0.5560],'linewidth',3)
m_plot([141.818 142.222],[39.506 39.506],'Color', [0.4940 0.1840 0.5560],'linewidth',3)
m_plot([141.818 141.818],[39.283 39.506],'Color', [0.4940 0.1840 0.5560],'linewidth',3)
m_plot([142.222 142.222],[39.283 39.506],'Color', [0.4940 0.1840 0.5560],'linewidth',3)
lon0 = 142.50;
lat0 = 39.1;
dlon_10km = 10 / (111.32 * cosd(lat0));
m_plot([lon0, lon0 + dlon_10km],[lat0, lat0],'m','linewidth',5);


figure(3)
m_proj('Mercator','long',[141.88, 142.01],'lat',[39.305, 39.39]);
m_contourf(gxx,gyy,TOPO3,[-1600 -1500 -1400 -1300 -1200 -1100 -1000 -900 -800 -700 -600 -500 -400 -300 -200 -100 ... 
    100 200 300 400 500 600 700 800 900 1000 1100 1200 1300 1400 1500 1600],'linecolor','none','linewidth',1);
colormap(mycmap)
caxis([-1500 1500])
hold on
m_contour(gxx,gyy,TOPO3,[20:20:500],'k','linewidth',1);
[c,h]=m_contour(gxx,gyy,TOPO3,[100 200 300 400 500],'k','linewidth',2);
clabel(c,h,'FontSize',25,'LabelSpacing',1000)
[c,h]=m_contour(gxx,gyy,TOPO3,[-200 -100],'k','linewidth',2);
clabel(c,h,'FontSize',25,'LabelSpacing',10000)
[c,h]=m_contour(gxx,gyy,TOPO3,[-120 -100 -80 -60 -40 -20],'k');
clabel(c,h,'FontSize',25,'LabelSpacing',500)

lon0 = 141.99;
lat0 = 39.32;
dlon_1km = 1 / (111.32 * cosd(lat0));
m_plot([lon0, lon0 + dlon_1km],[lat0, lat0],'m','linewidth',5);

m_patch(coast_m(:,2),coast_m(:,1),[0.5 0.5 0.5],'facealpha',0.5,'EdgeColor','k','linewidth',2)
m_patch(coast_island(59:73,2),coast_island(59:73,1),[0.5 0.5 0.5],'facealpha',0.5,'EdgeColor','k','linewidth',2)
m_patch(coast_island(594:645,2),coast_island(594:645,1),[0.5 0.5 0.5],'facealpha',0.5,'EdgeColor','k','linewidth',2)
mc=mycmap;
m_plot(kotsuchi(1,:),kotsuchi(2,:),'color',[mc(16,:)],'linewidth',6)
m_plot(kotsuchi(1,:)+0.0005,kotsuchi(2,:)+0.0006,'k','linewidth',2.2)
m_plot(kotsuchi(1,:),kotsuchi(2,:)-0.0006,'k','linewidth',2.2)

m_plot(otsuchi(1,:),otsuchi(2,:),'color',[mc(16,:)],'linewidth',6)
m_plot(otsuchi(1,:)+0.0006,otsuchi(2,:)+0.0005,'k','linewidth',2.2)
m_plot(otsuchi(1,:)-0.0006,otsuchi(2,:)-0.0005,'k','linewidth',2.2)

m_plot(unosumai(1,:)-0.0007,unosumai(2,:)+0.00001,'k','linewidth',1.5)
m_plot(unosumai(1,:)+0.0006,unosumai(2,:)-0.00001,'k','linewidth',1.5)
m_plot(unosumai(1,:),unosumai(2,:),'color',[mc(16,:)],'linewidth',3)

m_grid('linestyle',...
    '--','box','fancy','tickdir','out','tickstyle','dd','xtick',[141.9 141.95 142],'ytick',[39.3 39.35 39.4],'fontsize',40);
 daspect([1 1 1])
p=m_plot(141+54.3/60,39+21.6/60,'bp','linewidth',10,'Markersize',10)
wo=m_plot(141.932,39.348,'r^','linewidth',10,'Markersize',10)
obs=m_plot(sta_14(:,1),sta_14(:,2),'ko','linewidth',5,'Markersize',5)
legend([p wo obs],{'Precipitation','Wind','CTD (42 sta.)'},'Location','northeastoutside','fontsize',30)
