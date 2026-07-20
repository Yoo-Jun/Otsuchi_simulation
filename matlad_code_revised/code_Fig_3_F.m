clc; close all; clear all;
load('dataset_modeldomain')
load('dataset_river')
load('dataset_initialcondition')
load('dataset_grid_ocean')

figure(1)
nx=300;
ny=250;
lon_high_1d=linspace(min(min(XLON_O)),max(max(XLON_O)),nx*4);
lat_high_1d=linspace(min(min(XLAT_O)),max(max(XLAT_O)),ny*4);
[lon_high,lat_high]=meshgrid(lon_high_1d,lat_high_1d);

% Topo_high=interp2(XLON_sub,XLAT_sub,TOPO2_sub,lon_high,lat_high);


Mtop2d_land=double(isnan(Mtop2d));
Mtop2d_high=interp2(XLON_O,XLAT_O,Mtop2d_land,lon_high,lat_high);
Mtop2d_high(Mtop2d_high==0)=nan;
Mtop2d_high(1:ny*4-1,1:nx*4)=Mtop2d_high(1:ny*4-1,1:nx*4)+Mtop2d_high(2:ny*4,1:nx*4);
Mtop2d_high(2:ny*4,1:nx*4)=Mtop2d_high(1:ny*4-1,1:nx*4)+Mtop2d_high(2:ny*4,1:nx*4);
Mtop2d_high(1:ny*4,1:nx*4-1)=Mtop2d_high(1:ny*4,1:nx*4-1)+Mtop2d_high(1:ny*4,2:nx*4);
Mtop2d_high(1:ny*4,2:nx*4)=Mtop2d_high(1:ny*4,1:nx*4-1)+Mtop2d_high(1:ny*4,2:nx*4);
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
contourf(XLON_O,XLAT_O,geo,[-220:20:0],'linecolor','none')
colormap(mycmap)
caxis([-1500 1500])
hold on
[c,h]=contour(XLON_O,XLAT_O,geo,[-180 -160 -140 -120 -80 -60 -40 -20],'k','linewidth',1);
clabel(c,h,'FontSize',25,'LabelSpacing',4000)
xticks([141.9 142 142.1 142.2])
xticklabels({'141.9','142','142.1','142.2'})
yticks([39.3 39.35 39.4 39.45])
yticklabels({'39.3','39.35','39.4','39.45'})

[c,h]=contour(XLON_O,XLAT_O,geo,[-200 -100],'k','linewidth',2);
clabel(c,h,'FontSize',25,'LabelSpacing',500)
mesh(lon_high,lat_high,Mtop2d_high,'edgecolor','none','facecolor','[0.5 0.5 0.5]')
daspect([1 1 1])
grid on
ax = gca;
ax.FontSize = 40;
lon0 = 142.15;
lat0 = 39.3;
dlon_1km = 5 / (111.32 * cosd(lat0));
plot([lon0, lon0 + dlon_1km],[lat0, lat0],'m','linewidth',5);

figure(2)
tile = tiledlayout(1,1);
ax1 = axes(tile);
plot(ax1,sa11,-griddep,'-.','color',[0.8500 0.3250 0.0980],'linewidth',2)
ax1.XColor = [0.8500 0.3250 0.0980];
ax1.YColor = 'k';
xlabel('Salinity (PSU)')
ylim([-200 0])
yticks([-200 -150 -100 -50 0])
% xlim([33.2 34.1])
pbaspect([1 2 1])
ax2 = axes(tile);
plot(ax2,ta11,-griddep,'-k','linewidth',2)
ax2.XAxisLocation = 'top';
ax2.YAxisLocation = 'right';
ylim([-200 0])
yticks([-200 -175 -150 -125 -100 -75 -50 -25 0])
yticklabels({'','','','',''})
xlim([9 15])
xticks([9 10 11 12 13 14 15])
xticklabels({'9','','11','','13','','15'})
xlabel('Temperature (^oC)')
ax2.Color = 'none';
ax1.Box = 'off';
ax2.Box = 'off';
grid on
pbaspect([1 2 1])

figure(3)
colororder({'k','k'})
plot(0:360,ORHL,'b','linewidth',3)
hold on
plot(0:360,URHL,'r','linewidth',1)
grid on
ylabel('Q (m^3/s)')
ylim([0 40])
legend('Or + Kr','Ur')
xticks(linspace(0,360,31))
xticklabels({'0623','','0624','','0625','','0626','','0627','','0628','','0629' ...
    ,'','0630','','0701','','0702','','0703','','0704','','0705','','0706','','0707','','0708'})
xtickangle(90)
xlim([0 360-12-4])
ax = gca;
ax.FontSize = 40; 
pbaspect([5 2 1])

