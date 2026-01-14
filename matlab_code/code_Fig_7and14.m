clc; clear all; close all;

nx=300;
ny=250;
load('dataset_grid_ocean.mat')
load('dataset_SSS.mat')

%colormap
mcmap=[0.2422    0.1504    0.2
    0.2422    0.1504    0.6603
    0.2810    0.3228    0.9579
    0.1786    0.5289    0.9682
    0.2161    0.7843    0.5923
    0.6720    0.7793    0.2227
    0.9970    0.7659    0.2199
    0.9769    0.9839    0.0805];

%contour grid mismatch adjustment
lon_sub=XLON_O;
lat_sub=XLAT_O;
lon_high_1d=linspace(min(min(lon_sub)),max(max(lon_sub)),nx*4);
lat_high_1d=linspace(min(min(lat_sub)),max(max(lat_sub)),ny*4);
[lon_high,lat_high]=meshgrid(lon_high_1d,lat_high_1d);

Mtop2d_land=double(isnan(Mtop2d));

Mtop2d_high=interp2(lon_sub,lat_sub,Mtop2d_land,lon_high,lat_high);
Mtop2d_high(Mtop2d_high==0)=nan;
Mtop2d_high(1:ny*4-1,1:nx*4)=Mtop2d_high(1:ny*4-1,1:nx*4)+Mtop2d_high(2:ny*4,1:nx*4);
Mtop2d_high(2:ny*4,1:nx*4)=Mtop2d_high(1:ny*4-1,1:nx*4)+Mtop2d_high(2:ny*4,1:nx*4);
Mtop2d_high(1:ny*4,1:nx*4-1)=Mtop2d_high(1:ny*4,1:nx*4-1)+Mtop2d_high(1:ny*4,2:nx*4);
Mtop2d_high(1:ny*4,2:nx*4)=Mtop2d_high(1:ny*4,1:nx*4-1)+Mtop2d_high(1:ny*4,2:nx*4);


SSS_high_140830_U=interp2(lon_sub,lat_sub,SSS_140830_U',lon_high,lat_high);
SSS_high_140230_U=interp2(lon_sub,lat_sub,SSS_140230_U',lon_high,lat_high);
SSS_high_140830_H=interp2(lon_sub,lat_sub,SSS_140830_H',lon_high,lat_high);
SSS_high_140230_H=interp2(lon_sub,lat_sub,SSS_140230_H',lon_high,lat_high);
SSS_high_140830_ERA5=interp2(lon_sub,lat_sub,SSS_140830_ERA5',lon_high,lat_high);

%figures code

figure(1)
%Fig.7a
SSS_high_140830_H(SSS_high_140830_H<29.5)=29.5;
contourf(lon_high,lat_high,SSS_high_140830_H,linspace(29.5,33,8),'color',[0.4 0.4 0.4])
hold on
colormap(mcmap)
caxis([29.5,33.5])
axis([141.88, 142.01,39.32, 39.38])
mesh(lon_high,lat_high,Mtop2d_high,'edgecolor','none','facecolor','[0.5 0.5 0.5]')
colorbar('fontsize',50)
axis([141.885, 141.985,39.322, 39.368])
daspect([1.3 1 1])
xticks([141.9 141.95 142])
xticklabels({'141.9','141.95','142'})
yticks([39.35])
yticklabels({'39.35'})
grid on
ax = gca;
ax.FontSize = 40;

figure(2)
%Fig.7b
SSS_high_140830_U(SSS_high_140830_U<29.5)=29.5;
contourf(lon_high,lat_high,SSS_high_140830_U,linspace(29.5,33,8),'color',[0.4 0.4 0.4])
hold on
colormap(mcmap)
caxis([29.5,33.5])
axis([141.88, 142.01,39.32, 39.38])
mesh(lon_high,lat_high,Mtop2d_high,'edgecolor','none','facecolor','[0.5 0.5 0.5]')
colorbar('fontsize',50)
axis([141.885, 141.985,39.322, 39.368])
daspect([1.3 1 1])
xticks([141.9 141.95 142])
xticklabels({'141.9','141.95','142'})
yticks([39.35])
yticklabels({'39.35'})
grid on
ax = gca;
ax.FontSize = 40;

figure(3)
%Fig.7c
SSS_high_140230_H(SSS_high_140230_H<29.5)=29.5;
contourf(lon_high,lat_high,SSS_high_140230_H,linspace(29.5,33,8),'color',[0.4 0.4 0.4])
hold on
colormap(mcmap)
caxis([29.5,33.5])
axis([141.88, 142.01,39.32, 39.38])
mesh(lon_high,lat_high,Mtop2d_high,'edgecolor','none','facecolor','[0.5 0.5 0.5]')
colorbar('fontsize',50)
axis([141.885, 141.985,39.322, 39.368])
daspect([1.3 1 1])
xticks([141.9 141.95 142])
xticklabels({'141.9','141.95','142'})
yticks([39.35])
yticklabels({'39.35'})
grid on
ax = gca;
ax.FontSize = 40;

figure(4)
%Fig.7d
SSS_high_140230_U(SSS_high_140230_U<29.5)=29.5;
contourf(lon_high,lat_high,SSS_high_140230_U,linspace(29.5,33,8),'color',[0.4 0.4 0.4])
hold on
colormap(mcmap)
caxis([29.5,33.5])
axis([141.88, 142.01,39.32, 39.38])
mesh(lon_high,lat_high,Mtop2d_high,'edgecolor','none','facecolor','[0.5 0.5 0.5]')
colorbar('fontsize',50)
axis([141.885, 141.985,39.322, 39.368])
daspect([1.3 1 1])
xticks([141.9 141.95 142])
xticklabels({'141.9','141.95','142'})
yticks([39.35])
yticklabels({'39.35'})
grid on
ax = gca;
ax.FontSize = 40;

figure(5)
%Fig.12
SSS_high_140830_ERA5(SSS_high_140830_ERA5<29.5)=29.5;
contourf(lon_high,lat_high,SSS_high_140830_ERA5,linspace(29.5,33,8),'color',[0.4 0.4 0.4])
hold on
colormap(mcmap)
caxis([29.5,33.5])
axis([141.88, 142.01,39.32, 39.38])
mesh(lon_high,lat_high,Mtop2d_high,'edgecolor','none','facecolor','[0.5 0.5 0.5]')
colorbar('fontsize',50)
axis([141.885, 141.985,39.322, 39.368])
daspect([1.3 1 1])
xticks([141.9 141.95 142])
xticklabels({'141.9','141.95','142'})
yticks([39.35])
yticklabels({'39.35'})
grid on
ax = gca;
ax.FontSize = 40;

%