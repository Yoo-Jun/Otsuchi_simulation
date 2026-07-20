clc; clear all; close all;

load('dataset_coastline_river.mat')
load('dataset_wind_D4')
load('dataset_wind_D4_m')
load('dataset_grid_atmosphere')

LL=0.03;
%
hgt_r=hgt;
mask_r=mask;
mycmap=[0.3686    0.3098    0.6353
    0.2395    0.4158    0.7107
    0.1958    0.5226    0.7407
    0.2536    0.6369    0.7101
    0.3764    0.7438    0.6533
    0.5000    0.8105    0.6409
    0.6334    0.8532    0.6452
    0.7556    0.9009    0.6232
    0.8668    0.9465    0.5932
    0.9427    0.9807    0.6356
    0.9866    1.0000    0.7329
    1.0000    0.9838    0.7243
    1.0000    0.9233    0.6109
    0.9958    0.8466    0.5106
    0.9958    0.7548    0.4304
    0.9895    0.6484    0.3597
    0.9780    0.5220    0.2934
    0.9490    0.4074    0.2612
    0.8985    0.3223    0.2879
    0.8272    0.2330    0.3096
    0.7341    0.1247    0.2953
    0.6196    0.0039    0.2588];


%figure 5a
figure(1)
lat_raw=lat;
lat=lat_raw;
lon_raw=lon;
lon=lon_raw;



u10_raw=u10_hourly_sim(:,:,1);
v10_raw=v10_hourly_sim(:,:,1);
u10=u10_raw+maskct;
v10=v10_raw+maskct;

u10(isnan(u10))=0;
v10(isnan(v10))=0;

u10x=zeros(size(u10));
u10y=zeros(size(u10));
v10x=zeros(size(u10));
v10y=zeros(size(u10));
u10x(2:end,:,:)=1.225*(1.5*10^-3)*(abs(u10(2:end,:,:)).*u10(2:end,:,:)-abs(u10(1:end-1,:,:)).*u10(1:end-1,:,:))/100;
u10y(:,2:end,:)=1.225*(1.5*10^-3)*(abs(u10(:,2:end,:)).*u10(:,2:end,:)-abs(u10(:,1:end-1,:)).*u10(:,1:end-1,:))/100;
v10x(2:end,:,:)=1.225*(1.5*10^-3)*(abs(v10(2:end,:,:)).*v10(2:end,:,:)-abs(v10(1:end-1,:,:)).*v10(1:end-1,:,:))/100;
v10y(:,2:end,:)=1.225*(1.5*10^-3)*(abs(v10(:,2:end,:)).*v10(:,2:end,:)-abs(v10(:,1:end-1,:)).*v10(:,1:end-1,:))/100;

ws=u10x+v10y;
[xend,yend]=size(squeeze(u10));

lon_map(1:6:xend,1:6:yend)=lon(1:6:xend,1:6:yend);
lon_map(4:6:xend,4:6:yend)=lon(4:6:xend,4:6:yend);
lat_map(1:6:xend,1:6:yend)=lat(1:6:xend,1:6:yend);
lat_map(4:6:xend,4:6:yend)=lat(4:6:xend,4:6:yend);

u10_map(1:6:xend,1:6:yend)=u10(1:6:xend,1:6:yend);
u10_map(4:6:xend,4:6:yend)=u10(4:6:xend,4:6:yend);
v10_map(1:6:xend,1:6:yend)=v10(1:6:xend,1:6:yend);
v10_map(4:6:xend,4:6:yend)=v10(4:6:xend,4:6:yend);
u10_map_r(1:6:xend,1:6:yend)=u10_raw(1:6:xend,1:6:yend)+mask(1:6:xend,1:6:yend);
u10_map_r(4:6:xend,4:6:yend)=u10_raw(4:6:xend,4:6:yend)+mask(4:6:xend,4:6:yend);
v10_map_r(1:6:xend,1:6:yend)=v10_raw(1:6:xend,1:6:yend)+mask(1:6:xend,1:6:yend);
v10_map_r(4:6:xend,4:6:yend)=v10_raw(4:6:xend,4:6:yend)+mask(4:6:xend,4:6:yend);

ws(1:75,end-170:end)=nan;
hgt(1:71,end-165:end)=nan;
mask(1:71,end-165:end)=nan;
% 
ws(ws<-1*10^-4)=-1*10^-4-10^-5;
ws(ws>1*10^-4)=1*10^-4+10^-5;
contourf(lon,lat,ws,linspace(-1*10^-4-10^-5,1*10^-4+10^-5,23),'linecolor','none')

colormap(mycmap)
colorbar
caxis([-1*10^-4-10^-5,1*10^-4+10^-5])
hold on

mesh(lon,lat,mask,'edgecolor','none','facecolor','[0.8 0.8 0.8]','FaceAlpha','1')
contour(lon,lat,hgt,[20:20:500],'linecolor',[0.5 0.5 0.5],'linewidth',1);
[C,h]=contour(lon,lat,hgt,[100 200 300 400 500],'linecolor',[0.5 0.5 0.5],'linewidth',3);
clabel(C,h,'FontSize',20,'LabelSpacing',10000)
u10_map(1:71,end-165:end)=nan;
v10_map(1:71,end-165:end)=nan;
u10_map(10,end-15)=20;
v10_map(10,end-15)=0;
u10_map(62,end-156)=10;
v10_map(62,end-156)=0;
u10_map_r(1:71,end-165:end)=nan;
v10_map_r(1:71,end-165:end)=nan;
u10_map_r(10,end-15)=20;
v10_map_r(10,end-15)=0;
lon_map(62,end-156)=lon(62,end-156);
lat_map(62,end-156)=lat(62,end-156);
plot(141.9326,39.3488,'r^','markersize',15,'linewidth',10)
plot(141.907,39.360,'bp','linewidth',10,'Markersize',15)
quiver(lon_map,lat_map,u10_map_r,v10_map_r,LL,'color',[0.3 0.3 0.3],'linewidth',2)
quiver(lon_map,lat_map,u10_map,v10_map,LL,'k','linewidth',4)
plot(coast_m(:,2),coast_m(:,1),'m:','linewidth',5)
text(141.885,39.362,'10 m/s','fontsize',50)
lon0 = 141.97;
lat0 = 39.362;
dlon_1km = 1 / (111.32 * cosd(lat0));
plot([lon0, lon0 + dlon_1km],[lat0, lat0],'m','linewidth',5);
grid on
ax = gca;
ax.FontSize = 40;
axis([141.885, 141.985,39.322, 39.368])
daspect([1.3 1 1])
xticks([141.9 141.95 142])
xticklabels({'141.9','141.95','142'})
yticks([39.35])
yticklabels({'39.35'})

%figure 5b
figure(2)
lat_raw=lat;
lat=lat_raw;
lon_raw=lon;
lon=lon_raw;

u10_raw=u10_hourly_sim(:,:,2);
v10_raw=v10_hourly_sim(:,:,2);
u10=u10_raw+maskct;
v10=v10_raw+maskct;

u10(isnan(u10))=0;
v10(isnan(v10))=0;

u10x=zeros(size(u10));
u10y=zeros(size(u10));
v10x=zeros(size(u10));
v10y=zeros(size(u10));
u10x(2:end,:,:)=1.225*(1.5*10^-3)*(abs(u10(2:end,:,:)).*u10(2:end,:,:)-abs(u10(1:end-1,:,:)).*u10(1:end-1,:,:))/100;
u10y(:,2:end,:)=1.225*(1.5*10^-3)*(abs(u10(:,2:end,:)).*u10(:,2:end,:)-abs(u10(:,1:end-1,:)).*u10(:,1:end-1,:))/100;
v10x(2:end,:,:)=1.225*(1.5*10^-3)*(abs(v10(2:end,:,:)).*v10(2:end,:,:)-abs(v10(1:end-1,:,:)).*v10(1:end-1,:,:))/100;
v10y(:,2:end,:)=1.225*(1.5*10^-3)*(abs(v10(:,2:end,:)).*v10(:,2:end,:)-abs(v10(:,1:end-1,:)).*v10(:,1:end-1,:))/100;

ws=u10x+v10y;
[xend,yend]=size(squeeze(u10));
hgt=hgt_r;
mask=mask_r;
lon_map=zeros(size(xend,yend));
lat_map=zeros(size(xend,yend));

lon_map(1:6:xend,1:6:yend)=lon(1:6:xend,1:6:yend);
lon_map(4:6:xend,4:6:yend)=lon(4:6:xend,4:6:yend);
lat_map(1:6:xend,1:6:yend)=lat(1:6:xend,1:6:yend);
lat_map(4:6:xend,4:6:yend)=lat(4:6:xend,4:6:yend);

u10_map=zeros(size(xend,yend));
v10_map=zeros(size(xend,yend));
u10_map_r=zeros(size(xend,yend));
v10_map_r=zeros(size(xend,yend));

u10_map(1:6:xend,1:6:yend)=u10(1:6:xend,1:6:yend);
u10_map(4:6:xend,4:6:yend)=u10(4:6:xend,4:6:yend);
v10_map(1:6:xend,1:6:yend)=v10(1:6:xend,1:6:yend);
v10_map(4:6:xend,4:6:yend)=v10(4:6:xend,4:6:yend);
u10_map_r(1:6:xend,1:6:yend)=u10_raw(1:6:xend,1:6:yend)+mask(1:6:xend,1:6:yend);
u10_map_r(4:6:xend,4:6:yend)=u10_raw(4:6:xend,4:6:yend)+mask(4:6:xend,4:6:yend);
v10_map_r(1:6:xend,1:6:yend)=v10_raw(1:6:xend,1:6:yend)+mask(1:6:xend,1:6:yend);
v10_map_r(4:6:xend,4:6:yend)=v10_raw(4:6:xend,4:6:yend)+mask(4:6:xend,4:6:yend);

% ws(1:75,end-170:end)=nan;
% hgt(1:71,end-165:end)=nan;
% mask(1:71,end-165:end)=nan;
% 
ws(ws<-1*10^-4)=-1*10^-4-10^-5;
ws(ws>1*10^-4)=1*10^-4+10^-5;
contourf(lon,lat,ws,linspace(-1*10^-4-10^-5,1*10^-4+10^-5,23),'linecolor','none')

colormap(mycmap)
colorbar
caxis([-1*10^-4-10^-5,1*10^-4+10^-5])
hold on

mesh(lon,lat,mask,'edgecolor','none','facecolor','[0.8 0.8 0.8]','FaceAlpha','1')
contour(lon,lat,hgt,[20:20:500],'linecolor',[0.5 0.5 0.5],'linewidth',1);
[C,h]=contour(lon,lat,hgt,[100 200 300 400 500],'linecolor',[0.5 0.5 0.5],'linewidth',3);
clabel(C,h,'FontSize',20,'LabelSpacing',10000)
% u10_map(1:71,end-165:end)=nan;
% v10_map(1:71,end-165:end)=nan;
u10_map(10,end-15)=20;
v10_map(10,end-15)=0;
% u10_map(62,end-156)=10;
% v10_map(62,end-156)=0;
% u10_map_r(1:71,end-165:end)=nan;
% v10_map_r(1:71,end-165:end)=nan;
u10_map_r(10,end-15)=20;
v10_map_r(10,end-15)=0;
% lon_map(62,end-156)=lon(62,end-156);
% lat_map(62,end-156)=lat(62,end-156);
plot(141.9326,39.3488,'r^','markersize',15,'linewidth',10)
plot(141.907,39.360,'bp','linewidth',10,'Markersize',15)
quiver(lon_map,lat_map,u10_map_r,v10_map_r,LL,'color',[0.3 0.3 0.3],'linewidth',2)
quiver(lon_map,lat_map,u10_map,v10_map,LL,'k','linewidth',4)
plot(coast_m(:,2),coast_m(:,1),'m:','linewidth',5)
% text(141.885,39.362,'10 m/s','fontsize',50)
grid on
ax = gca;
ax.FontSize = 40;
axis([141.885, 141.985,39.322, 39.368])
daspect([1.3 1 1])
xticks([141.9 141.95 142])
xticklabels({'141.9','141.95','142'})
yticks([39.35])
yticklabels({'39.35'})


figure(3)
lat_raw=lat;
lat=lat_raw;
lon_raw=lon;
lon=lon_raw;

u10_raw=m_u10_hour;
v10_raw=m_v10_hour;
u10=u10_raw+maskct;
v10=v10_raw+maskct;

u10(isnan(u10))=0;
v10(isnan(v10))=0;

u10x=zeros(size(u10));
u10y=zeros(size(u10));
v10x=zeros(size(u10));
v10y=zeros(size(u10));
u10x(2:end,:,:)=1.225*(1.5*10^-3)*(abs(u10(2:end,:,:)).*u10(2:end,:,:)-abs(u10(1:end-1,:,:)).*u10(1:end-1,:,:))/100;
u10y(:,2:end,:)=1.225*(1.5*10^-3)*(abs(u10(:,2:end,:)).*u10(:,2:end,:)-abs(u10(:,1:end-1,:)).*u10(:,1:end-1,:))/100;
v10x(2:end,:,:)=1.225*(1.5*10^-3)*(abs(v10(2:end,:,:)).*v10(2:end,:,:)-abs(v10(1:end-1,:,:)).*v10(1:end-1,:,:))/100;
v10y(:,2:end,:)=1.225*(1.5*10^-3)*(abs(v10(:,2:end,:)).*v10(:,2:end,:)-abs(v10(:,1:end-1,:)).*v10(:,1:end-1,:))/100;

ws=u10x+v10y;
[xend,yend]=size(squeeze(u10));
hgt=hgt_r;
mask=mask_r;
lon_map=zeros(size(xend,yend));
lat_map=zeros(size(xend,yend));

lon_map(1:6:xend,1:6:yend)=lon(1:6:xend,1:6:yend);
lon_map(4:6:xend,4:6:yend)=lon(4:6:xend,4:6:yend);
lat_map(1:6:xend,1:6:yend)=lat(1:6:xend,1:6:yend);
lat_map(4:6:xend,4:6:yend)=lat(4:6:xend,4:6:yend);

u10_map=zeros(size(xend,yend));
v10_map=zeros(size(xend,yend));
u10_map_r=zeros(size(xend,yend));
v10_map_r=zeros(size(xend,yend));

u10_map(1:6:xend,1:6:yend)=u10(1:6:xend,1:6:yend);
u10_map(4:6:xend,4:6:yend)=u10(4:6:xend,4:6:yend);
v10_map(1:6:xend,1:6:yend)=v10(1:6:xend,1:6:yend);
v10_map(4:6:xend,4:6:yend)=v10(4:6:xend,4:6:yend);
u10_map_r(1:6:xend,1:6:yend)=u10_raw(1:6:xend,1:6:yend)+mask(1:6:xend,1:6:yend);
u10_map_r(4:6:xend,4:6:yend)=u10_raw(4:6:xend,4:6:yend)+mask(4:6:xend,4:6:yend);
v10_map_r(1:6:xend,1:6:yend)=v10_raw(1:6:xend,1:6:yend)+mask(1:6:xend,1:6:yend);
v10_map_r(4:6:xend,4:6:yend)=v10_raw(4:6:xend,4:6:yend)+mask(4:6:xend,4:6:yend);

% ws(1:75,end-170:end)=nan;
% hgt(1:71,end-165:end)=nan;
% mask(1:71,end-165:end)=nan;
% 
ws(ws<-1*10^-4)=-1*10^-4-10^-5;
ws(ws>1*10^-4)=1*10^-4+10^-5;
contourf(lon,lat,ws,linspace(-1*10^-4-10^-5,1*10^-4+10^-5,23),'linecolor','none')

colormap(mycmap)
colorbar
caxis([-1*10^-4-10^-5,1*10^-4+10^-5])
hold on

mesh(lon,lat,mask,'edgecolor','none','facecolor','[0.8 0.8 0.8]','FaceAlpha','1')
contour(lon,lat,hgt,[20:20:500],'linecolor',[0.5 0.5 0.5],'linewidth',1);
[C,h]=contour(lon,lat,hgt,[100 200 300 400 500],'linecolor',[0.5 0.5 0.5],'linewidth',3);
clabel(C,h,'FontSize',20,'LabelSpacing',10000)
% u10_map(1:71,end-165:end)=nan;
% v10_map(1:71,end-165:end)=nan;
u10_map(10,end-15)=20;
v10_map(10,end-15)=0;
% u10_map(62,end-156)=10;
% v10_map(62,end-156)=0;
% u10_map_r(1:71,end-165:end)=nan;
% v10_map_r(1:71,end-165:end)=nan;
u10_map_r(10,end-15)=20;
v10_map_r(10,end-15)=0;
% lon_map(62,end-156)=lon(62,end-156);
% lat_map(62,end-156)=lat(62,end-156);
plot(141.9326,39.3488,'r^','markersize',15,'linewidth',10)
plot(141.907,39.360,'bp','linewidth',10,'Markersize',15)
quiver(lon_map,lat_map,u10_map_r,v10_map_r,LL,'color',[0.3 0.3 0.3],'linewidth',2)
quiver(lon_map,lat_map,u10_map,v10_map,LL,'k','linewidth',4)
plot(coast_m(:,2),coast_m(:,1),'m:','linewidth',5)
% text(141.885,39.362,'10 m/s','fontsize',50)
grid on
ax = gca;
ax.FontSize = 40;
axis([141.885, 141.985,39.322, 39.368])
daspect([1.3 1 1])
xticks([141.9 141.95 142])
xticklabels({'141.9','141.95','142'})
yticks([39.35])
yticklabels({'39.35'})