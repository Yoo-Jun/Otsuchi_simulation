clc; clear all; close all;

load('dataset_coastline_river.mat')
load('dataset_grid_atmosphere')
load('dataset_grid_atmosphere_ERA')
load('dataset_grid_atmosphere_D1')
load('dataset_grid_atmosphere_D2')
load('dataset_grid_atmosphere_D3')
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

load('dataset_wind_ERA')
load('dataset_wind_D1')
load('dataset_wind_D2')
load('dataset_wind_D3')

LL=0.03;

% figure 6d
figure(1)

u10=u10_hourly_ERA;
v10=v10_hourly_ERA;

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
[xend,yend]=size(u10);


lon_map(1:9:xend,1:9:yend)=lon(1:9:xend,1:9:yend);
lon_map(5:9:xend,5:9:yend)=lon(5:9:xend,5:9:yend);
lat_map(1:9:xend,1:9:yend)=lat(1:9:xend,1:9:yend);
lat_map(5:9:xend,5:9:yend)=lat(5:9:xend,5:9:yend);

u10_map_r=zeros(size(xend,yend));
v10_map_r=zeros(size(xend,yend));
u10_map_r(1:9:xend,1:9:yend)=u10_hourly_ERA(1:9:xend,1:9:yend);
u10_map_r(5:9:xend,5:9:yend)=u10_hourly_ERA(5:9:xend,5:9:yend);
v10_map_r(1:9:xend,1:9:yend)=v10_hourly_ERA(1:9:xend,1:9:yend);
v10_map_r(5:9:xend,5:9:yend)=v10_hourly_ERA(5:9:xend,5:9:yend);

ws(1:75,end-170:end)=nan;

ws(ws<-1*10^-4)=-1*10^-4-10^-5;
ws(ws>1*10^-4)=1*10^-4+10^-5;

hold on

contourf(lon,lat,ws,linspace(-1*10^-4-10^-5,1*10^-4+10^-5,23),'linecolor','none');
colormap(mycmap)
colorbar
caxis([-1*10^-4-10^-5,1*10^-4+10^-5])
axis equal tight

mask1 = ~isnan(mask_ERA);  
mask1 = logical(mask1);

[nj, ni] = size(mask1);
lonc = nan(nj+1, ni+1);
latc = nan(nj+1, ni+1);

% interior corners
lonc(2:nj, 2:ni) = 0.25 * ( ...
    lon_era(1:end-1,1:end-1) + lon_era(2:end,1:end-1) + ...
    lon_era(1:end-1,2:end)   + lon_era(2:end,2:end) );

latc(2:nj, 2:ni) = 0.25 * ( ...
    lat_era(1:end-1,1:end-1) + lat_era(2:end,1:end-1) + ...
    lat_era(1:end-1,2:end)   + lat_era(2:end,2:end) );


plot(lon_era+mask_ERA-mask_ERA,lat_era+mask_ERA-mask_ERA,'k')

mesh(lon,lat,mask,'edgecolor','none','facecolor','[0.8 0.8 0.8]','FaceAlpha','1')
u10_map_r(10,end-15)=20;
v10_map_r(10,end-15)=0;
plot(141.9326,39.3488,'r^','markersize',15,'linewidth',10)
plot(141.907,39.360,'bp','linewidth',10,'Markersize',15)
xp=[141.8 141.8 141.901 141.901];
yp=[39.5041 39.359 39.359 39.5041];
quiver(lon_map,lat_map,u10_map_r,v10_map_r,LL,'color','k','linewidth',4)
grid on

ax = gca;
ax.FontSize = 40;
axis([141.885, 141.985,39.322, 39.368])
daspect([1.3 1 1])
xticks([141.9 141.95 142])
xticklabels({'141.9','141.95','142'})
yticks([39.35])
yticklabels({'39.35'})

% % figure 6c
figure(3)

u10=u10_hourly_D1;
v10=v10_hourly_D1;

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
[xend,yend]=size(u10);

lon_map(1:9:xend,1:9:yend)=lon(1:9:xend,1:9:yend);
lon_map(5:9:xend,5:9:yend)=lon(5:9:xend,5:9:yend);
lat_map(1:9:xend,1:9:yend)=lat(1:9:xend,1:9:yend);
lat_map(5:9:xend,5:9:yend)=lat(5:9:xend,5:9:yend);
u10_map_r=zeros(size(xend,yend));
v10_map_r=zeros(size(xend,yend));
u10_map_r(1:9:xend,1:9:yend)=u10_hourly_D1(1:9:xend,1:9:yend);
u10_map_r(5:9:xend,5:9:yend)=u10_hourly_D1(5:9:xend,5:9:yend);
v10_map_r(1:9:xend,1:9:yend)=v10_hourly_D1(1:9:xend,1:9:yend);
v10_map_r(5:9:xend,5:9:yend)=v10_hourly_D1(5:9:xend,5:9:yend);

ws(1:75,end-170:end)=nan;

ws(ws<-1*10^-4)=-1*10^-4-10^-5;
ws(ws>1*10^-4)=1*10^-4+10^-5;

hold on

contourf(lon,lat,ws,linspace(-1*10^-4-10^-5,1*10^-4+10^-5,23),'linecolor','none');
colormap(mycmap)
colorbar
caxis([-1*10^-4-10^-5,1*10^-4+10^-5])
axis equal tight

mask1_d1 = ~isnan(mask_D1);  
mask1_d1 = logical(mask1_d1);

[nj, ni] = size(mask1_d1);
lonc = nan(nj+1, ni+1);
latc = nan(nj+1, ni+1);

% interior corners
lonc(2:nj, 2:ni) = 0.25 * ( ...
    lon_d1(1:end-1,1:end-1) + lon_d1(2:end,1:end-1) + ...
    lon_d1(1:end-1,2:end)   + lon_d1(2:end,2:end) );

latc(2:nj, 2:ni) = 0.25 * ( ...
    lat_d1(1:end-1,1:end-1) + lat_d1(2:end,1:end-1) + ...
    lat_d1(1:end-1,2:end)   + lat_d1(2:end,2:end) );

mesh(lon,lat,mask,'edgecolor','none','facecolor','[0.8 0.8 0.8]','FaceAlpha','1')
hgt_mask=hgt_D1+mask_D1-mask_D1;
hgt_mask(isnan(hgt_mask))=-1;

u10_map_r(10,end-15)=20;
v10_map_r(10,end-15)=0;
plot(141.9326,39.3488,'r^','markersize',15,'linewidth',10)
plot(141.907,39.360,'bp','linewidth',10,'Markersize',15)
xp=[141.8 141.8 141.901 141.901];
yp=[39.5041 39.359 39.359 39.5041];
quiver(lon_map,lat_map,u10_map_r,v10_map_r,LL,'color','k','linewidth',4)
grid on

ax = gca;
ax.FontSize = 40;
axis([141.885, 141.985,39.322, 39.368])
daspect([1.3 1 1])
xticks([141.9 141.95 142])
xticklabels({'141.9','141.95','142'})
yticks([39.35])
yticklabels({'39.35'})

% % figure 6b
figure(5)

u10=u10_hourly_D2;
v10=v10_hourly_D2;

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
[xend,yend]=size(u10);

lon_map(1:9:xend,1:9:yend)=lon(1:9:xend,1:9:yend);
lon_map(5:9:xend,5:9:yend)=lon(5:9:xend,5:9:yend);
lat_map(1:9:xend,1:9:yend)=lat(1:9:xend,1:9:yend);
lat_map(5:9:xend,5:9:yend)=lat(5:9:xend,5:9:yend);
u10_map_r=zeros(size(xend,yend));
v10_map_r=zeros(size(xend,yend));
u10_map_r(1:9:xend,1:9:yend)=u10_hourly_D2(1:9:xend,1:9:yend);
u10_map_r(5:9:xend,5:9:yend)=u10_hourly_D2(5:9:xend,5:9:yend);
v10_map_r(1:9:xend,1:9:yend)=v10_hourly_D2(1:9:xend,1:9:yend);
v10_map_r(5:9:xend,5:9:yend)=v10_hourly_D2(5:9:xend,5:9:yend);

ws(1:75,end-170:end)=nan;

ws(ws<-1*10^-4)=-1*10^-4-10^-5;
ws(ws>1*10^-4)=1*10^-4+10^-5;

hold on

contourf(lon,lat,ws,linspace(-1*10^-4-10^-5,1*10^-4+10^-5,23),'linecolor','none');
colormap(mycmap)
colorbar
caxis([-1*10^-4-10^-5,1*10^-4+10^-5])
axis equal tight

mask1_d2 = ~isnan(mask_D2);  
mask1_d2 = logical(mask1_d2);


[nj, ni] = size(mask1_d2);
lonc = nan(nj+1, ni+1);
latc = nan(nj+1, ni+1);

% interior corners
lonc(2:nj, 2:ni) = 0.25 * ( ...
    lon_d2(1:end-1,1:end-1) + lon_d2(2:end,1:end-1) + ...
    lon_d2(1:end-1,2:end)   + lon_d2(2:end,2:end) );

latc(2:nj, 2:ni) = 0.25 * ( ...
    lat_d2(1:end-1,1:end-1) + lat_d2(2:end,1:end-1) + ...
    lat_d2(1:end-1,2:end)   + lat_d2(2:end,2:end) );

mesh(lon,lat,mask,'edgecolor','none','facecolor','[0.8 0.8 0.8]','FaceAlpha','1')

u10_map_r(10,end-15)=20;
v10_map_r(10,end-15)=0;
plot(141.9326,39.3488,'r^','markersize',15,'linewidth',10)
plot(141.907,39.360,'bp','linewidth',10,'Markersize',15)
xp=[141.8 141.8 141.901 141.901];
yp=[39.5041 39.359 39.359 39.5041];
quiver(lon_map,lat_map,u10_map_r,v10_map_r,LL,'color','k','linewidth',4)
grid on

ax = gca;
ax.FontSize = 40;
axis([141.885, 141.985,39.322, 39.368])
daspect([1.3 1 1])
xticks([141.9 141.95 142])
xticklabels({'141.9','141.95','142'})
yticks([39.35])
yticklabels({'39.35'})

% % figure 6a
figure(7)

u10=u10_hourly_D3;
v10=v10_hourly_D3;

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
[xend,yend]=size(u10);

lon_map(1:9:xend,1:9:yend)=lon(1:9:xend,1:9:yend);
lon_map(5:9:xend,5:9:yend)=lon(5:9:xend,5:9:yend);
lat_map(1:9:xend,1:9:yend)=lat(1:9:xend,1:9:yend);
lat_map(5:9:xend,5:9:yend)=lat(5:9:xend,5:9:yend);
u10_map_r=zeros(size(xend,yend));
v10_map_r=zeros(size(xend,yend));
u10_map_r(1:9:xend,1:9:yend)=u10_hourly_D3(1:9:xend,1:9:yend);
u10_map_r(5:9:xend,5:9:yend)=u10_hourly_D3(5:9:xend,5:9:yend);
v10_map_r(1:9:xend,1:9:yend)=v10_hourly_D3(1:9:xend,1:9:yend);
v10_map_r(5:9:xend,5:9:yend)=v10_hourly_D3(5:9:xend,5:9:yend);

ws(1:75,end-170:end)=nan;

ws(ws<-1*10^-4)=-1*10^-4-10^-5;
ws(ws>1*10^-4)=1*10^-4+10^-5;

hold on

contourf(lon,lat,ws,linspace(-1*10^-4-10^-5,1*10^-4+10^-5,23),'linecolor','none');
colormap(mycmap)
colorbar
caxis([-1*10^-4-10^-5,1*10^-4+10^-5])
axis equal tight

mask1_d3 = ~isnan(mask_D3);  
mask1_d3 = logical(mask1_d3);

[nj, ni] = size(mask1_d3);
lonc = nan(nj+1, ni+1);
latc = nan(nj+1, ni+1);

% interior corners
lonc(2:nj, 2:ni) = 0.25 * ( ...
    lon_d3(1:end-1,1:end-1) + lon_d3(2:end,1:end-1) + ...
    lon_d3(1:end-1,2:end)   + lon_d3(2:end,2:end) );

latc(2:nj, 2:ni) = 0.25 * ( ...
    lat_d3(1:end-1,1:end-1) + lat_d3(2:end,1:end-1) + ...
    lat_d3(1:end-1,2:end)   + lat_d3(2:end,2:end) );

mesh(lon,lat,mask,'edgecolor','none','facecolor','[0.8 0.8 0.8]','FaceAlpha','1')

u10_map_r(1:71,end-165:end)=nan;
v10_map_r(1:71,end-165:end)=nan;
u10_map_r(10,end-15)=20;
v10_map_r(10,end-15)=0;
u10_map_r(62,end-157)=10;
v10_map_r(62,end-157)=0;
lon_map(62,end-157)=lon(62,end-157);
lat_map(62,end-157)=lat(62,end-157);
plot(141.9326,39.3488,'r^','markersize',15,'linewidth',10)
plot(141.907,39.360,'bp','linewidth',10,'Markersize',15)
xp=[141.8 141.8 141.901 141.901];
yp=[39.5041 39.359 39.359 39.5041];
patch(xp,yp,'w')
quiver(lon_map,lat_map,u10_map_r,v10_map_r,LL,'color','k','linewidth',4)
text(141.885,39.362,'10 m/s','fontsize',50)
grid on

ax = gca;
ax.FontSize = 40;
axis([141.885, 141.985,39.322, 39.368])
daspect([1.3 1 1])
xticks([141.9 141.95 142])
xticklabels({'141.9','141.95','142'})
yticks([39.35])
yticklabels({'39.35'})
lon0 = 141.97;
lat0 = 39.362;
dlon_1km = 1 / (111.32 * cosd(lat0));
plot([lon0, lon0 + dlon_1km],[lat0, lat0],'m','linewidth',5);
