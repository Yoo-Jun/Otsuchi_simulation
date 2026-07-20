clc; clear all; close all;


nx=300;
ny=250;
dt=60*30;

load('dataset_grid_ocean.mat')
load('dataset_w_estimation.mat')


% box
S_end=double(XLAT_O(16+64,1));
N_end=double(XLAT_O(25+64,1));
W_end=double(XLON_O(1,40));
E_end=double(XLON_O(1,52));

% grid for figure
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

% divergence / curl 
dd(:,:,:)=Div_H;
cd(:,:,:)=Curl_H;

% estimation of w based on divergence
w_est_div=zeros(300,250);
for i=1:length(Div_H(1,1,:))-1
w_est_div(:,:,i+1)=w_est_div(:,:,i)+dd(:,:,i+1)*dt;
end

% time mean of estimated w by divergence and curl
mean_w_est_div(:,:)=trapz(w_est_div(:,:,end-12:end),3)*dt/(60*60*12);
mean_w_est_curl(:,:)=trapz(cd(:,:,end-12:end),3)*dt/(60*60*12);

A_lon=lon_sub;
A_lon(lon_sub<W_end)=nan;
A_lon(lon_sub>E_end)=nan;
A_lat=lat_sub;
A_lat(lat_sub<S_end)=nan;
A_lat(lat_sub>N_end)=nan;

AREA=double(A_lon+A_lat);
AREA(isnan(AREA)~=1)=0;
% isopycnal shoal
% by divervence
right1_A=mean(mean(mean_w_est_div+AREA',2,'omitnan'),1,'omitnan')*60*60*6
% by curl
right2_A=mean(mean(mean_w_est_curl+AREA',2,'omitnan'),1,'omitnan')*60*60*6

% figure grid
SFVERT1=interp2(lon_sub,lat_sub,mean_w_est_curl',lon_high,lat_high);
SFVERT2=interp2(lon_sub,lat_sub,mean_w_est_div',lon_high,lat_high);

mycmap=[0.0196    0.1882    0.3804
    0.0690    0.2935    0.5459
    0.1236    0.3907    0.6651
    0.1756    0.4771    0.7266
    0.2431    0.5590    0.7563
    0.3707    0.6529    0.8055
    0.5300    0.7481    0.8572
    0.6632    0.8216    0.8983
    0.7796    0.8786    0.9309
    0.8730    0.9275    0.9585
    0.9439    0.9655    0.9779
    0.9842    0.9572    0.9405
    0.9934    0.9029    0.8465
    0.9903    0.8257    0.7331
    0.9769    0.7258    0.5984
    0.9446    0.6116    0.4762
    0.8895    0.4831    0.3729
    0.8273    0.3493    0.2860
    0.7705    0.2026    0.2158
    0.6875    0.0854    0.1644
    0.5631    0.0213    0.1332
    0.4039         0    0.1216];


figure(1)
%Fig.10a
scale=2.2*10^-4;
scale2=scale*10^-2;
SFVERT1(SFVERT1>scale)=scale-scale2;
SFVERT1(SFVERT1<-scale)=-scale+scale2;
contourf(lon_high,lat_high,SFVERT1,linspace(-scale, scale,23),'color',[0.4 0.4 0.4])
hold on
colormap(mycmap)
caxis([-scale,scale])
mesh(lon_high,lat_high,Mtop2d_high,'edgecolor','none','facecolor','[0.5 0.5 0.5]')
colorbar
axis([141.885, 141.985,39.322, 39.368])
daspect([1.3 1 1])
xticks([141.9 141.95 142])
xticklabels({'','',''})
yticks([39.35])
yticklabels({''})
plot([W_end W_end],[S_end N_end],'k','linewidth',3)
plot([E_end E_end],[S_end N_end],'k','linewidth',3)
plot([W_end E_end],[S_end S_end],'k','linewidth',3)
plot([W_end E_end],[N_end N_end],'k','linewidth',3)
grid on
ax = gca;
ax.FontSize = 40;
lon0 = 141.97;
lat0 = 39.362;
dlon_1km = 1 / (111.32 * cosd(lat0));
plot([lon0, lon0 + dlon_1km],[lat0, lat0],'m','linewidth',5);

figure(2)
%Fig.10b
SFVERT2(SFVERT2>scale)=scale-scale2;
SFVERT2(SFVERT2<-scale)=-scale+scale2;
contourf(lon_high,lat_high,SFVERT2,linspace(-scale, scale,23),'color',[0.4 0.4 0.4])
hold on
colormap(mycmap)
caxis([-scale,scale])

mesh(lon_high,lat_high,Mtop2d_high,'edgecolor','none','facecolor','[0.5 0.5 0.5]')
colorbar
axis([141.885, 141.985,39.322, 39.368])
daspect([1.3 1 1])
xticks([141.9 141.95 142])
xticklabels({'','',''})
yticks([39.35])
yticklabels({''})
plot([W_end W_end],[S_end N_end],'k','linewidth',3)
plot([E_end E_end],[S_end N_end],'k','linewidth',3)
plot([W_end E_end],[S_end S_end],'k','linewidth',3)
plot([W_end E_end],[N_end N_end],'k','linewidth',3)
grid on
ax = gca;
ax.FontSize = 40;

