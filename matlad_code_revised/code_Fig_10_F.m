clc; clear all; close all;

nx=300;
ny=250;
surfacelayerthickness=2;
dt=30*60;


load('dataset_salt_budget_sprime.mat')
load('dataset_salt_budget_SSS.mat')
load('dataset_salt_budget_w.mat')
load('dataset_grid_ocean.mat')

%boxed region
S_end=double(XLAT_O(16+64,1));
N_end=double(XLAT_O(25+64,1));
W_end=double(XLON_O(1,40));
E_end=double(XLON_O(1,52));

%salinity tendency 
dsdt_H=(SSS_H(:,:,2:end)-SSS_H(:,:,1:end-1))./(dt);% Case_H salinity tendency
dsdt_U=(SSS_U(:,:,2:end)-SSS_U(:,:,1:end-1))./(dt);% Case_U salinity tendency

%eq.6 time-integration
left(:,:)=trapz(dsdt_H(:,:,end-12:end)-dsdt_U(:,:,end-12:end),3)*dt;
right_1(:,:)=-(trapz((W2m_H(:,:,end-12:end)-W2m_U(:,:,end-12:end)).*Sprime_U(:,:,end-12:end),3)*dt)./surfacelayerthickness;
right_2(:,:)=-(trapz(W2m_U(:,:,end-12:end).*(Sprime_H(:,:,end-12:end)-Sprime_U(:,:,end-12:end)),3)*dt)./surfacelayerthickness;

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

SFVERT1=interp2(lon_sub,lat_sub,left',lon_high,lat_high);%left in eq.6
SFVERT2=interp2(lon_sub,lat_sub,(right_1'),lon_high,lat_high);%right 1 in eq.6
SFVERT3=interp2(lon_sub,lat_sub,(right_2'),lon_high,lat_high);%right 2 in eq.6
SFVERT4=interp2(lon_sub,lat_sub,left'-(right_1'+right_2'),lon_high,lat_high);%left - right in eq.6


%boxed region average
A_lon=lon_sub;
A_lon(lon_sub<W_end)=nan;
A_lon(lon_sub>E_end)=nan;
A_lat=lat_sub;
A_lat(lat_sub<S_end)=nan;
A_lat(lat_sub>N_end)=nan;
AREA=double(A_lon+A_lat);
AREA(isnan(AREA)~=1)=0;
left_A=mean(mean(left+AREA',2,'omitnan'),1,'omitnan');
right1_A=mean(mean(right_1+AREA',2,'omitnan'),1,'omitnan');
right2_A=mean(mean(right_2+AREA',2,'omitnan'),1,'omitnan');

%figures code
figure(1)
%Fig.8a
scale=1.1;
SFVERT1(SFVERT1>scale)=scale-0.0001;
SFVERT1(SFVERT1<-scale)=-scale+0.0001;
contourf(lon_high,lat_high,SFVERT1,linspace(-scale, scale,23),'color',[0.4 0.4 0.4])
hold on
colormap(jet(22))
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

% cross-section line in Fig.9
xxx=XLON_O(1,1:60);
yyy(20:60)=XLAT_O(64+20,1);
yyy(1:20)=nan;
for i=21:40
    yyy(i)=XLAT_O(64+20-round((40-i)/1.2),1);   
end
plot([xxx(21),xxx(40),xxx(60)],[yyy(21),yyy(40),yyy(60)],':w','linewidth',10)
lon0 = 141.97;
lat0 = 39.362;
dlon_1km = 1 / (111.32 * cosd(lat0));
plot([lon0, lon0 + dlon_1km],[lat0, lat0],'m','linewidth',5);

figure(2)
%Fig.8b
SFVERT2(SFVERT2>scale)=scale-0.01;
SFVERT2(SFVERT2<-scale)=-scale+0.01;
contourf(lon_high,lat_high,SFVERT2,linspace(-scale, scale,23),'color',[0.4 0.4 0.4])
hold on
colormap(jet(22))
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

figure(3)
%Fig.8c
SFVERT3(SFVERT3>scale)=scale-0.01;
SFVERT3(SFVERT3<-scale)=-scale+0.01;
contourf(lon_high,lat_high,SFVERT3,linspace(-scale, scale,23),'color',[0.4 0.4 0.4])
hold on
colormap(jet(22))
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

figure(4)
%Fig.8d
SFVERT4(SFVERT4>scale)=scale-0.01;
SFVERT4(SFVERT4<-scale)=-scale+0.01;
contourf(lon_high,lat_high,SFVERT4,linspace(-scale, scale,23),'color',[0.4 0.4 0.4])
hold on
colormap(jet(22))
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