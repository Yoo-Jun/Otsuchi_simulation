clc; clear all; close all;
%need m_map

load('dataset_grid_ocean.mat')
load('dataset_OI_coef.mat')
load('dataset_obs_TS.mat')
load('dataset_obs_station.mat')
load('dataset_coastline_river.mat')
obs_period=10;

%20140626;obs_period=1
%20140728;obs_period=2
%20141030;obs_period=3
%20150216;obs_period=4
%20150317;obs_period=5
%20150423;obs_period=6
%20150527;obs_period=7
%20150612;obs_period=8
%20150624;obs_period=9
%20150707;obs_period=10
%20150722;obs_period=11
%20150819;obs_period=12
%20150903;obs_period=13
%20151023;obs_period=14
%20151106;obs_period=15

south_end=39.30;
west_end=141.87;
north_end=39.4;
east_end=142.03;

long_s(1:42)=sta_14(:,1);
lat_s(1:42)=sta_14(:,2);

Gx=double(XLON_O);
Gy=double(XLAT_O);
Py=double((XLAT_O(:,:)-XLAT_O(1,1))*111.04);
Px=double((XLON_O(:,:)-XLON_O(1,1))*111.04*cos(deg2rad(39.345)));

gX=reshape(Gx,[],1);
gY=reshape(Gy,[],1);
pX=reshape(Px,[],1);
pY=reshape(Py,[],1);

sraw=obs_sal(:,obs_period);

lonraw=sta_14(:,1);
latraw=sta_14(:,2);


dist_xraw=(lonraw-Gx(1,1))*111.04*cos(deg2rad(39.345));
dist_yraw=(latraw-Gy(1,1))*111.04;

% sraw(1)=nan;
s=rmmissing(sraw);


lonobs=rmmissing(lonraw+sraw-sraw);
latobs=rmmissing(latraw+sraw-sraw);
dist_x=rmmissing(dist_xraw+sraw-sraw);
dist_y=rmmissing(dist_yraw+sraw-sraw);

for i=1:length(s)
    for j=1:length(s)
xobs(i,j)=(lonobs(i)-lonobs(j))*111.04*cos(deg2rad(39.345));
yobs(i,j)=(latobs(i)-latobs(j))*111.04;
    end
end

cf=exp(A(1)*xobs.^2+A(2)*yobs.^2+A(3)*xobs.*yobs);
SNL=(1-exp(A(4)))/exp(A(4));
SNLm=eye(size(xobs))*SNL;
Cm=cf+SNLm;


mend=1;

Cok=ones(length(s)+mend,length(s)+mend,1);
Cok(1:length(s),1:length(s),1)=Cm;
Cok(length(s)+1:end,length(s)+1:end,1)=0;


Cstaok=ones(length(s)+mend,1);
for i=1:length(gX)
for k=1:length(s)
    dxx(k,1)=(lonobs(k)-gX(i))*111.04*cos(deg2rad(39.345));
    dyy(k,1)=(latobs(k)-gY(i))*111.04;
end
Csta=exp(A(1)*dxx.^2+A(2)*dyy.^2+A(3)*dxx.*dyy);
Cstaok(1:length(s))=Csta;

w_ukq=Cok(1:length(s)+mend,1:length(s)+mend)\Cstaok(1:length(s)+mend);

oi_uks(i)=w_ukq(1:end-mend)'*s;
ei_uk(i)=1-((Cstaok(1:length(s)+mend)'/Cok(1:length(s)+mend,1:length(s)+mend))*Cstaok(1:length(s)+mend));
end

li=griddata(lonobs,latobs,s,Gx,Gy,'linear');
OI_uks(:,:,obs_period)=reshape(oi_uks,[],300);
OI_UKS(:,:,obs_period)=OI_uks(:,:,obs_period)+li-li;
ei_uks(:,:,obs_period)=reshape(ei_uk,[],300);
E_UK(:,:,obs_period)=ei_uks(:,:,obs_period)+li-li;

figure(obs_period)
m_proj('Mercator','long',[141.894, 141.985],'lat',[39.322, 39.368]);
[cs, h]=m_contourf(Gx,Gy,OI_UKS(:,:,obs_period),linspace(30,33,7),'color',[0.5 0.5 0.5]);
OBS(:,:)=OI_UKS(:,:,obs_period);
mcmap=[0.2422    0.1504    0.2
    0.2422    0.1504    0.6603
    0.2810    0.3228    0.9579
    0.1786    0.5289    0.9682
    0.2161    0.7843    0.5923
    0.6720    0.7793    0.2227
    0.9970    0.7659    0.2199
    0.9769    0.9839    0.0805];
colormap(mcmap)
caxis([29.5,33.5])
hold on
colorbar('Ticks',[30 31 32 33],...
         'TickLabels',{'30','31','32','33'},'fontsize',50)
m_grid('box','fancy','tickdir','in','tickstyle','dd','xtick',[141.9 141.95 142],'ytick',[39.3 39.35 39.4],'fontsize',40);
 daspect([1 1 1])
hold on
m_plot(sta_14(:,1),sta_14(:,2),'ko','linewidth',5,'Markersize',5)
m_patch(coast_m(:,2),coast_m(:,1),[0.5 0.5 0.5],'facealpha',1,'EdgeColor','k','linewidth',2)
m_patch(coast_island(59:73,2),coast_island(59:73,1),[0.5 0.5 0.5],'facealpha',1,'EdgeColor','k','linewidth',2)
m_patch(coast_island(594:645,2),coast_island(594:645,1),[0.5 0.5 0.5],'facealpha',1,'EdgeColor','k','linewidth',2)
m_plot(kotsuchi(1,:)+0.0001,kotsuchi(2,:),'w','linewidth',10)
m_plot(kotsuchi(1,:),kotsuchi(2,:)+0.0001,'w','linewidth',10)
m_plot(kotsuchi(1,:)+0.0001,kotsuchi(2,:)+0.0006,'k','linewidth',2.2)
m_plot(kotsuchi(1,:)-0.0001,kotsuchi(2,:)-0.0005,'k','linewidth',2.2)
m_plot([141.9083 141.9086],[39.35545 39.35545],'w','linewidth',3)
m_plot(otsuchi(1,:),otsuchi(2,:),'w','linewidth',15)
m_plot(otsuchi(1,:)+0.0007,otsuchi(2,:)+0.0005,'k','linewidth',2.2)
m_plot(otsuchi(1,:)-0.0007,otsuchi(2,:)-0.0005,'k','linewidth',2.2)
lon0 = 141.97;
lat0 = 39.362;
dlon_1km = 1 / (111.32 * cosd(lat0));
m_plot([lon0, lon0 + dlon_1km],[lat0, lat0],'m','linewidth',5);

clear all;
%need m_map
addpath('H:\Otsuchi\code_9\Figure_codes\dataset','H:\Otsuchi\m_map1.4\m_map')
%
load('dataset_grid_ocean.mat')
load('dataset_OI_coef_iso.mat')
load('dataset_obs_TS.mat')
load('dataset_obs_station.mat')
load('dataset_coastline_river.mat')
obs_period=10;

%20140626;obs_period=1
%20140728;obs_period=2
%20141030;obs_period=3
%20150216;obs_period=4
%20150317;obs_period=5
%20150423;obs_period=6
%20150527;obs_period=7
%20150612;obs_period=8
%20150624;obs_period=9
%20150707;obs_period=10
%20150722;obs_period=11
%20150819;obs_period=12
%20150903;obs_period=13
%20151023;obs_period=14
%20151106;obs_period=15

south_end=39.30;
west_end=141.87;
north_end=39.4;
east_end=142.03;

long_s(1:42)=sta_14(:,1);
lat_s(1:42)=sta_14(:,2);

Gx=double(XLON_O);
Gy=double(XLAT_O);
Py=double((XLAT_O(:,:)-XLAT_O(1,1))*111.04);
Px=double((XLON_O(:,:)-XLON_O(1,1))*111.04*cos(deg2rad(39.345)));

gX=reshape(Gx,[],1);
gY=reshape(Gy,[],1);
pX=reshape(Px,[],1);
pY=reshape(Py,[],1);

sraw=obs_sal(:,obs_period);

lonraw=sta_14(:,1);
latraw=sta_14(:,2);


dist_xraw=(lonraw-Gx(1,1))*111.04*cos(deg2rad(39.345));
dist_yraw=(latraw-Gy(1,1))*111.04;

% sraw(1)=nan;
s=rmmissing(sraw);


lonobs=rmmissing(lonraw+sraw-sraw);
latobs=rmmissing(latraw+sraw-sraw);
dist_x=rmmissing(dist_xraw+sraw-sraw);
dist_y=rmmissing(dist_yraw+sraw-sraw);

for i=1:length(s)
    for j=1:length(s)
xobs(i,j)=(lonobs(i)-lonobs(j))*111.04*cos(deg2rad(39.345));
yobs(i,j)=(latobs(i)-latobs(j))*111.04;
    end
end

cf=exp(A(1)*(xobs.^2+yobs.^2));
SNL=(1-exp(A(2)))/exp(A(2));
SNLm=eye(size(xobs))*SNL;
Cm=cf+SNLm;

mend=1;

Cok=ones(length(s)+mend,length(s)+mend,1);
Cok(1:length(s),1:length(s),1)=Cm;
Cok(length(s)+1:end,length(s)+1:end,1)=0;


Cstaok=ones(length(s)+mend,1);
for i=1:length(gX)
for k=1:length(s)
    dxx(k,1)=(lonobs(k)-gX(i))*111.04*cos(deg2rad(39.345));
    dyy(k,1)=(latobs(k)-gY(i))*111.04;
end
Csta=exp(A(1)*(dxx.^2+dyy.^2));
Cstaok(1:length(s))=Csta;

w_ukq=Cok(1:length(s)+mend,1:length(s)+mend)\Cstaok(1:length(s)+mend);

oi_uks(i)=w_ukq(1:end-mend)'*s;
ei_uk(i)=1-((Cstaok(1:length(s)+mend)'/Cok(1:length(s)+mend,1:length(s)+mend))*Cstaok(1:length(s)+mend));
end

li=griddata(lonobs,latobs,s,Gx,Gy,'linear');
OI_uks(:,:,obs_period)=reshape(oi_uks,[],300);
OI_UKS(:,:,obs_period)=OI_uks(:,:,obs_period)+li-li;
ei_uks(:,:,obs_period)=reshape(ei_uk,[],300);
E_UK(:,:,obs_period)=ei_uks(:,:,obs_period)+li-li;

figure(obs_period+1)
m_proj('Mercator','long',[141.894, 141.985],'lat',[39.322, 39.368]);
[cs, h]=m_contourf(Gx,Gy,OI_UKS(:,:,obs_period),linspace(30,33,7),'color',[0.5 0.5 0.5]);
OBS(:,:)=OI_UKS(:,:,obs_period);
mcmap=[0.2422    0.1504    0.2
    0.2422    0.1504    0.6603
    0.2810    0.3228    0.9579
    0.1786    0.5289    0.9682
    0.2161    0.7843    0.5923
    0.6720    0.7793    0.2227
    0.9970    0.7659    0.2199
    0.9769    0.9839    0.0805];
colormap(mcmap)
caxis([29.5,33.5])
hold on
colorbar('Ticks',[30 31 32 33],...
         'TickLabels',{'30','31','32','33'},'fontsize',50)
m_grid('box','fancy','tickdir','in','tickstyle','dd','xtick',[141.9 141.95 142],'ytick',[39.3 39.35 39.4],'fontsize',40);
 daspect([1 1 1])
hold on
m_plot(sta_14(:,1),sta_14(:,2),'ko','linewidth',5,'Markersize',5)
m_patch(coast_m(:,2),coast_m(:,1),[0.5 0.5 0.5],'facealpha',1,'EdgeColor','k','linewidth',2)
m_patch(coast_island(59:73,2),coast_island(59:73,1),[0.5 0.5 0.5],'facealpha',1,'EdgeColor','k','linewidth',2)
m_patch(coast_island(594:645,2),coast_island(594:645,1),[0.5 0.5 0.5],'facealpha',1,'EdgeColor','k','linewidth',2)
m_plot(kotsuchi(1,:)+0.0001,kotsuchi(2,:),'w','linewidth',10)
m_plot(kotsuchi(1,:),kotsuchi(2,:)+0.0001,'w','linewidth',10)
m_plot(kotsuchi(1,:)+0.0001,kotsuchi(2,:)+0.0006,'k','linewidth',2.2)
m_plot(kotsuchi(1,:)-0.0001,kotsuchi(2,:)-0.0005,'k','linewidth',2.2)
m_plot([141.9083 141.9086],[39.35545 39.35545],'w','linewidth',3)
m_plot(otsuchi(1,:),otsuchi(2,:),'w','linewidth',15)
m_plot(otsuchi(1,:)+0.0007,otsuchi(2,:)+0.0005,'k','linewidth',2.2)
m_plot(otsuchi(1,:)-0.0007,otsuchi(2,:)-0.0005,'k','linewidth',2.2)


