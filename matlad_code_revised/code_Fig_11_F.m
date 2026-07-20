clc; clear all; close all;

nx=300;
ny=250;
nz=80;
load('dataset_vertical_snapshot.mat')
load('dataset_grid_ocean.mat')

% East & West end of box in Fig.8, dashed line in Fig.9
W_end=double(XLON_O(1,40));
E_end=double(XLON_O(1,52));

% grid for figure
dz(1:20)=4;
dz(21:40)=3;
dz(41:60)=2;
dz(61:80)=1;
dep=zeros(1,nz);
for i=1:79
dep(i)=-sum(dz(i:end))+dz(i)/2;
end
lon_sub=XLON_O(1,1:108)-XLON_O(1,1);
[gd,gx]=meshgrid(dep(40:80),(double(XLON_O(1,1:60))-XLON_O(1,1))*111.04*cos(deg2rad(39.345))*10,60);
[gddh,gxxh]=meshgrid(linspace(dep(40),dep(80),161), ...
    linspace(0,double(XLON_O(1,60)-XLON_O(1,1))*111.04*cos(deg2rad(39.345))*10,60*3));

WWbay_H=interp2(double(gx)',gd',Wvert_H',gxxh',gddh');%Case_H vertical velocity
UUbay_H=interp2(double(gx)',gd',Uvert_H',gxxh',gddh');%Case_H horizontal velocity

WWbay_U=interp2(double(gx)',gd',Wvert_U',gxxh',gddh');%Case_U vertical velocity
UUbay_U=interp2(double(gx)',gd',Uvert_U',gxxh',gddh');%Case_U horizontal velocity

WWbay_A=interp2(double(gx)',gd',Wvert_H'-Wvert_U',gxxh',gddh');%Case_H-Case_U vertical velocity
UUbay_A=interp2(double(gx)',gd',Uvert_H'-Uvert_U',gxxh',gddh');%Case_H-Case_U horizontal velocity

%subsample for vector
spacex=10;
spacey=10;

uvech=squeeze(UUbay_H(1:spacex:end,1:spacey:end));
wvech=squeeze(WWbay_H(1:spacex:end,1:spacey:end))*200;

uvecu=squeeze(UUbay_U(1:spacex:end,1:spacey:end));
wvecu=squeeze(WWbay_U(1:spacex:end,1:spacey:end))*200;

uveca=squeeze(UUbay_A(1:spacex:end,1:spacey:end));
wveca=squeeze(WWbay_A(1:spacex:end,1:spacey:end))*200;

% scale vector
uvech(8,8)=0.1;
wvech(8,8)=0;
wvech(9,8)=0.1;
uvech(9,8)=0;

uvecu(8,8)=0.1;
wvecu(8,8)=0;
wvecu(9,8)=0.1;
uvecu(9,8)=0;

uveca(8,8)=0.1;
wveca(8,8)=0;
wveca(9,8)=0.1;
uveca(9,8)=0;

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
%Fig.9a

maskvert=double(isnan(WWbay_H));
maskvert(maskvert==0)=nan;
maskvert(maskvert==1)=0;
maskvert(60:110,65:100)=nan;
mesh(gxxh',gddh',maskvert,'edgecolor','none','facecolor','[0.5 0.5 0.5]')
view([0 0 1])
scale=2.2*10^-4;
WWbay_H(WWbay_H<-scale)=-scale+scale*10^-2;
WWbay_H(WWbay_H>scale)=scale-scale*10^-2;
hold on
contourf(gxxh',gddh',WWbay_H,linspace(-scale,scale,23))
plot([(W_end-XLON_O(1,1))*111.04*cos(deg2rad(39.345))*10 (W_end-XLON_O(1,1))*111.04*cos(deg2rad(39.345))*10], [-60 0],':k','linewidth',3)
plot([(E_end-XLON_O(1,1))*111.04*cos(deg2rad(39.345))*10 (E_end-XLON_O(1,1))*111.04*cos(deg2rad(39.345))*10], [-60 0],':k','linewidth',3)
colormap(mycmap)
hold on
quiver(gxxh(1:spacex:end,1:spacey:end)',gddh(1:spacex:end,1:spacey:end)',uvech,wvech,6.5,'w','linewidth',5)
quiver(gxxh(1:spacex:end,1:spacey:end)',gddh(1:spacex:end,1:spacey:end)',uvech,wvech,6.5,'k','linewidth',2.5)
text(23,-37,'0.1 m/s','fontsize',40)
text(24,-25,'5 x 10^-^4','fontsize',40)
text(28,-29,'m/s','fontsize',40)
xticks([22, 38.7587, 50.6825])
xticklabels({'141.9','141.9212','141.9351'})
grid on
caxis([-scale,scale])
axis([15 58 -50 1])
pbaspect([1.2 1 1])
ax = gca;
ax.FontSize = 40;
c=colorbar;

plot([39.75 49.75],[-45, -45],'m','linewidth',5);

figure(2)
%Fig.9b

maskvert=double(isnan(WWbay_U));
maskvert(maskvert==0)=nan;
maskvert(maskvert==1)=0;
maskvert(60:110,65:100)=nan;
mesh(gxxh',gddh',maskvert,'edgecolor','none','facecolor','[0.5 0.5 0.5]')
view([0 0 1])
scale=2.2*10^-4;
WWbay_U(WWbay_U<-scale)=-scale+scale*10^-2;
WWbay_U(WWbay_U>scale)=scale-scale*10^-2;
hold on
contourf(gxxh',gddh',WWbay_U,linspace(-scale,scale,23))
plot([(W_end-XLON_O(1,1))*111.04*cos(deg2rad(39.345))*10 (W_end-XLON_O(1,1))*111.04*cos(deg2rad(39.345))*10], [-60 0],':k','linewidth',3)
plot([(E_end-XLON_O(1,1))*111.04*cos(deg2rad(39.345))*10 (E_end-XLON_O(1,1))*111.04*cos(deg2rad(39.345))*10], [-60 0],':k','linewidth',3)
colormap(mycmap)
hold on
quiver(gxxh(1:spacex:end,1:spacey:end)',gddh(1:spacex:end,1:spacey:end)',uvecu,wvecu,6.61,'w','linewidth',5)
quiver(gxxh(1:spacex:end,1:spacey:end)',gddh(1:spacex:end,1:spacey:end)',uvecu,wvecu,6.61,'k','linewidth',2.5)
text(23,-37,'0.1 m/s','fontsize',40)
text(24,-25,'5 x 10^-^4','fontsize',40)
text(28,-29,'m/s','fontsize',40)
xticks([22, 38.7587, 50.6825])
xticklabels({'141.9','141.9212','141.9351'})
grid on
caxis([-scale,scale])
axis([15 58 -50 1])
pbaspect([1.2 1 1])
ax = gca;
ax.FontSize = 40;
c=colorbar;


figure(3)
%Fig.9c
maskvert=double(isnan(WWbay_A));
maskvert(maskvert==0)=nan;
maskvert(maskvert==1)=0;
maskvert(60:110,65:100)=nan;
mesh(gxxh',gddh',maskvert,'edgecolor','none','facecolor','[0.5 0.5 0.5]')
view([0 0 1])
scale=2.2*10^-4;
WWbay_A(WWbay_A<-scale)=-scale+scale*10^-2;
WWbay_A(WWbay_A>scale)=scale-scale*10^-2;
hold on
contourf(gxxh',gddh',WWbay_A,linspace(-scale,scale,23))
plot([(W_end-XLON_O(1,1))*111.04*cos(deg2rad(39.345))*10 (W_end-XLON_O(1,1))*111.04*cos(deg2rad(39.345))*10], [-60 0],':k','linewidth',3)
plot([(E_end-XLON_O(1,1))*111.04*cos(deg2rad(39.345))*10 (E_end-XLON_O(1,1))*111.04*cos(deg2rad(39.345))*10], [-60 0],':k','linewidth',3)
colormap(mycmap)
quiver(gxxh(1:spacex:end,1:spacey:end)',gddh(1:spacex:end,1:spacey:end)',uveca,wveca,3.07,'w','linewidth',5)
quiver(gxxh(1:spacex:end,1:spacey:end)',gddh(1:spacex:end,1:spacey:end)',uveca,wveca,3.07,'k','linewidth',2.5)
text(23,-37,'0.1 m/s','fontsize',40)
text(24,-25,'5 x 10^-^4','fontsize',40)
text(28,-29,'m/s','fontsize',40)
xticks([22, 38.7587, 50.6825])
xticklabels({'141.9','141.9212','141.9351'})
grid on
caxis([-scale,scale])
axis([15 58 -50 1])
pbaspect([1.2 1 1])
ax = gca;
ax.FontSize = 40;
c=colorbar;