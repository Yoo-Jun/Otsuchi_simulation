clc; clear all; close all;

load('dataset_obs_wind_land')
load('dataset_obs_wind_sea')
load('dataset_simulated_wind_land')
load('dataset_simulated_wind_sea')
time_1=linspace(0,24*(((length(L_U10)-1)/6)/24),length(L_U10));
time_2=linspace(0,24*(((length(L_U10)-1)/6)/24),361);

figure(1)
% wind over the sea
subplot(2,1,1)
plot(time_1,S_U10,'k','linewidth',3)
hold on
plot(time_2,S_u10_hour,'r','linewidth',3)
xticks(linspace(0,360,31))
xticklabels({'','','','','','','','','','','','','' ...
    ,'','','','','','','','','','','','','','','','',''})
axis([0 time_1(end)-12-4 -10 10])
legend({'OBS','Simulation'},'Location','best')
grid on
ylabel('U10 (m/s)')
ax = gca;
ax.FontSize = 20;
mean(S_U10(end-12-4-6:end-12-4))
mean(S_U10(end-12-3:end-12-2))

subplot(2,1,2)
plot(time_1,S_V10,'k','linewidth',3)
hold on
plot(time_2,S_v10_hour,'r','linewidth',3)
xticks(linspace(0,480,41))
xticklabels({'0623','','0624','','0625','','0626','','0627','','0628','','0629' ...
    ,'','0630','','0701','','0702','','0703','','0704','','0705','','0706','','0707','','0708' ...
    ,'','0709','','0710','','0711','','0712','','0713'})
xtickangle(90)
axis([0 time_1(end)-12-4 -10 10])
grid on
ylabel('V10 (m/s)')
ax = gca;
ax.FontSize = 20; 

figure(2)
% wind over the land
subplot(2,1,1)
plot(time_1,L_U10,'k','linewidth',3)
hold on
plot(time_2,L_u10_hour,'r','linewidth',3)
xticks(linspace(0,360,31))
xticklabels({'','','','','','','','','','','','','' ...
    ,'','','','','','','','','','','','','','','','',''})
axis([0 time_1(end)-12-4 -10 10])
legend({'OBS','Simulation'},'Location','best')
grid on
ylabel('U10 (m/s)')
ax = gca;
ax.FontSize = 20; 

subplot(2,1,2)
plot(time_1,L_V10,'k','linewidth',3)
hold on
plot(time_2,L_v10_hour,'r','linewidth',3)
xticks(linspace(0,360,31))
xticklabels({'0623','','0624','','0625','','0626','','0627','','0628','','0629' ...
    ,'','0630','','0701','','0702','','0703','','0704','','0705','','0706','','0707',''})
xtickangle(90)
axis([0 time_1(end)-12-4 -10 10])
grid on
ylabel('V10 (m/s)')
ax = gca;
ax.FontSize = 20; 



