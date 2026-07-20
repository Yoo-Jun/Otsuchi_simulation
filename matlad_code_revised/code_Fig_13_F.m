clc; clear all; close all

addpath('H:\Otsuchi\code_9\Figure_codes\dataset')
%
load('dataset_psd_windstress_D4')
load('dataset_psd_W')

dt=5*60;
latitude = 39.35; % Latitude in degrees
omega = 7.2921e-5; % Earth's angular velocity in rad/s
f_cor = 2 * omega * sind(latitude); % Coriolis frequency in Hz
inertial_frequency = 1/(2*pi/f_cor)*3600; % Convert to cycles/hour

% wind sress divergence & curl D4
txx=zeros(size(Ftx_h));
txy=zeros(size(Ftx_h));
tyx=zeros(size(Ftx_h));
tyy=zeros(size(Ftx_h));
txx(2:end-1,2:end-1,:)=((Ftx_h(2:end-1,2:end-1,:))-(Ftx_h(1:end-2,2:end-1,:)))/100;
txy(2:end-1,2:end-1,:)=((Ftx_h(2:end-1,2:end-1,:))-(Ftx_h(2:end-1,1:end-2,:)))/100;
tyx(2:end-1,2:end-1,:)=((Fty_h(2:end-1,2:end-1,:))-(Fty_h(1:end-2,2:end-1,:)))/100;
tyy(2:end-1,2:end-1,:)=((Fty_h(2:end-1,2:end-1,:))-(Fty_h(2:end-1,1:end-2,:)))/100;
div=(txx(2:end-1,2:end-1,:)+tyy(2:end-1,2:end-1,:))/1025;
curl=(tyx(2:end-1,2:end-1,:)-txy(2:end-1,2:end-1,:))*(f_cor/1025);

% Ekman response
f2W=f_cor*f_cor*W;

% Adjustment response
dw2dt2 = zeros(13,10,4178);
dw2dt2(:, :, 2:4178-1) = (W(:, :, 3:4178) - 2 * W(:, :, 2:4178-1) + W(:, :, 1:4178-2)) / dt^2;

% Divergence forcing
divdt = zeros(13,10,4178);
divdt(:,:,2:4178-1)=(div(:,:,3:4178)-div(:,:,1:4178-2))/(2*dt);



fs = 1 / dt;
n_days = 14.5; % Number of days
n_samples_per_day = fs * 24 * 60 * 60; % Samples per day
total_samples = round(n_samples_per_day * n_days); % Total samples (rounded)
time = (0:total_samples-1) / fs;

window = hamming(1024); % Window for pwelch
noverlap=512;
nfft=2048;

sigw(:)=sum(sum(dw2dt2(:,:,2:end-1),2,'omitnan'),1,'omitnan');
sigf(:)=sum(sum(f2W(:,:,2:end-1),2,'omitnan'),1,'omitnan');
sigd(:)=sum(sum(divdt(:,:,2:end-1),2,'omitnan'),1,'omitnan');
sigc(:)=sum(sum(curl(:,:,2:end-1),2,'omitnan'),1,'omitnan');





[pxxws fq pxxwc]= pwelch(sigw, window, noverlap, nfft, fs,'confidencelevel',0.99);
[pxxds fq pxxdc]= pwelch(sigd, window, noverlap, nfft, fs,'confidencelevel',0.99);
[pxxcs fq pxxcc]= pwelch(sigc, window, noverlap, nfft, fs,'confidencelevel',0.99);
[pxxfs fq pxxfc]= pwelch(sigf, window, noverlap, nfft, fs,'confidencelevel',0.99);

f_cph = fq*3600; % Convert Hz to cph
f_cph(f_cph<1/(24*5))=nan; %mask for timescale longer than 5 day

% Figure 11
figure(11)
loglog(f_cph, pxxws, 'k','linewidth',2.5); 
hold on
loglog(f_cph, pxxfs, 'm','linewidth',2.5); 
loglog(f_cph, pxxds, 'b','linewidth',4);
loglog(f_cph, pxxcs, 'r','linewidth',4);

% confidence interval
for i=2:length(f_cph)
patch([f_cph(i),f_cph(i),f_cph(i-1),f_cph(i-1)],[pxxwc(i,1),pxxwc(i,2),pxxwc(i-1,2),pxxwc(i-1,1)],'k','Facealpha',0.3,'Edgecolor','none')
end
hold on
for i=2:length(f_cph)
patch([f_cph(i),f_cph(i),f_cph(i-1),f_cph(i-1)],[pxxdc(i,1),pxxdc(i,2),pxxdc(i-1,2),pxxdc(i-1,1)],'b','Facealpha',0.3,'Edgecolor','none')
end
for i=2:length(f_cph)
patch([f_cph(i),f_cph(i),f_cph(i-1),f_cph(i-1)],[pxxcc(i,1),pxxcc(i,2),pxxcc(i-1,2),pxxcc(i-1,1)],'r','Facealpha',0.3,'Edgecolor','none')
end
for i=2:length(f_cph)
patch([f_cph(i),f_cph(i),f_cph(i-1),f_cph(i-1)],[pxxfc(i,1),pxxfc(i,2),pxxfc(i-1,2),pxxfc(i-1,1)],'m','Facealpha',0.3,'Edgecolor','none')
end

xlabel('Frequency (cph)');
ylabel('Power Spectral Density ({m^2}{s^{-5}})');

ax = gca; 
grid(ax, 'on');  
ax.GridColor = [0, 0, 0]; 
ax.GridLineStyle = '-'; 
ax.GridAlpha = 1; 
ax.LineWidth = 1;

% Annotate key frequencies using vertical lines
p24_freq = 1/24; % 1 per 24 h
xline(p24_freq, 'k--', 'LineWidth', 1.5);
p48_freq = 1/48; % 1 per 48 h
xline(p48_freq, 'k--', 'LineWidth', 1.5);
p72_freq = 1/72; % 1 per 72 h
xline(p72_freq, 'k--', 'LineWidth', 1.5);
p12_freq = 1/12; % 1 per 12 h
xline(p12_freq, 'k--', 'LineWidth', 1.5); 
p6_freq = 1/6;  % 1 per 6 h
xline(p6_freq, 'k--', 'LineWidth', 1.5); 
p2_freq = 1/2; % 1 per 2 h
xline(p2_freq , 'k--', 'LineWidth', 1.5); 
p1_freq = 1/1; % 1 per 1 h
xline(p1_freq , 'k--', 'LineWidth', 1.5); 

xline(inertial_frequency/2, 'color',[0.8500 0.3250 0.0980], 'LineWidth', 1.5); % 1 per inertial period
xline(inertial_frequency, 'color',[0.8500 0.3250 0.0980], 'LineWidth', 1.5); % 1 per inertial period
xline(inertial_frequency*2, 'color',[0.8500 0.3250 0.0980], 'LineWidth', 1.5); % 1 per half inertial period

xlim([10^-2.5 10^1])
ylim([10^-23 10^-9])
legend('$\frac{\partial^2<w_H>}{\partial t^2}$','{$f^2$}{$<w_H>$}','{$\frac{1}{\rho}$}{$\frac{\partial}{\partial t}$}{$(\nabla\bullet\tau)$}',...
    '{$\frac{f}{\rho}$}{$(\nabla\times\tau)$}','Interpreter', 'latex','Location','eastoutside','fontsize',30);
ax = gca;
ax.FontSize = 30;
pbaspect([1 1 1])


