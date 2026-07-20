clc; clear all; close all;

load('dataset_scatter_s')

figure(1)
plot(sm_h,so,'or','linewidth',5)
corrcoef(sm_h,so)
grid on
hold on
plot(sm_u,so,'oc','linewidth',5)
corrcoef(sm_u,so)
xlabel('Simulated SSS (psu)')
ylabel('Observed SSS (psu)')
xticks(29:34)
axis([29 34.5 29 34.5])
pbaspect([1 1 1])
plot([29 34.5],[29 34.5],'k','linewidth',3)
plot([29 34.5],[29+0.5 34.5+0.5],'k:','linewidth',1.5)
plot([29 34.5],[29+1 34.5+1],'k:','linewidth',3)
plot([29 34.5],[29+1.5 34.5+1.5],'k:','linewidth',1.5)
plot([29 34.5],[29+2 34.5+2],'k:','linewidth',3)
plot([29 34.5],[29-0.5 34.5-0.5],'k:','linewidth',1.5)
plot([29 34.5],[29-1 34.5-1],'k:','linewidth',3)
plot([29 34.5],[29-1.5 34.5-1.5],'k:','linewidth',1.5)
plot([29 34.5],[29-2 34.5-2],'k:','linewidth',3)
legend('CaseH','CaseU')
ax = gca;
ax.FontSize = 40;


