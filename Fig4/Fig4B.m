clc
clear
figure
Rv=1;
load(['.\Data\Stoch_R1_' num2str(Rv) '_RC_woFluc.mat'])
ymax=max([GFP RFP]);
ymax=50*ceil(ymax/50)-50;
xmax=500;
subplot(2,1,1)
yline(g_locs(gp1))
hold on
yline(g_locs(gp2))
hold on
plot(t_tot(1:199831),GFP(1:199831),'g','linewidth',1.5)
title(['c_{p_1}=' num2str(Rv)])
ylabel('GFP Expression')
xlabel('Time')
xlim([0 xmax])
ylim([0 200])
box on
subplot(2,1,2)
yline(r_locs(rp1))
hold on
yline(r_locs(rp2))
hold on
plot(t_tot(1:199831),RFP(1:199831),'r','linewidth',1.5)
title(['c_{p_1}=' num2str(Rv)])
ylabel('RFP Expression')
xlabel('Time')
xlim([0 xmax])
ylim([0 150])
box on 