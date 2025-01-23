close all
clc
clear
Kchange=[13 17 23];
for i=1:length(Kchange)
    figure
    Ks=Kchange(i);
    load(['.\Data\Stoch_Kg_' num2str(Ks) '_RC.mat'])
    ymax=max([GFP RFP]);
    ymax=50*ceil(ymax/50)-50;
    xmax=500;
    subplot(2,1,1)
    yline(g_locs(gp1))
    hold on
    yline(g_locs(gp2))
    hold on
    plot(t_tot(1:199831),GFP(1:199831),'g','linewidth',1.5)
    title(['K_g=' num2str(Ks)])
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
    title(['K_g=' num2str(Ks)])
    ylabel('RFP Expression')
    xlabel('Time')
    xlim([0 xmax])
    ylim([0 150])
    box on
    savefig(['SupFig2c' num2str(i) 'traj'])
end