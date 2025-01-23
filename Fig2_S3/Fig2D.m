close all
clc
clear
Rs_all=[0.5 0.8 1 1.2 1.5];
startindex=35800000;
for i=1:length(Rs_all)
    figure
    Rv=Rs_all(i);
    load(['.\Data\Stoch_R1_' num2str(Rv) '_RC_wFluc.mat'])
    ymax=max([GFP RFP]);
    ymax=50*ceil(ymax/50)-50;
    xmax=500;
    subplot(2,1,1)
    yline(g_locs(gp1))
    hold on
    yline(g_locs(gp2))
    hold on
    if i==2 | i==4
        plot(t_tot(startindex:(startindex+199831)),GFP(startindex:(startindex+199831)),'g','linewidth',1.5)
        xlim([t_tot(startindex) t_tot(startindex)+xmax])
    else
        plot(t_tot(1:199831),GFP(1:199831),'g','linewidth',1.5)
        xlim([0 xmax])
    end
    title(['c_{p_1}=' num2str(Rv)])
    ylabel('GFP Expression')
    xlabel('Time')
    ylim([0 200])
    box on
    subplot(2,1,2)
    yline(r_locs(rp1))
    hold on
    yline(r_locs(rp2))
    hold on
    if i==2 | i==4
        plot(t_tot(startindex:(startindex+199831)),RFP(startindex:(startindex+199831)),'r','linewidth',1.5)
        xlim([t_tot(startindex) t_tot(startindex)+xmax])
    else
        plot(t_tot(1:199831),RFP(1:199831),'r','linewidth',1.5)
        xlim([0 xmax])
    end
    title(['c_{p_1}=' num2str(Rv)])
    ylabel('RFP Expression')
    xlabel('Time')
    ylim([0 150])
    box on
    savefig(['Fig2d' num2str(i) 'traj'])
end