close all
clc
clear
Rs_all=[0.5 0.8 1 1.2 1.5];
for i=1:length(Rs_all)
    figure
    Rv=Rs_all(i);
    load(['.\Data\Stoch_R1_' num2str(Rv) '_RC_wFluc.mat'])
    subplot(1,2,1)
    h=histogram(GFP,0:220,'Normalization','probability','DisplayStyle','stairs')
    area(h.BinEdges(1:end-1),h.Values,'linewidth',2)
    ylabel('Probability')
    xlabel('GFP Level')
    title(['c_{p_1}=' num2str(Rv)])
    if i==1
        ylim([0 0.12])
    elseif i==2
        ylim([0 0.08])
    elseif i==3
        ylim([0 0.035])
    elseif i==4
        ylim([0 0.03])
    else
        ylim([0 0.04])
    end
    xlim([0 200])
    axis square
    box on
    subplot(1,2,2)
    h=histogram(RFP,0:150,'Normalization','probability','DisplayStyle','stairs')
    area(h.BinEdges(1:end-1),h.Values,'linewidth',2)
    ylabel('Probability')
    xlabel('RFP Level')
    title(['c_{p_1}=' num2str(Rv)])
    if i==1
        ylim([0 0.12])
    elseif i==2
        ylim([0 0.08])
    elseif i==3
        ylim([0 0.035])
    elseif i==4
        ylim([0 0.03])
    else
        ylim([0 0.04])
    end
    xlim([0 150])
    axis square
    box on
    savefig(['Fig2d' num2str(i) 'hist'])
end