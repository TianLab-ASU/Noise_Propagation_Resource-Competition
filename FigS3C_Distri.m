close all
clc
clear
Ks_all=[13 17 23];
for i=1:length(Ks_all)
    figure
    Kv=Ks_all(i);
    load(['.\Data\Stoch_Kg_' num2str(Kv) '_RC.mat'])
    subplot(1,2,1)
    h=histogram(GFP,0:220,'Normalization','probability','DisplayStyle','stairs')
    area(h.BinEdges(1:end-1),h.Values,'linewidth',2)
    ylabel('Probability')
    xlabel('GFP Level')
    title(['K_g=' num2str(Kv)])
    if i==3
        ylim([0 0.08])
    else
        ylim([0 0.035])
    end
    xlim([0 200])
    axis square
    box on
    subplot(1,2,2)
    h=histogram(RFP,0:150,'Normalization','probability','DisplayStyle','stairs')
    area(h.BinEdges(1:end-1),h.Values,'linewidth',2)
    ylabel('Probability')
    xlabel('RFP Level')
    title(['K_g=' num2str(Kv)])
    if i==3
        ylim([0 0.08])
    else
        ylim([0 0.035])
    end
    xlim([0 150])
    axis square
    box on 
end