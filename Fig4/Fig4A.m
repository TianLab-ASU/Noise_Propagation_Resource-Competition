close all
clc
clear
Rs_all=[0.1:0.1:0.7 0.75:0.05:1.1 1.2:0.1:2];
m1M=[];
p1M=[];
m2M=[];
p2M=[];
gmN=[];
gpN=[];
rmN=[];
rpN=[];
mCV=[];
pCV=[];
mCC=[];
pCC=[];
for i=1:length(Rs_all)
    Rs=Rs_all(i);
    load(['.\Data\Stoch_R1_' num2str(Rs) '_RC_woFluc.mat'])
    % Means
    m1M=[m1M mean(gfp)];
    p1M=[p1M mean(GFP)];
    m2M=[m2M mean(rfp)];
    p2M=[p2M mean(RFP)];
    % Noise
    gmN=[gmN std(gfp)/mean(gfp)];
    gpN=[gpN std(GFP)/mean(GFP)];
    rmN=[rmN std(rfp)/mean(rfp)];
    rpN=[rpN std(RFP)/mean(RFP)];
end
%%
% hold on
% yyaxis left
% plot(Rs_all,p1M,'--','color',[0.47,0.67,0.19],'linewidth',3)
% ylabel('GFP Mean')
% yyaxis right
% plot(Rs_all,p2M,'--','color',[0.64,0.08,0.18],'linewidth',3)
% ylabel('RFP Mean')
% hold on
% xline(0.6)
% hold on
% xline(1.4)
% xlabel('GFP Copy Number (c_{p_1})')
% xlim([0 2])
% set(gca,'xscale','log')
% axis square
% box on

figure
yyaxis left
plot(Rs_all,gpN,'g','linewidth',3)
ylabel('GFP Noise')
yyaxis right
plot(Rs_all,rpN,'r','linewidth',3)
ylabel('RFP Noise')
xlabel('GFP Copy Number (c_{p_1})')
hold on
xline(0.6)
hold on
xline(1.4)
xlim([0 2])
set(gca,'xscale','log')
% hold on
% xline(0.87)
% hold on
% xline(1.13)
% xlim([0 2])
axis square
box on 