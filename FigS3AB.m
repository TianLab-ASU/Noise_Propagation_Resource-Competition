%%% In the bistable region.
close all
clear
clc
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
Ks_all=1:2:40;
for i=1:length(Ks_all)
    Ks=Ks_all(i);
    load(['.\Data\Stoch_Kg_' num2str(Ks) '_RC.mat'])
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

subplot(2,2,1)
plot(Ks_all,p1M,'g','linewidth',3)
hold on
plot(Ks_all,p2M,'r','linewidth',3)
hold on
xline(15)
hold on
xline(19)
ylabel('Protein Mean')
xlabel('Inhibition Threshold (K_g)')
axis square
box on 

subplot(2,2,2)
plot(Ks_all,gpN,'g','linewidth',3)
hold on
plot(Ks_all,rpN,'r','linewidth',3)
hold on
xline(15)
hold on
xline(19)
ylabel('Protein Noise')
xlabel('Inhibition Threshold (K_g)')
box on
axis square
savefig('FigS3ab')