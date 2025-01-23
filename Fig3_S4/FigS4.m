clc
clear
close all
Parameters
 
RC_Data=load('RC_Bifur_K.txt'); 
BifurcationPlot(RC_Data(:,1),RC_Data(:,4),RC_Data(:,6)*Omega,'g')
BifurcationPlot(RC_Data(:,1),RC_Data(:,4),RC_Data(:,8)*Omega,'r')
xlim([0 40])
xlabel('Inhibition Threshold (K_g)') 
ylabel('Protein Expression Level')
axis square
box on
savefig('SupFig3a')

