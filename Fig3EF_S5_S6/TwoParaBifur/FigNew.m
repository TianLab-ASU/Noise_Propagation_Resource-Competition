clc
clear
close all
Parameters
RC_Data=load('.\OneParaBifur\km01b.txt');  
hold on
BifurcationPlot(RC_Data(:,1),RC_Data(:,4),RC_Data(:,5)*Omega,'g')  
% xlim([0 1.5])
xlabel('GFP Inducer (I1)') 
ylabel('Protein Numbers') 
axis square
box on 
 