clc
clear
close all
Parameters
RC_Data=load('RC_Bifur.txt'); 
subplot(2,2,2)
hold on
BifurcationPlot(RC_Data(:,1),RC_Data(:,4),RC_Data(:,6)*Omega,'g')
BifurcationPlot(RC_Data(:,1),RC_Data(:,4),RC_Data(:,8)*Omega,'r')

xline(1)
xlim([0 1.5])
xlabel('GFP Inducer (I1)') 
ylabel('Protein Numbers') 
axis square
box on 
%%
% clc
clear
subplot(2,2,3)
Parameters
% p2 or m2 nullcline with p1 as input
p1=0:0.01:70;
R1=cp1;
R2=cp2*Kg^n./(p1.^n+Kg^n);
PFm=1+R1/Jm1+R2/Jm2;
m1=(km01+km1*R1./PFm)/dm1;
m2=(km02+km2*R2./PFm)/dm2;
PFp=1+m1/Jp1+m2/Jp2;
p2=kp2*m2./PFp/dp2;
plot(p1*Omega,p2*Omega,'color',[0.92,0.58,0.29],'LineWidth',3) 
hold on
clc
clear
Parameters
% p1 nullcline with m2 as input 
R1=cp1;
R2=0:0.001:3;
PFm=1+R1/Jm1+R2/Jm2; 
m1=(km01+km1*R1./PFm)/dm1;
m2=(km02+km2*R2./PFm)/dm2;
PFp=1+m1/Jp1+m2/Jp2;
p1=kp1*m1./PFp/dp1;
p2=kp2*m2./PFp/dp2;
plot(p1*Omega,p2*Omega,'color',[0.00,0.45,0.74],'LineWidth',3) 
hold on
plot(10.25*Omega,52.85*Omega,'.','MarkerEdgeColor','k','MarkerSize',30)
hold on
plot(23.7*Omega,43.65*Omega,'o','MarkerEdgeColor','k','MarkerSize',8)
hold on
plot(53.75*Omega,23*Omega,'.','MarkerEdgeColor','k','MarkerSize',30)
xlim([0 60*Omega])
ylim([10*Omega 60*Omega])
xlabel('GFP Numbers')
ylabel('RFP Numbers')
legend('dp_2/dt = 0','dp_1/dt = 0')
box on
axis square 

%%
% %%
subplot(2,2,4)
clc
clear
Parameters
RC_Data=load('RC_Bifur_Jp2_1.txt');
plot(RC_Data(:,1),RC_Data(:,2)/Jp2,'k','LineWidth',3) 
RC_Data=load('RC_Bifur_Jp2_2.txt'); 
hold on
plot(RC_Data(:,1),RC_Data(:,2)/Jp2,'k','LineWidth',3) 
xlabel('GFP Copy Number (cp_1)') 
ylabel('RFP Resource Competitivity (1/Jp2)') 
hold on
yline(0.5)
ylim([0 1])
xlim([0 2])
axis square
box on 