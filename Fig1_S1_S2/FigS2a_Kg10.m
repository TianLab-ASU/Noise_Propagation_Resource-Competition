clc
clear

syms H11 H22 H33 H44
syms H21 H23 H32 H41 H43
syms tau1 tau2 tau3 tau4
syms m1 m2 p1 p2

M=-[H11/tau1     0         0         0;
    H21/tau2  H22/tau2     0         0;
    0         H32/tau3  H33/tau3     0;
    0            0      H43/tau4  H44/tau4;
    ];

D=[2/m1/tau1     0          0          0 ;
    0        2/p1/tau2      0          0 ;
    0           0       2/m2/tau3      0 ;
    0           0           0      2/p2/tau4  ;];


noise=FlucDis(M,D);
noise1=simplify(noise(1,1));
noise2=simplify(noise(2,2));
noise3=simplify(noise(3,3));
noise4=simplify(noise(4,4));

%%

syms H11 H22 H33 H44
syms H21 H23 H32 H43
syms tau1 tau2 tau3 tau4
syms m1 m2 p1 p2

noise1_AE=1/m1;

noise12=noise1_AE*H21^2*1/tau2/(1/tau1+1/tau2);
noise22=1/p1;
noise2_AE=noise22+noise12;


noise123=noise12*(1+1/tau2/(1/tau1+1/tau3))*H32^2*1/tau3/(1/tau3+1/tau2);
noise23=noise22*H32^2/tau3/(1/tau3+1/tau2);
noise33=1/m2;
noise3_AE=noise123+noise23+noise33;

noise1234=noise123*(1+1/tau3/(1/tau2+1/tau4)+1/tau3/(1/tau2+1/tau4)*1/tau2/(1/tau1+1/tau2+1/tau3)*(1/tau3+1/tau2)/(1/tau1+1/tau4))*H43^2*1/tau4/(1/tau3+1/tau4); 
noise234 =noise23*(1+1/tau3/(1/tau2+1/tau4))*H43^2*1/tau4/(1/tau3+1/tau4);
noise34  =noise33                     *H43^2*1/tau4/(1/tau3+1/tau4);
noise44=1/p2;
noise4_AE=noise44+noise1234+noise234+noise34;

simplify(noise1-noise1_AE)
simplify(noise2-noise2_AE)
simplify(noise3-noise3_AE)
simplify(noise4-noise4_AE)

%%
Parameters
Kg=10; %% FigS2a
syms cp1
R1=cp1;
R2=cp2*Kg^n/(p1^n+Kg^n);

g1=km01*Omega+km1*R1*Omega;
f1=m1*dm1;

g2=kp1*m1;
f2=p1*dp1;

g3=km02*Omega+km2*R2*Omega;
f3=m2*dm2;

g4=kp2*m2;
f4=p2*dp2;

eqns=[g1-f1==0,g2-f2==0,g3-f3==0,g4-f4==0];
S = solve(eqns,[m1 p1 m2 p2]);
mean_All=[S.m1 S.p1 S.m2 S.p2]; 

Matrix=[f1/g1,f2/g2,f3/g3,f4/g4];

J=jacobian(Matrix,[m1,p1,m2,p2]);

H=J./repmat(Matrix.',1,4).*repmat(([m1 p1 m2 p2]),4,1);

H=subs(H,{m1,p1,m2,p2},{mean_All(1),mean_All(2),mean_All(3),mean_All(4)});

% to calcute overall noise
RFP_noise_total=subs(noise4,{H11, H22, H33, H44, H21, H23, H32, H43},...
    {H(1,1),H(2,2),H(3,3),H(4,4),H(2,1),H(2,3),H(3,2),H(4,3)});
RFP_noise_total=subs(RFP_noise_total,{tau1, tau2, tau3, tau4},{1/dm1,1/dp1,1/dm2,1/dp2});
RFP_noise_total=subs(RFP_noise_total,{m1,p1,m2,p2},{mean_All(1),mean_All(2),mean_All(3),mean_All(4)});
RFP_noise_total=simplify(sqrt(RFP_noise_total))

% 
% 1/p2+noise1234+noise234+noise34;

RFP_noise1234_total=subs(noise1234,{H11, H22, H33, H44, H21, H23, H32, H43},...
    {H(1,1),H(2,2),H(3,3),H(4,4),H(2,1),H(2,3),H(3,2),H(4,3)});
RFP_noise1234_total=subs(RFP_noise1234_total,{tau1, tau2, tau3, tau4},{1/dm1,1/dp1,1/dm2,1/dp2});
RFP_noise1234_total=subs(RFP_noise1234_total,{m1,p1,m2,p2},{mean_All(1),mean_All(2),mean_All(3),mean_All(4)});
RFP_noise1234_total=simplify(sqrt(RFP_noise1234_total))

RFP_noise234_total=subs(noise234,{H11, H22, H33, H44, H21, H23, H32, H43},...
    {H(1,1),H(2,2),H(3,3),H(4,4),H(2,1),H(2,3),H(3,2),H(4,3)});
RFP_noise234_total=subs(RFP_noise234_total,{tau1, tau2, tau3, tau4},{1/dm1,1/dp1,1/dm2,1/dp2});
RFP_noise234_total=subs(RFP_noise234_total,{m1,p1,m2,p2},{mean_All(1),mean_All(2),mean_All(3),mean_All(4)});
RFP_noise234_total=simplify(sqrt(RFP_noise234_total))

RFP_noise34_total=subs(noise34,{H11, H22, H33, H44, H21, H23, H32, H43},...
    {H(1,1),H(2,2),H(3,3),H(4,4),H(2,1),H(2,3),H(3,2),H(4,3)});
RFP_noise34_total=subs(RFP_noise34_total,{tau1, tau2, tau3, tau4},{1/dm1,1/dp1,1/dm2,1/dp2});
RFP_noise34_total=subs(RFP_noise34_total,{m1,p1,m2,p2},{mean_All(1),mean_All(2),mean_All(3),mean_All(4)});
RFP_noise34_total=simplify(sqrt(RFP_noise34_total))

RFP_noise4_total=subs(1/p2,{H11, H22, H33, H44, H21, H23, H32, H43},...
    {H(1,1),H(2,2),H(3,3),H(4,4),H(2,1),H(2,3),H(3,2),H(4,3)});
RFP_noise4_total=subs(RFP_noise4_total,{tau1, tau2, tau3, tau4},{1/dm1,1/dp1,1/dm2,1/dp2});
RFP_noise4_total=subs(RFP_noise4_total,{m1,p1,m2,p2},{mean_All(1),mean_All(2),mean_All(3),mean_All(4)});
RFP_noise4_total=simplify(sqrt(RFP_noise4_total))

rfp_noise_total=subs(noise3,{H11, H22, H33, H44, H21, H23, H32, H43},...
    {H(1,1),H(2,2),H(3,3),H(4,4),H(2,1),H(2,3),H(3,2),H(4,3)});
rfp_noise_total=subs(rfp_noise_total,{tau1, tau2, tau3, tau4},{1/dm1,1/dp1,1/dm2,1/dp2});
rfp_noise_total=subs(rfp_noise_total,{m1,p1,m2,p2},{mean_All(1),mean_All(2),mean_All(3),mean_All(4)});
rfp_noise_total=simplify(sqrt(rfp_noise_total))

gfp_noise_total=subs(noise1,{H11, H22, H33, H44, H21, H23, H32, H43},...
    {H(1,1),H(2,2),H(3,3),H(4,4),H(2,1),H(2,3),H(3,2),H(4,3)});
gfp_noise_total=subs(gfp_noise_total,{tau1, tau2, tau3, tau4},{1/dm1,1/dp1,1/dm2,1/dp2});
gfp_noise_total=subs(gfp_noise_total,{m1,p1,m2,p2},{mean_All(1),mean_All(2),mean_All(3),mean_All(4)});
gfp_noise_total=simplify(sqrt(gfp_noise_total))
   
GFP_noise_total=subs(noise2,{H11, H22, H33, H44, H21, H23, H32, H43},...
    {H(1,1),H(2,2),H(3,3),H(4,4),H(2,1),H(2,3),H(3,2),H(4,3)});
GFP_noise_total=subs(GFP_noise_total,{tau1, tau2, tau3, tau4},{1/dm1,1/dp1,1/dm2,1/dp2});
GFP_noise_total=subs(GFP_noise_total,{m1,p1,m2,p2},{mean_All(1),mean_All(2),mean_All(3),mean_All(4)});
GFP_noise_total=simplify(sqrt(GFP_noise_total))

GFP_noise12_total=subs(noise12,{H11, H22, H33, H44, H21, H23, H32, H43},...
    {H(1,1),H(2,2),H(3,3),H(4,4),H(2,1),H(2,3),H(3,2),H(4,3)});
GFP_noise12_total=subs(GFP_noise12_total,{tau1, tau2, tau3, tau4},{1/dm1,1/dp1,1/dm2,1/dp2});
GFP_noise12_total=subs(GFP_noise12_total,{m1,p1,m2,p2},{mean_All(1),mean_All(2),mean_All(3),mean_All(4)});
GFP_noise12_total=simplify(sqrt(GFP_noise12_total))

GFP_noise2_total=subs(1/p1,{H11, H22, H33, H44, H21, H23, H32, H43},...
    {H(1,1),H(2,2),H(3,3),H(4,4),H(2,1),H(2,3),H(3,2),H(4,3)});
GFP_noise2_total=subs(GFP_noise2_total,{tau1, tau2, tau3, tau4},{1/dm1,1/dp1,1/dm2,1/dp2});
GFP_noise2_total=subs(GFP_noise2_total,{m1,p1,m2,p2},{mean_All(1),mean_All(2),mean_All(3),mean_All(4)});
GFP_noise2_total=simplify(sqrt(GFP_noise2_total))
  
 
cp1_ALL=10.^(-5:.01:.3);
 
RFP_noise_ALL_num =(double(subs(RFP_noise_total, cp1, cp1_ALL))); 
RFP_noise1234_num =(double(subs(RFP_noise1234_total, cp1, cp1_ALL))); 
RFP_noise234_num =(double(subs(RFP_noise234_total, cp1, cp1_ALL))); 
RFP_noise34_num =(double(subs(RFP_noise34_total, cp1, cp1_ALL))); 
RFP_noise4_num =(double(subs(RFP_noise4_total, cp1, cp1_ALL))); 
GFP_noise_ALL_num =(double(subs(GFP_noise_total, cp1, cp1_ALL))); 
GFP_noise12_num =(double(subs(GFP_noise12_total, cp1, cp1_ALL))); 
GFP_noise2_num =(double(subs(GFP_noise2_total, cp1, cp1_ALL))); 
rfp_noise_ALL_num =(double(subs(rfp_noise_total, cp1, cp1_ALL))); 
gfp_noise_ALL_num =(double(subs(gfp_noise_total, cp1, cp1_ALL))); 

GFP_Mean_num=(double(subs(mean_All(2), cp1, cp1_ALL))); 
RFP_Mean_num=(double(subs(mean_All(4), cp1, cp1_ALL))); 
gfp_Mean_num=(double(subs(mean_All(1), cp1, cp1_ALL))); 
rfp_Mean_num=(double(subs(mean_All(3), cp1, cp1_ALL))); 
%%

subplot(2,2,1)
hold on
plot(cp1_ALL,RFP_noise_ALL_num,'color','r','linewidth',2)
plot(cp1_ALL,RFP_noise1234_num,'color',[0.93,0.69,0.13],'linewidth',2)
plot(cp1_ALL,RFP_noise234_num,'color',[0.07,0.62,1.00],'linewidth',2)
plot(cp1_ALL,RFP_noise34_num,'color',[0.94,0.50,0.77],'linewidth',2)
plot(cp1_ALL,RFP_noise4_num,'color',[0.49,0.18,0.56],'linewidth',2)
xlabel('GFP Copy Number (c_{p_1})')
ylabel('RFP Noise Level')
% set(gca,'yscale','log')
legend('\eta_{p_2 total}','\eta_{{p_2}<-{m_2}<-{p_1}<-{m_1}}','\eta_{{p_2}<-{m_2}<-{p_1}}','\eta_{{p_2}<-{m_2}}','\eta_{{p_2}<-{p_2}}','location','northeast')
box on
axis square
set(gca,'xscale','log')
xlim([0.01 2])
ylim([0 1.2])
% savefig('SupFig1b')

 