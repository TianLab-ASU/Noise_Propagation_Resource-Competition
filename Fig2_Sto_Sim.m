Rs_all=[0.1:0.1:0.7 0.75:0.05:1.1 1.2:0.1:2];

parfor i=1:length(Rs_all)
    Rs=Rs_all(i);
    StochInitial(Rs);
end

function StochInitial(Rs)

ICleft=[1 1 30 60];
ICright=[15 75 0 5];
IC0=[7 32 6 44];

% GFP // Gene 1
global km01 km1 kp1 dm1 dp1 cp1 Jm1 Jp1
km01=.4; km1=8; kp1=15; dm1=1; dp1=1; cp1=Rs; Jm1=40; Jp1=20;
% RFP // Gene 2
global km02 km2 kp2 dm2 dp2 cp2 Jm2 Jp2
km02=1; km2=25; kp2=30; dm2=1; dp2=1; cp2=1; Jm2=40; Jp2=2;
global n Kg Omega  
n=3; Kg=17; Omega=1.5;

[T,Y] = ode23s(@(t,y)ODE(t,y),[0 5000],IC0,[]);
IC(:)=ceil(Y(end,:)*Omega);
m1Mean=IC(1);
p1Mean=IC(2);
m2Mean=IC(3);
p2Mean=IC(4);

[T,Y] = ode23s(@(t,y)ODE(t,y),[0 5000],ICleft,[]);
IC_left(:)=ceil(Y(end,:)*Omega);
m1Ml=IC_left(1);
p1Ml=IC_left(2);
m2Ml=IC_left(3);
p2Ml=IC_left(4);

[T,Y] = ode23s(@(t,y)ODE(t,y),[0 5000],ICright,[]);
IC_right(:)=ceil(Y(end,:)*Omega);
m1Mr=IC_right(1);
p1Mr=IC_right(2);
m2Mr=IC_right(3);
p2Mr=IC_right(4);

% Stochastic Model:
t0=0;
gfp=[];
GFP=[];
rfp=[];
RFP=[];
t_tot=[];

for cell=1:1000
    [t,y] = ExactStochWendy(@(t,y) SDE(t,y),[t0 t0+200],IC);
    gfp=[gfp y(:,1).'];
    GFP=[GFP y(:,2).'];
    rfp=[rfp y(:,3).'];
    RFP=[RFP y(:,4).'];
    t_tot=[t_tot t.'];
    t0=t(end);
    IC=y(end,:);
    save(['Stoch_R1_' num2str(Rs) '_RC_wFluc.mat'])
end
 
end
%%

function dydt=ODE(~,y)
global km01 km1 kp1 dm1 dp1 cp1 Jm1 Jp1
global km02 km2 kp2 dm2 dp2 cp2 Jm2 Jp2
global n Kg    

m1=y(1);
p1=y(2);
m2=y(3);
p2=y(4);

R1=cp1;
R2=cp2*Kg^n/(p1^n+Kg^n);
PFm=1+R1/Jm1+R2/Jm2; % Resource Competition for mRNA
PFp=1+m1/Jp1+m2/Jp2; % Resource Competition for protein
 
dydt = [km01+km1*R1/PFm-dm1*m1;
    kp1*m1/PFp-dp1*p1;
    km02+km2*R2/PFm-dm2*m2;
    kp2*m2/PFp-dp2*p2;];
end
 

%%

function [t,X]=ExactStochWendy(DefineReactions,tspan, IC)
%Implementation of Gillespie Exact Stochastic Algorithm.
%BY WENDY THOMAS
% S = stoichiometry of C substrates in R reactions
% P = stoichiometry of products
% K = vector of reaction rates
% INITIALIZE by calling the passed function to define the reactions:
if ~isa(DefineReactions, 'function_handle') ...
        disp('pass a function handle with your reaction defs'); return; end

% if size(S) ~= size(P); disp('reaction defs inconsistent'); return;end
% if size(IC,1) ~= 1; IC = IC'; end  % force row vector
% if length(IC) ~= C; disp('ICs inconsistent with reactions'); return; end
% if size(K,2) ~= 1; K = K'; end  % force row vector
% if length(K) ~= R; disp ('parameters inconsistent with reactions');return;end
if length(tspan) ~= 2; disp('Time vector: [dt;maxtime] needed'); return;end
%set the uniform random number generator
rand('state',sum(100*clock));
nRC  = 0;  %reaction counter
t(1) = tspan(1);
maxT = tspan(2);
X = IC;  % initialize the matrix of chemicals
%------------------------------------------------------------ 
while (t(end)<maxT && any(X(end,:)))
    
    
    [S, P, K] = DefineReactions(t,X(end,:)); % defines reaction stoichiometry
    [R, C] = size(S);
    %step 1: Calculate a's (reaction rates given system state) 
    a = K;
    for r = 1:R
        for c = 1:C
            if S(r,c) == 1; a(r) = a(r)*X(end,c);
            elseif S(r,c) == 2
                a(r) = a(r)*X(end,c)*(X(end, c)-1)/2;
            elseif S(r,c) == 3
                a(r) = a(r)*X(end,c)*(X(end, c)-1)/2*(X(end,c)-2)/3;
            end
        end
    end
    a0 = sum(a); % a0 is the total rate of change of system
    if a0 == 0 % system can't change; finish and exit;
        X(end+1,:) = X(end,:);
        t =[t; maxT];
        break;
    end
    
    %Step 2: calculate tau and r using random number generators
    % determine time of next reaction:
    p1  = rand;  tau = (1/a0)*log(1/p1);
    % determine which next reaction is:
    p2 = rand;
    for r=1:R
        if (sum(a(1:r)) >= p2*a0); break; end
    end
    %Step 3: carry out the reaction
    t =[t; t(end) + tau];  % t is time array; add last entry to it.
    nRC = nRC + 1       ; % nRC is number of reactions so far.
    X(end+1,:)=X(end,:)-S(r,:)+P(r,:);
end %end of while (t(end)<tspan(2) & any(X))
end % end of function

%%

function [S, P, K] = SDE(~,y)
global km01 km1 kp1 dm1 dp1 cp1 Jm1 Jp1
global km02 km2 kp2 dm2 dp2 cp2 Jm2 Jp2
global n Kg Omega  

S = [  % How many of each chemical is used as substractes
    0 0 0 0% DNA  --> mRNA, GFP
    1 0 0 0% mRNA --> nothing, GFP
    1 0 0 0% mRNA --> P, GFP
    0 1 0 0% P --> nothing,  GFP
    0 0 0 0% DNA  --> mRNA, RFP
    0 0 1 0% mRNA --> nothing, RFP
    0 0 1 0% mRNA --> P, RFP
    0 0 0 1% P --> nothing, RFP
    ]  ;
P = [  % How many does each chemical species is made as products
    1 0 0 0% DNA  --> mRNA, GFP
    0 0 0 0% mRNA --> nothing, GFP
    1 1 0 0% mRNA --> P, GFP
    0 0 0 0% P --> nothing,  GFP
    0 0 1 0% DNA  --> mRNA, RFP
    0 0 0 0% mRNA --> nothing, RFP
    0 0 1 1% mRNA --> P, RFP
    0 0 0 0% P --> nothing, RFP
    ]  ;

m1=y(1);
p1=y(2);
m2=y(3);
p2=y(4);

R1=cp1;
R2=cp2*Kg^n/((p1/Omega)^n+Kg^n);
PFm=1+R1/Jm1+R2/Jm2; % Resource Competition for mRNA
PFp=1+m1/Omega/Jp1+m2/Omega/Jp2; % Resource Competition for protein

K  = [km01*Omega+km1*R1/PFm*Omega;
    dm1;
    kp1/PFp;
    dp1;
    km02*Omega+km2*R2/PFm*Omega;
    dm2;
    kp2/PFp;
    dp2;];
end