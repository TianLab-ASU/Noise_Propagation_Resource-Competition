clc
clear
close all

file_list = dir('.\OneParaBifur\*.txt');
paraname = {file_list.name};
paraname = strrep(paraname, '.txt', '');
paraname = reshape(paraname,2,[]);
paraname = paraname';
paraname = paraname(:,1);
paraname = strrep(paraname, 'a', '');

cor=['r';'b'];
for i=3:length(file_list)
    Parameters
    subplot(4,5,floor((i-1)/2))
    RC_Data=load(['.\OneParaBifur\' file_list(i).name]);
    hold on 
    SN(i,:)=BifurcationPlot(RC_Data(:,1),RC_Data(:,4),RC_Data(:,5)*Omega,cor(mod(i,2)+1));

    RC_Data=load('.\OneParaBifur\I1a.txt');
    hold on 
    BifurcationPlot(RC_Data(:,1),RC_Data(:,4),RC_Data(:,5)*Omega,'k');
  

    xlim([0.5 1.5])
    title(paraname(floor((i+1)/2)))
    % xlabel('GFP Inducer (I1)')
    % ylabel('Protein Numbers')
    axis square
    box on 
end
    xlabel('GFP Inducer (I1)')
    ylabel('Protein Numbers')
%%
figure
subplot(2,2,3) 
SN1=SN(:,1);
SN1=reshape(SN1,2,[])/1.1294;
SN2=SN(:,2);
SN2=reshape(SN2,2,[])/0.8778;
[~, neworder]=sort(max(abs(SN1-1)));
neworder=neworder(1:end-1); 

hold on
b=bar(SN1(1,neworder));
b=bar(SN1(2,neworder));
b.BaseValue =1;
xticks(1:18) 
xticklabels(paraname(neworder))
box on
ylim([0.7 1.4])
 
subplot(2,2,4)
hold on
b=bar(SN2(1,neworder));
b=bar(SN2(2,neworder));
b.BaseValue =1;
box on
ylim([0.7 1.4])
xticks(1:18) 
xticklabels(paraname(neworder))
%%
figure
file_list = dir('.\TwoParaBifur\*.txt');
subplot(2,2,3)

cor=['r';'b']; 
for i=5:length(file_list)
    subplot(4,5,floor((i-1)/4))
    Parameters
    RC_Data=load(['.\TwoParaBifur\' file_list(i).name]);
    hold on 
    plot(RC_Data(:,1),RC_Data(:,2),cor(mod(ceil(i/ 2), 2) + 1),'LineWidth',1)

    RC_Data=load('.\TwoParaBifur\I1_2p_a1.txt' );
    hold on 
    plot(RC_Data(:,1),RC_Data(:,2),'--k','LineWidth',1)

    RC_Data=load('.\TwoParaBifur\I1_2p_a2.txt');
    hold on 
    plot(RC_Data(:,1),RC_Data(:,2),'--k','LineWidth',1)


    xlim([0 2])
    ylim([0 1])
    title(paraname(floor((i+3)/4)))
    axis square
    box on 
end
    xlabel('GFP Inducer (I1)')
ylabel('RFP Resource Competitivity (1/Jp2)') 