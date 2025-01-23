function SN=BifurcationPlot(k,x,y,cor)
[a,~]=find(abs(diff(k))>1);
index=[1 a' length(x)]; 
LineWidth=[0.5 2];
linewid = ["--" ,"-"];
for i=1:length(index)-1
   plot(x(index(i):index(i+1)),y(index(i):index(i+1)),'LineWidth',LineWidth(mod(i,2)+1),'color',cor,'LineStyle',linewid(mod(i,2)+1))
   hold on
end 
index=index(2:end-1); 
index=sort(index);
bpoint=find(diff(index)==2);

index=index(~ismember(index,index([bpoint,bpoint+1])));
SN=x(index);
plot(x(index),y(index),'o','MarkerEdgeColor',cor)
% plot(x(index(i+1):end),y(index(i+1):end),'LineWidth',4)