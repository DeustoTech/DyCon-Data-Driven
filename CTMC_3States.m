%%
clear all 
close all
%%


figure('unit','norm','position',[0 0 0.6 1])
l1 = plot(xt,'o-','LineWidth',4)
yticks([1 2])
yticklabels({'Casa','Trabajo'})
ylim([0.75 2.25])
xlabel('Días')
ylabel('Estado')
l1.Parent.FontSize = 22;
title('Evolución temporal del estado')
saveas(l1.Parent.Parent,'T12.png');

%%
ixa = xt(1);
count = 0;
Ns = 0;
iter = 1;
for ixc = xt(2:end)'
   iter = iter + 1;
   if ixa == ixc
       count = count + 1;
   else
       count = 0;
   end
   Ns(iter) = count;
   ixa = ixc;
end

data = [xt Ns'];
%%
data = [data(2:end,1) data(1:end-1,1) data(2:end,2)];
%%

data1 = data(data(:,1) == 1,:);
data11 = data1(data1(:,2) == 1,:);
data12 = data1(data1(:,2) == 2,:);

P12 = length(data12(:,1))/length(data1(:,1));
P11 = length(data11(:,1))/length(data1(:,1));

%%

data2 = data(data(:,1) == 2,:);
data21 = data2(data2(:,2) == 1,:);
data22 = data2(data2(:,2) == 2,:);
%

P21 = length(data21(:,1))/length(data2(:,1));
P22 = length(data22(:,1))/length(data2(:,1));

%%

A = [P11 P12;
     P21 P22];
Px0 = [1 0]
Pxt(1,:) = Px0;
%
xtex = [];
for it = 2:50
    ixt = double(data(it,1) == [1 2]);
    Pxt(it,:) = A'*Pxt(it-1,:)'; 
    xtex(it-1) = randsample(2,1,true,A*ixt');
end

fig = figure('unit','norm','position',[0 0 0.6 1]);
subplot(2,1,1)
l1 = plot(xt,'o-','LineWidth',4);
yticks([1 2])
yticklabels({'Casa','Trabajo'})
ylim([0.75 2.25])
xlabel('Días')
ylabel('Estado')
l1.Parent.FontSize = 22;
title('Evolución temporal real del estado')
subplot(2,1,2)
l1 = plot(xtex,'o-','LineWidth',4);
yticks([1 2])
yticklabels({'Casa','Trabajo'})
ylim([0.75 2.25])
xlabel('Días')
ylabel('Estado')
l1.Parent.FontSize = 22;
title('Evolución temporal estimada del estado')

saveas(fig,'TEstimada.png')


%%
P12N = [];
P11N = [];
P21N = [];
P22N = [];

for iN = 0:10
    data1iN = data1(data1(:,3) == iN,:);
    
    N11 = length(data1iN(data1iN(:,2) ==1,2));
    N12 = length(data1iN(data1iN(:,2) ==2,2));

    P11N(iN+1) = (N11+1)/(N11+N12+2);
    P12N(iN+1) = (N12+1)/(N11+N12+2);

    data2iN = data2(data2(:,3) == iN,:);
    N21 = length(data2iN(data2iN(:,2) ==1,2));
    N22 = length(data2iN(data2iN(:,2) ==2,2));
    
    P21N(iN+1) = (N21+1)/(N21+N22+2);
    P22N(iN+1) = (N22+1)/(N21+N22+2);
end

%%
xtex = [];
for it = 2:68
%
    ixt = double(data(it,1) == [1 2]);
    iN =  data(it,3) + 1;
    A = [P11N(iN) P12N(iN);
         P21N(iN) P22N(iN)];
    %
    xtex(it-1) = median(randsample(2,9,true,A'*ixt'));
end

fig = figure('unit','norm','position',[0 0 0.6 1]);
subplot(2,1,1)
l1 = plot(xt(2:end-1),'o-','LineWidth',4);
yticks([1 2])
yticklabels({'Casa','Trabajo'})
ylim([0.75 2.25])
xlabel('Días')
ylabel('Estado')
l1.Parent.FontSize = 22;
title('Evolución temporal real del estado')
subplot(2,1,2)
l1 = plot(xtex,'o-','LineWidth',4);
yticks([1 2])
yticklabels({'Casa','Trabajo'})
ylim([0.75 2.25])
xlabel('Días')
ylabel('Estado')
l1.Parent.FontSize = 22;
title('Evolución temporal estimada del estado')
saveas(fig,'Memory.png')