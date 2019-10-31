clear all;close all

ImportGolden_users;
out = out(out.Speed == 0.0,:);
time = datenum(out.Date);
wk   = weekday(out.Date);


%% Order by time 
[t, indexord] = sort(time);
%t = t - t(1);
x   = out.Lat(indexord);
y   = out.Long(indexord);
wk = wk(indexord);
%%
X = [x y];
[idx,C] = kmeans(X,3,'Distance','cityblock');
 
%%
data = [t C(idx,:) idx X];
%%
%
RepIndex = diff(t) == 0;
data(RepIndex,:) = [];
t(RepIndex) = [];
idx(RepIndex) = [];
X(RepIndex,:) = [];


%%
newtspan = t(1):(1/2):t(end);
newidx = interp1(t,idx,newtspan,'nearest');
newX   = interp1(t,X,newtspan,'nearest');
data = [newtspan' C(newidx,:) newidx' newX ];
%%
subplot(1,2,1)
plot([newX(:,1)-mean(newX(:,1)) newX(:,2)-mean(newX(:,2))])
title('Real Model')
subplot(1,2,2)
plot([data(:,2)-mean(data(:,2)) data(:,3)-mean(data(:,3))])
title('Reduce Model')
%%
figure('unit','norm','pos',[0 0 1 1])
%%
clf
gscatter(X(:,1),X(:,2),idx,'rgbyck',[],[],'off')
hold on
plot(C(:,1),C(:,2),'ko','MarkerSize',15)

% l1 = line(x(1),y(1),'Marker','.','MarkerSize',40,'Color','r');

l2 = line(x(1),y(1),'Marker','.','MarkerSize',15,'Color','b');

grid on
xlabel time
ylabel latitud
zlabel longitud
t1 = title("t = "+t(1));

for it = 2:1:length(newtspan)
%     l1.XData = data(it,2);
%     l1.YData = data(it,3);
    
    l2.XData = data(it,5);
    l2.YData = data(it,6);
    %l1.ZData = data(it,3);

    pause(0.1)
%     l1.Color = 'r';
    t1.String = datestr(t(it), 'ddd-dd-mmm HH:MM:SS');
    pause(0.1)
%     l1.Color = 'k';

end

xt = data(:,4)