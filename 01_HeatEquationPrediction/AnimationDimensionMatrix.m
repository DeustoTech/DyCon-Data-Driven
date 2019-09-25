close all
fig = figure;
ax  = axes('Parent',fig);
%%
i = 5;
[A ,xms,yms] = Nbymatrix(i);
isurf  = surf(xms,yms,A/i^2,'Parent',ax);

%% 

%%
ititle = title(ax,"dim(A) ="+i+"x"+i);
view(ax,-20,60)

for i = 10:1:80
   [A ,xms,yms] = Nbymatrix(i);
   
   isurf.XData = xms;
   isurf.YData = yms;
   isurf.ZData = A/i^2;

   ititle.String = "dim(A) ="+i+"x"+i;
   pause(0.1)

end