%%
%%
clear all
% 
load('/home/djoroya/Documentos/GitHub/Data/DarioData.mat')
load('/home/djoroya/Documentos/GitHub/Data/tspan.mat')
% load('/home/djoroya/Documentos/GitHub/DyConData/code/DarioData/data/DarioData2Reduce.mat')
% load('/home/djoroya/Documentos/GitHub/DyConData/code/DarioData/data/tspan2.mat')
solution;
%solution = solution(1:1500,:,:);

%% See data


%%

[Nexp, Nt, Nvar] = size(solution);

%
sym_solution = size(solution);
for i = 1:Nexp
    soliy = reshape(solution(i,:,:),Nt,Nvar);
    soliy = fliplr(soliy);
    solution(Nexp+i,:,:) = reshape(soliy,1,Nt,Nvar);
end
%%
[Nexp, Nt, Nvar] = size(solution);
%% 
nn =2;

tspanfn = linspace(tspan(1),tspan(end),nn*Nt);

[xms   , tms  ] = meshgrid(1:Nvar,tspan);
[xmsfn , tmsfn] = meshgrid(1:Nvar,tspanfn);

%%
PreAllData.Y = zeros(Nexp*Nt,Nvar);
PreAllData.dYdt =  zeros(Nexp*Nt,Nvar);
PreAllData.Yfn =  zeros(Nexp*Nt*nn,Nvar);
PreAllData.dYdtfn = zeros(Nexp*Nt*nn,Nvar);


for i = 1:Nexp
    Data(i).Y     = reshape(solution(i,:,:),Nt,Nvar);
    Data(i).Yfine = interp2(xms,tms,Data(i).Y,xmsfn,tmsfn,'linear');
    %Data(i).Yfine = Data(i).Y;
    %
    dYdt = zeros(Nvar,Nt);
    for j = 1:Nvar
        dYdt(j,:) = gradient(Data(i).Y(:,j),tspan);
    end
    dYdt = dYdt';

    Data(i).dYdt  = dYdt;
    %
    dYdtfine = zeros(Nvar,nn*Nt);
    for j = 1:Nvar
        dYdtfine(j,:) = gradient(Data(i).Yfine(:,j),tspanfn);
    end
    dYdtfine = dYdtfine';
    Data(i).dYdtfine  = dYdtfine;
    %Data(i).dYdtfine  = interp2(xms,tms,dYdt,xmsfn,tmsfn,'spline');

    %
    PreAllData.Y((i-1)*Nt + 1:i*Nt,:)            = Data(i).Y;
    PreAllData.dYdt((i-1)*Nt + 1:i*Nt,:)         = Data(i).dYdt;
    PreAllData.Yfn((i-1)*Nt*nn + 1:i*Nt*nn,:)    = Data(i).Yfine;
    PreAllData.dYdtfn((i-1)*Nt*nn + 1:i*Nt*nn,:) = Data(i).dYdtfine;

    %
end

%% 


%% Remove Zero Data

AllData.Y    = zeros(size(PreAllData.Y));
AllData.dYdt = zeros(size(PreAllData.dYdt));

eps = 1e-5;
iter = 0;
for i = 1:Nt*Nexp*nn
    if norm(PreAllData.Yfn(i,:)) > eps
        iter = iter + 1;
        AllData.Y(iter,:)    =  PreAllData.Yfn(i,:)   ;
        AllData.dYdt(iter,:) =  PreAllData.dYdtfn(i,:);
    end
end

AllData.Y = AllData.Y(1:iter,:);
AllData.dYdt = AllData.dYdt(1:iter,:);

%%
% We can multiply the number of data by interpolation 

%%
%%
dYdt = AllData.dYdt';
Y    = AllData.Y';

% funcbasis = {@(Y) Y };
%%
% funcbasis = {    @(Y) Y       , ...
%                  @(Y) Y.^3    ,...
%                  @(Y) sin(Y)};
             
funcbasis = {    @(Y) Y       , ...
 @(Y) (Y).^3};
% %%
% funcbasis = {@(Y) Y,@(Y) sin(Y),@(Y) cos(Y) } ;
% % %              
 Theta = [];
for ifcn = funcbasis
    Theta = [Theta ;ifcn{:}(Y)];
end

%%
%A = (Theta')\(dYdt');
A = lsqminnorm(Theta',dYdt');

A = A';
%
figure
surf(A)


view(0,-90)
% daspect([1 1 1e-4])
colorbar



%%
number_of_exp = 3;
Y0 = Data(number_of_exp).Y(1,:);     

newtspan = tspan;

Ntsmall = 200;
newtspan = tspan(1:Ntsmall);
%
[~ ,Yestimation] = ode23(@(t,Y) estifnc(t,Y,funcbasis,A,Nvar), newtspan,Y0);
%
figure
subplot(1,3,1)
surf(Yestimation)
title('Estimation')
zlim([-2 2])
subplot(1,3,2)
surf(Data(number_of_exp).Y(1:Ntsmall,:))
title('Real')
zlim([-2 2])

subplot(1,3,3)
surf(Data(number_of_exp).Y(1:Ntsmall,:) - Yestimation)
title('Diff')
zlim([-2 2])

%%
AOpt = lsqminnorm(Theta',dYdt')
surf(AOpt-A)
%%
A0 = zeros(Nvar,2*Nvar);
A0 = A;
options = optimoptions(@fmincon,'display','iter','MaxFunctionEvaluations',1e6,'SpecifyObjectiveGradient',false);
AOpt = fmincon(@(A)Abfunction(A,dYdt,Theta),A0, [],[], ...
                                                [],[], ...
                                                [],[], ...
                                                @nlc,options);
%%
A0 = [0 0 0];

options = optimoptions(@fminunc,'display','iter','MaxFunctionEvaluations',1e4,'SpecifyObjectiveGradient',false);
AOpt = fminunc(@(A)Abfunctionfew(A,dYdt,Theta),A0,options);
%%
%%

function result = estifnc(t,Y,funcbasis,A,Nvar)
result = zeros(Nvar,1);
iter = 1;
for ifcn = funcbasis
    result = result + A(:,1+(iter-1)*Nvar:iter*Nvar)*ifcn{:}(Y);
    iter = iter +1;
end
end


function [res ,dres] = Abfunction(A,dYdt,Theta)

    
    L2  = norm(A*Theta - dYdt,2);
    %AL1 = norm(A,1);
    res = L2; %+ AL1; 
    if nargout > 1
        dres = Theta;
    end
end

function [res ,dres] = Abfunctionfew(params,dYdt,Theta)
    n=50;
    e = params(1)*ones(n,1);
    d = params(2)*ones(n,1);
    A = spdiags([e d e],-1:1,n,n);
    %
    A(:,1)   = 0;
    A(1,:)   = 0;
    A(:,end) = 0;
    A(end,:) = 0;
    
    A3 = params(3)*eye(n);
    A3(:,1)   = 0;
    A3(1,:)   = 0;
    A3(:,end) = 0;
    A3(end,:) = 0;    
    
    A = [A;A3]';
    
    L2  = norm(A*Theta - dYdt,2);
    %AL1 = norm(A,1);
    res = L2; %+ AL1; 
    if nargout > 1
        dres = Theta;
    end
end

function [c ,ceq] = nlc(A)
veclineal = A(:,1:50)  - A(:,1:50)';
vecsin    = A(:,51:100)- A(:,51:100)';

ceq =[veclineal(:) vecsin(:)];
c = [];
end

