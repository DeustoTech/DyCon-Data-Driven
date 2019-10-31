%%
%%
clear all
% 
load('/home/djoroya/Documentos/GitHub/Data/DarioData01.mat')
load('/home/djoroya/Documentos/GitHub/Data/tspan01.mat')
tspan = time;
solution = solution(:,1:50,2:end-1);
tspan = tspan(1:50);
% solution1 = repmat(solution(1,:,:),10,1,1);
% solution1 = solution1 +  0.1*rand(size(solution1));
%  solution = [solution1 ;solution];
%solution = solution(1:1500,:,:);

%% See data
%%
% [Nexp, Nt, Nvar] = size(solution);
% %
% sym_solution = size(solution);
% for i = 1:Nexp
%     soliy = reshape(solution(i,:,:),Nt,Nvar);
%     soliy = fliplr(soliy);
%     solution(Nexp+i,:,:) = reshape(soliy,1,Nt,Nvar);
% end
%%
[Nexp, Nt, Nvar] = size(solution);
%% 
nn =1;

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
    Data(i).Yfine = interp2(xms,tms,Data(i).Y,xmsfn,tmsfn,'spline');
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
%
xline = linspace(-1,1,48);
dx = xline(2) - xline(1);
ALaplacian = FDLaplacian(xline);
%
Dxx = @(Y) reshape(cell2mat(arrayfun(@(i)gradient(xline,gradient(xline,Y(:,i))),1:length(Y(1,:)),'UniformOutput',0)),Nvar,length(Y(1,:)));
Dx  = @(Y) reshape(cell2mat(arrayfun(@(i) gradient(xline,Y(:,i)),1:length(Y(1,:)),'UniformOutput',0)),Nvar,length(Y(1,:)));
%
%%
Yxx = ArrayGradient(xline,ArrayGradient(xline,Y));

%funcbasis = {@(Y) Grad(Y)};

%%
funcbasis = {};

ind = 0;
for iter = 1:2:7
    ind = ind + 1;
%     funcbasis{iter} = @(Y) Y.^(2*iter-1);
    funcbasis{ind} = @(Y) chebfunbasis(Y,iter);

end
%%
%funcbasis = {    @(Y) Y };
%
% funcbasis = {    @(Y) Y       , ...
%                  @(Y) Y.^3    ,...
%                  @(Y) sin(Y)};
%
%%
%
% funcbasis = {  @(Y) dx^2*ALaplacian*Y, @(Y) sin(Y),@(Y) Y.^3 };
%%
%funcbasis = {@(Y) ArrayGradient(xline,ArrayGradient(xline,Y)),@(Y) sin(Y) } ;
funcbasis = {@(Y) -ALaplacian*Y ,@(Y) sin(pi*Y) } ;
%%              
 Theta = [];
for ifcn = funcbasis
    Theta = [Theta ;ifcn{:}(Y)];
end

%%
%A = (Theta')\(dYdt');
A = lsqminnorm(Theta',dYdt');
A = A';
%%
opti = casadi.Opti();  % CasADi function
%
Nx = 48;
Acas = opti.variable(Nx,2*Nx); % state trajectory
opti.minimize( sum(sum((Acas*Theta - dYdt).^2) + sum(sum(abs(Acas))))); % minimizing time

% opti.subject_to(Acas(1:Nx,Nx+1:2*Nx)==Acas(1:Nx,Nx+1:2*Nx)');
%opti.subject_to(Acas(1:Nx,1:Nx)==Acas(1:Nx,1:Nx)');

opti.solver('ipopt'); % set numerical backend
sol = opti.solve();   % actual solve
A = sol.value(Acas);
%%

%A = A';
%
figure
surf(A)


view(0,-90)
% daspect([1 1 1e-4])
colorbar



%%
number_of_exp =30;
Y0 = Data(number_of_exp).Y(1,:);     
%
Ntmax = 20;
close all
figure('unit','norm','pos',[0 0 1 1])
L = 1.5859;
L = 1.0 % <-- best
xline = linspace(0,L,48);
ALap = FDLaplacian(xline);
Ysym = sym('y',[48 1]);
Usym = sym('u',[1 1]);
F = @(t,Y,U,Params) ALap*Y + sin(pi*Y);
idynamics = pde(F,Ysym,Usym);
idynamics.InitialCondition = Y0';
idynamics.Solver = @ode23;
idynamics.Nt = 200;
idynamics.FinalTime = 0.1;
[~,YDyCon] = solve(idynamics);
subplot(2,2,1)
surf(YDyCon(1:Ntmax,:))
zlim([-2 2])
title('DyCon')
%
[~ ,Yestimation] = ode45(@(t,Y) estifnc(t,Y,funcbasis,A,Nvar), tspan,Y0);
%
subplot(2,2,2)
surf(Yestimation(1:Ntmax,:))
title('Estimation')
zlim([-2 2])
subplot(2,2,3)
surf(Data(number_of_exp).Y(1:Ntmax,:))
title('Real')
zlim([-2 2])

subplot(2,2,4)
surf(Data(number_of_exp).Y(1:Ntmax,:) - Yestimation(1:Ntmax,:))
title('Difference')
%zlim([-2 2])

%%

function result = estifnc(t,Y,funcbasis,A,Nvar)
result = zeros(Nvar,1);
iter = 1;
for ifcn = funcbasis
    result = result + A(:,1+(iter-1)*Nvar:iter*Nvar)*ifcn{:}(Y);
    iter = iter +1;
end
end
