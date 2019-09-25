%%
clear;close all
load('/home/djoroya/Documentos/GitHub/DyConData/code/DarioData/data/DarioData.mat')
load('/home/djoroya/Documentos/GitHub/DyConData/code/DarioData/data/tspan.mat')

%%
xline = linspace(-1,1,48);
A = FDLaplacian(xline);
idyn = ode('A',A);
idyn.InitialCondition = cos(0.5*pi*xline');
Nt = 200;
Nvar = 48
idyn.Nt = Nt;
[~ ,exp1] = solve(idyn);
%%
dt = tspan(2) - tspan(1);

Ntmax = 200;
%solution = solution(:,1:Ntmax,2:end-1);
tspan = tspan(1:Ntmax);
%%
%[Nexp, Nt ,Nvar] = size(solution);
% figure('unit','norm','pos',[0 0.5 0.5 0.5])
% Nsqrt = 3;
% for iter = 1:Nsqrt^2
%     subplot(Nsqrt,Nsqrt,iter)
%     surf(reshape(solution(iter,:,:),Nt,Nvar))
%     title("Experiment "+iter)
% end
%%
%exp1 = reshape(solution(1,:,:),Nt,Nvar);
%%
dXdt = zeros(Nvar,Nt);
for iter = 1:Nvar
    dXdt(iter,:) = gradient(exp1(:,iter),tspan);
end
dXdt = dXdt';
%%
X = exp1;
order = 1;
Theta = [];
for iter = 1:2:order
    Theta = [Theta X.^iter];
end
Theta = X;
%%
