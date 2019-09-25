%%
clear
load('/home/djoroya/Documentos/GitHub/DyConData/code/data/DarioData.mat')
load('/home/djoroya/Documentos/GitHub/DyConData/code/data/tspan.mat')

Ntmax = 200;
dt = tspan(2) - tspan(1);
tspan = tspan(1:Ntmax);
solution = solution(:,1:Ntmax,1:end);

%%
[Nexp, Nt ,Nvar] = size(solution);
clf
Nsqrt = 4;
for iter = 1:Nsqrt^2
    subplot(Nsqrt,Nsqrt,iter)
    surf(reshape(solution(iter,:,:),Nt,Nvar))
    title("Experiment "+iter)

end
%%
exp1 = reshape(solution(1,:,:),Nt,Nvar);
%%
[xms, tms] = meshgrid(1:Nvar,tspan);

tspan_fine = linspace(tspan(1),tspan(end),4*Nt);
[xms_fine, tms_fine] = meshgrid(1:1:Nvar,tspan_fine);
%%
exp1_fine = interp2(xms,tms,exp1,xms_fine,tms_fine);

%%
dxdt = zeros(Nvar,Nt);
for iter = 1:Nvar
    dxdt(iter,:) = gradient(exp1(:,iter),tspan);
end
dxdt = dxdt';
%%
X = exp1;
surf(dxdtX)
