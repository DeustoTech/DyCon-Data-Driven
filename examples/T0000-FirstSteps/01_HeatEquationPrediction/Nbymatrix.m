function [AConEsti,xms,yms] = Nbymatrix(Nx)

xline = linspace(-1,1,Nx);
A = -FEFractionalLaplacian(0.75,1,Nx);
%A = FDLaplacian(xline);
a = -0.5;
b = +0.5;
B = BInterior(xline,a,b);
%% 
% Number of degree of fredom of matrix A is: N^2
% Total Numer of unknowns quatities; N^2
%% 
% For this reason, we need almost 2N^2 mesaurements
Nt = 2*Nx^2;
dt = 0.5/Nt;
ADiscrete = dt*A + eye(Nx);
%%
tspan = linspace(0,5,Nt);
[xms, tms] = meshgrid(xline,tspan);
%%
u = rand(Nx,Nt);
%%
x(:,1) = sin(pi*xline');
for k = 1:Nt
    %u(:,k) = -K*x(:,k);
    x(:,k+1) = ADiscrete*x(:,k) + (dt*B)*u(:,k);
end
%%
Upsilon = u;
X  = x(:,1:end-1);
X2 = x(:,2:end);

ADisEsti   = (X2 - (dt*B)*Upsilon)*pinv(X);
AConEsti = (ADisEsti  - eye(Nx))/dt;

%% mesh

[xms, yms] = meshgrid(xline,xline);

end

