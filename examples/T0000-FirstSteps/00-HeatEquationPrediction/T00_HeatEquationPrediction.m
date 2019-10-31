clear all;
Nx = 1000;
xline = linspace(-1,1,Nx);
%A = -FEFractionalLaplacian(0.5,1,Nx);
A = FDLaplacian(xline);

%% 
% Number of degree of fredom of matrix A is: N^2
% Number of degree of fredom of matrix B is: N^2
% Total Numer of unknowns quatities; 2N^2
%% 
% For this reason, we need almost 2N^2 mesaurements
Nt = 50*Nx^2;
dt = 0.1/Nt;
ADiscrete = dt*A + eye(Nx);
%%
tspan = linspace(0,5,Nt);
[xms, tms] = meshgrid(xline,tspan);
%%

%%
%x(:,1) = cos(0.5*pi*xline');
x(:,1) = rand(size(xline'));
%x(:,1) = xline;

for k = 1:Nt
    x(:,k+1) = ADiscrete*x(:,k); 
end
%%

X  = x(:,1:end-1);
X2 = x(:,2:end);

ADisEsti   = (X2)*pinv(X);
AConEsti = (ADisEsti  - eye(Nx))/dt;
%%
figure
subplot(1,3,1)
surf(AConEsti)%%
subplot(1,3,2)
surf(A)%%
subplot(1,3,3)
surf(A-AConEsti)%%
