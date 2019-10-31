clear all;
Nx = 50;
xline = linspace(-1,1,Nx);
%A = -FEFractionalLaplacian(0.5,1,Nx);
A = FDLaplacian(xline);
a = -0.5;
b = +0.5;
B = BInterior(xline,a,b);
%% 
% Number of degree of fredom of matrix A is: N^2
% Number of degree of fredom of matrix B is: N^2
% Total Numer of unknowns quatities; 2N^2
%% 
% For this reason, we need almost 2N^2 mesaurements
Nt = 2*Nx^2;
dt = 5/Nt;
ADiscrete = dt*A + eye(Nx);
K = eye(Nx);
%%
tspan = linspace(0,5,Nt);
[xms, tms] = meshgrid(xline,tspan);
%%
u = tms.*sin(3*pi*xms).*sin(3*pi*tms);

u = u';
%u = 0.001*u'
%u = ones(Nx,Nt);
u = 100*rand(Nx,Nt);
%%
x(:,1) = sin(pi*xline');
for k = 1:Nt
    %u(:,k) = -K*x(:,k);
    x(:,k+1) = ADiscrete*x(:,k) + (dt*B)*u(:,k); 
end
%%
tline = linspace(0,10,Nt)
dXdt = ArrayGradient(tline,x);

%% Unknow B matrix, we create 

XUps = [x ;Upsilon];
%%
%ABmatrix = X2*pinv(XUps);
ABmatrix = XUps'\X2';
%ABmatrix = lsqminnorm(XUps',X2');
ABmatrix = ABmatrix';

%%
ADisEsti = ABmatrix(:,1:Nx);
BDisEsti = ABmatrix(:,Nx+1:end);

AConEsti = (ADisEsti  - eye(Nx))/dt;
BConEsti = BDisEsti/dt;

%%
figure
subplot(2,3,1)
surf(AConEsti)%%
title('A Estimation')

subplot(2,3,2)
surf(A)%%
title('A matrix')

subplot(2,3,3)
surf(A-AConEsti)%%
title('\Delta A')

subplot(2,3,4)
surf(BConEsti)%%
title('B Estimation')

subplot(2,3,5)
surf(B)%%
title('B matrix')

subplot(2,3,6)
surf(B-BConEsti)%%
title('\Delta B')