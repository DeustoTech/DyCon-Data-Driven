clear all;
Nx = 40;
xline = linspace(-1,1,Nx);
A = -FEFractionalLaplacian(0.8,1,Nx);
M =  MassMatrix(xline);
A = M\A;
%A = FDLaplacian(xline);
a = -0.5;
b = +0.5;
B = BInterior(xline,a,b);
B = M\B;

%% 
% Number of degree of fredom of matrix A is: N^2
% Number of degree of fredom of matrix B is: N^2
% Total Numer of unknowns quatities; 2N^2
%% 
% For this reason, we need almost 2N^2 mesaurements
Nt = 2*Nx^2;
dt = 0.5/Nt;
ADiscrete = dt*A + eye(Nx);
K = eye(Nx);
%%
tspan = linspace(0,5,Nt);
[xms, tms] = meshgrid(xline,tspan);
%% Zero Control
u = zeros(Nx,Nt);
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

%% Singular Value Descomposition

[U,S,V] = svd(X,'econ');

%% truncate at r modes
r = 5;  

U = U(:,1:r);
S = S(1:r,1:r);
V = V(:,1:r);

% Obtain new matrix
Atilde = U'*X2*V*inv(S);
[W,eigs] = eig(Atilde);
Phi = X2*V*inv(S)*W;

%%