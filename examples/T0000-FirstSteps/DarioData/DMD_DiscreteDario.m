%%
%%
clear all
% 
% load('/home/djoroya/Documentos/GitHub/DyConData/code/DarioData/data/DarioData.mat')
% load('/home/djoroya/Documentos/GitHub/DyConData/code/DarioData/data/tspan.mat')
load('/home/djoroya/Documentos/GitHub/DyConData/code/DarioData/data/DarioData2Reduce.mat')
load('/home/djoroya/Documentos/GitHub/DyConData/code/DarioData/data/tspan2.mat')
allsolution = solution;
solution = solution(1:200,:,:);
%%
[Nexp, Nt, Nvar] = size(solution);

X  = reshape(solution(10,1:end-1,:),Nt-1,Nvar );
X2 = reshape(solution(10,2:end,:),Nt-1,Nvar );
[U,S,V] = svd(X,'econ');

%%  Compute DMD (Phi are eigenvectors)
r = 10;  % truncate at 21 modes
U = U(:,1:r);
S = S(1:r,1:r);
V = V(:,1:r);
Atilde = U'*X2*V*inv(S);
[W,eigs] = eig(Atilde);
Phi = X2*V*inv(S)*W;

%%  Plot DMD spectrum
figure
theta = (0:1:100)*2*pi/100;
plot(cos(theta),sin(theta),'k--') % plot unit circle
hold on, grid on
scatter(real(diag(eigs)),imag(diag(eigs)),'ok')
axis([-1.1 1.1 -1.1 1.1]);