%%
%%
clear all
% 
% load('/home/djoroya/Documentos/Data/DarioData.mat')
% load('/home/djoroya/Documentos/Data/tspan.mat')
load('/home/djoroya/Documentos/GitHub/Data/DarioData2Reduce.mat')
load('/home/djoroya/Documentos/GitHub/Data/tspan2.mat')
allsolution = solution;
[Nexp, Nt, Nvar] = size(solution);
indexspace = [2:(0.5*Nvar-1) (0.5*Nvar+2):Nvar-1];
%indexspace = 1:Nvar;
solution = solution(1000:2000,:,indexspace);
%
[Nexp, Nt, Nvar] = size(solution);
% 

%
AllData.Y = zeros(Nexp*Nt,Nvar);
AllData.dYdt =  zeros(Nexp*Nt,Nvar);
AllData.Yfn =  zeros(Nexp*Nt,Nvar);
AllData.dYdtfn = zeros(Nexp*Nt,Nvar);

for i = 1:Nexp
    Yi     = reshape(solution(i,:,:),Nt,Nvar);
    %
    dYdt = zeros(Nvar,Nt);
    for j = 1:Nvar
        dYdt(j,:) = gradient(Yi(:,j),tspan);
    end

    %Data(i).dYdt  = dYdt;
    %
    dYdtfine = zeros(Nvar,nn*Nt);
    %
    AllData.Y((i-1)*Nt + 1:i*Nt,:)            = Yi;
    AllData.dYdt((i-1)*Nt + 1:i*Nt,:)         = dYdt';
    %
end
%
% We can multiply the number of data by interpolation 
%
dYdt = AllData.dYdt';
Y    = AllData.Y';
save('/home/djoroya/Documentos/GitHub/Data/dydt.mat','dYdt','Y','tspan','solution')
%%
clear all
load('/home/djoroya/Documentos/GitHub/Data/dydt.mat')
[Nexp, Nt, Nvar] = size(solution);
%%
 %funcbasis = {@(Y) Y };
%
funcbasis = {@(Y) Y   , ...
           @(Y) (Y).^3};
% %%
% funcbasis = {@(Y) Y, ...
%                  @(Y) sin(Y)};
% %%
% funcbasis = {@(Y) Y     , ...
%              @(Y) Y.^3  , ...
%              @(Y) Y.^5} ;
% %              
 Theta = [];
for ifcn = funcbasis
    Theta = [Theta ;ifcn{:}(Y)];
end
Theta(abs(Theta)<1e-3) = 0;
Theta = sparse(Theta);
%%
A = (Theta')\(dYdt');
%A = pinv(Theta')*(dYdt');
%A = lsqminnorm(Theta',dYdt');
A = A';
%%
figure
spy(abs(A)>10)

%%
surf(A)
view(0,-90)
%%
number_of_exp = 5;
Y0 = reshape(solution(number_of_exp,1,:),1,Nvar);  
Yt = reshape(solution(number_of_exp,:,:),Nt,Nvar);

%
[~ ,Yestimation] = ode23tb(@(t,Y) estifnc(t,Y,funcbasis,A,Nvar), tspan,Y0);
%%
figure
subplot(1,3,1)
surf(Yestimation)
title('Estimation','FontSize',15)
shading interp

zlim([-3 3])
subplot(1,3,2)
surf(Yt)
title('Real','FontSize',15)
zlim([-3 3])
shading interp
subplot(1,3,3)
surf(Yt - Yestimation)
title('Difference','FontSize',15)
shading interp

zlim([-3 3])

%%
AOpt = lsqminnorm(sparse(Theta'),sparse(dYdt'),1e-1);
surf(AOpt'-A)
%%
A0 = zeros(Nvar,2*Nvar);
A0 = A';
options = optimoptions(@fmincon,'display','iter','MaxFunctionEvaluations',1e6,'SpecifyObjectiveGradient',false);
AOpt = fmincon(@(A)Abfunction(A,dYdt,Theta),A0, [],[], ...
                                                Theta',dYdt', ...
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

    
    %L2  = norm(A*Theta - dYdt,2);
    res = norm(A,1);
    %res = L2; %+ AL1; 
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

