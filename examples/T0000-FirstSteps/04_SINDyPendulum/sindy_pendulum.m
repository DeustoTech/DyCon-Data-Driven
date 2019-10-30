clear all, close all, clc

%% parameters for simulation

N = 1; % number of simulations
dt = 0.1;
T = 20;
% d = pi/100; % initial perturbation parameter
% x0 = [pi+d; 1]; % initial condition around unstable equilibria


%% data simulation: damped nonlinear pendulum with random sampling

f = @(t,x)[x(2) ; -sin(x(1))];

new_sol = [];
new_difX = [];

for i = 1:N
    difX = [];
    init_theta = 2*pi*rand;
    init_vel = 10*rand;
    x0 = [init_theta; init_vel];
    [t, sol] = ode45(f, 0:dt:T, x0);
    for j = 1 : 2
        for i = 1 : length(sol(:,j))-1
            difX(i,j) = (sol(i+1,j)-sol(i,j))/dt;
        end
    end
    
    for j = 1 : 2
        difX(length(sol(:,j)),j) = (sol(length(sol(:,j)),j)-sol(length(sol(:,j))-1,j))/dt;
    end
  
    new_sol = [new_sol ; sol];
    new_difX = [new_difX ; difX];
end

%% defining matrix of observables

X = [new_sol sin(new_sol)]; % we just take sin()

%% lasso regression

reg_1 = lasso(X, new_difX(:,1));
reg_2 = lasso(X, new_difX(:,2));