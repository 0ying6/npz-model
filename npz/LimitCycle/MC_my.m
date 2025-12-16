% Monte Carlo simulation for a 3D dynamical system
clc; clear;
tic

epsilon=0.2;
a=0.2;
    b=0.2;
    c=0.4;
    d=1.396;
    e=0.03;
    k=0.05;
    r=0.15;
    s=0.04;
    N_0=0.6;
    alpha=0.25;
    beta=0.33;
    gamma=0.5;
    lamda=0.6;
    mu=0.035;
    
% Parameters
num_particles = 1000; % Number of Monte Carlo particles
N = 5;
dt = 0.0025; % Time step
%radius = 1; % Radius of the sphere
%D = 0.9;% Diffusion coefficient
% T = 100;
% t = 0 : dt: T;
% Drift term 
f = @(x,y,z)    -(x/(e+x))*(a/(b+c*y))*y+r*y+(beta*lamda*(y^2)/(mu^2+y^2))*z+gamma*d*(z^2)+k*(N_0-x); 
g = @(x,y,z)    (x/(e+x))*(a/(b+c*y))*y-r*y-(lamda*(y^2)/(mu^2+y^2))*z-(s+k)*y;
h = @(x,y,z)    (alpha*lamda*(y^2)/(mu^2+y^2))*z-d*(z^2);

% Initialize array to store escape times
escape_times = zeros(num_particles, 1);
LCV = load('LimitCycle_1.txt');
%  x(1) =0.8;
%  y(1) =0.3;
%  z(1) =0.31;
% Monte Carlo loop
for i = 1:num_particles
    % Initialize particle position at the origin
    u_L=0.2258;
    v_L=0.0469;
    w_L=0.0690;
    position=[u_L,v_L,w_L];
    time = 0;
    
    % Simulate until the particle exits the sphere
    while norm(position)<norm(LCV(:,i))
        % Update time
        time = time + dt;
        
        % Compute deterministic drift
    
    deterministic_drift1 = f(u_L,v_L,w_L) * dt;
    deterministic_drift2 = g(u_L,v_L,w_L) * dt;
    deterministic_drift3 = h(u_L,v_L,w_L) * dt;
    deterministic_drift = [deterministic_drift1,deterministic_drift2,deterministic_drift3];
%     x(i) = u_L + f(u_L,v_L,w_L)*dt+W1(num_particles)*epsilon*x(i);
%     y(i) = v_L + g(u_L,v_L,w_L)*dt+W2(num_particles)*epsilon*y(i);
%     z(i) = w_L + h(u_L,v_L,w_L)*dt+W3(num_particles)*epsilon*z(i);

        % Generate random displacement due to Brownian motion
        W1 = sqrt(2*0.8*epsilon*dt)*randn(1,3); 
        W2 = sqrt(2*0.3*epsilon*dt)*randn(1,3); 
        W3 = sqrt(2*0.31*epsilon*dt)*randn(1,3); 
        if abs(position-LCV(:,i))<=0.001
            break;
        end
        
        % Update position
        position = position + deterministic_drift + u_L * W1;
%         position(1) = position(1) + deterministic_drift1 + u_L * W1(1);
%         position(2) = position(2) + deterministic_drift2 + v_L * W1(2);
%         position(3) = position(3) + deterministic_drift3 + w_L * W1(3);
    end
    
    % Record escape time
    escape_times(i) = time;
end

% Compute mean escape time
mean_escape_time = mean(escape_times);

% Display results
fprintf('Mean escape time: %.4f\n', mean_escape_time);

% Optional: plot histogram of escape times
figure;
histogram(escape_times, 30, 'Normalization', 'pdf');
xlabel('Transition Time');
ylabel('Probability Distribution');
title('Distribution of transition times with \epsilon=0.8');
toc
