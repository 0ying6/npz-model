function [o] = ActionValue(u,v,w,L,T)
%% Compute the action functional 
%% L: Lagrangian; u,v: paths; T: Time
%T = 4;
n = length(u);
dt = T/n;

o = 0;
for i = 2 : n-1
    o = o + L(u(i),v(i),w(i),(u(i+1)-u(i-1))/(2*dt),(v(i+1)-v(i-1))/(2*dt),(w(i+1)-w(i-1))/(2*dt))*dt/3;
end