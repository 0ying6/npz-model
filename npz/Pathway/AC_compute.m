clear,clc

tic
%% initail parameters
epsilon=0.9;
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
%% Lagrangian
f = @(x,y,z)    -(x/(e+x))*(a/(b+c*y))*y+r*y+(beta*lamda*(y^2)/(mu^2+y^2))*z+gamma*d*(z^2)+k*(N_0-x);
g = @(x,y,z)    (x/(e+x))*(a/(b+c*y))*y-r*y-(lamda*(y^2)/(mu^2+y^2))*z-(s+k)*y;
h = @(x,y,z)    (alpha*lamda*(y^2)/(mu^2+y^2))*z-d*(z^2);
 
fx = @(x,y,z) (-e/((e+x)^2))*(a/(b+c*y))*y-k;
fxx = @(x,y,z) 2*e*(a*y/(b+c*y))/((e+x)^3);
fxy = @(x,y,z) (-e/((e+x)^2))*a*b/((b+c*y)^2);
fy = @(x,y,z) -(x/(e+x))*a*b/((b+c*y)^2)+r+2*beta*lamda*y*z*(mu^2)/((mu^2+y^2)^2);
fz = @(x,y,z)  beta*lamda*(y^2)/(mu^2+y^2)+gamma*2*d*z;
fxz = @(x,y,z) 0;

gx = @(x,y,z) (a*y/(b+c*y))*e/((e+x)^2);
gy = @(x,y,z) x/(e+x)*a*b/((b+c*y)^2)-r-2*lamda*y*z*mu^2/((mu^2+y^2)^2)-s-k;
gz = @(x,y,z)  -lamda*y^2/(mu^2+y^2);
gyx = @(x,y,z)  (a*b/(b+c*y)^2)*e/((e+x)^2);
gyy = @(x,y,z)  -2*a*b*c*x/(e+x)/((b+c*y)^3)-2*lamda*z*(mu^4-3*mu^2*y^2)/((mu^2+y^2)^3);
gyz = @(x,y,z)   -2*lamda*mu^2*y/((mu^2+y^2)^2);

hx = @(x,y,z)  0;
hy = @(x,y,z)  2*lamda*alpha*z*y*mu^2/((mu^2+y^2)^2);
hz = @(x,y,z)  alpha*lamda*(y^2)/(mu^2+y^2)-2*d*z;
hzx = @(x,y,z)  0;
hzy = @(x,y,z)  2*lamda*alpha*y*mu^2/((mu^2+y^2)^2);
hzz = @(x,y,z)  -2*d;

L = @(x,y,z,u,v,w) (1/2)*(epsilon*x)^(-2)*(u-f(x,y,z)-(1/2)*epsilon^2*x)^2+(1/2)*(epsilon*y)^(-2)*(v-g(x,y,z)-(1/2)*epsilon^2*y)^2+(1/2)*(epsilon*z)^(-2)*(w-h(x,y,z)-(1/2)*epsilon^2*z)^2+(1/2)*(fx(x,y,z)+gy(x,y,z)+hz(x,y,z)-f(x,y,z)/x-g(x,y,z)/y-h(x,y,z)/z);
%% Compute the action functional
T = 1;
uu = load('subdatau_opt.mat');
vv = load('subdatav_opt.mat');
ww = load('subdataw_opt.mat');
u = uu.x1n;
v = vv.x3n;
w = ww.x5n;
% 
[m n] = size(u);
for i = 1 : m
    Act(i) = ActionValue(u(i,:),v(i,:),w(i,:),L,T);
end

%% Find the minimizer of the action functional
ind = find(Act==min(min(Act)));
%ind = find(Act==max(max(Act)));

u_L =0.2258;
v_L = 0.0469;
w_L = 0.0690;
for j = 1 : n
    u_opt(j) = u(ind,j);
    v_opt(j) = v(ind,j);
    w_opt(j) = w(ind,j);
end

L = load('LimitCycle_1.txt');


plot3(u_opt,v_opt,w_opt,'b'); hold on
plot3(L(1,:),L(2,:),L(3,:),'r'); hold on 
plot3(u_L,v_L,w_L,'*'); hold on
plot3(L(1,553),L(2,553),L(3,553),'o');hold on
plot3(u_opt(1001),v_opt(1001),w_opt(1001),'*');hold on
grid on
% plot(L(1,j+1),L(2,j+1),'o');

toc