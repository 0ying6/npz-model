%% carbon cycle sysytem
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
%% functions 
f = @(x,y,z)    -(x/(e+x))*(a/(b+c*y))*y+r*y+beta*lamda*(y^2)*z/(mu^2+y^2)+gamma*d*(z^2)+k*(N_0-x);
g = @(x,y,z)    (x/(e+x))*(a/(b+c*y))*y-r*y-lamda*(y^2)*z/(mu^2+y^2)-(s+k)*y;
h = @(x,y,z)    alpha*lamda*(y^2)*z/(mu^2+y^2)-d*(z^2);
 
fx = @(x,y,z) (-e/((e+x)^2))*a*y/(b+c*y)-k;
fxx = @(x,y,z) 2*e*(a*y/(b+c*y))/((e+x)^3);
fxy = @(x,y,z) (-e/((e+x)^2))*a*b/((b+c*y)^2);
fy = @(x,y,z) (-x/(e+x))*a*b/((b+c*y)^2)+r+2*beta*lamda*y*z*(mu^2)/((mu^2+y^2)^2);
fz = @(x,y,z)  beta*lamda*(y^2)/(mu^2+y^2)+gamma*2*d*z;
fxz = @(x,y,z) 0;

gx = @(x,y,z) (a*y/(b+c*y))*e/((e+x)^2);
gy = @(x,y,z) (x/(e+x))*a*b/((b+c*y)^2)-r-2*lamda*y*z*mu^2/((mu^2+y^2)^2)-s-k;
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

%% stable point
u_L =0.2258;
v_L = 0.0469;
w_L = 0.0690;
%% 
dis=@(x,y,z,x1,y1,z1) ((x-x1)^2+(y-y1)^2+(z-z1)^2)^(1/2);
T = 0.0344;
N = 10^4;
dt = T/N;
t = 0 : dt : T;
F = @(x,y,z,u,v,w) x^(-1)*y^2+fy(x,z,v)*u+fz(x,z,v)*w-x^(-1)*f(x,z,v)^2+f(x,z,v)*fx(x,z,v)-x^2/z^2*(u-g(x,z,v))*gx(x,z,v)-x^2/v^2*(w-h(x,z,v))*hx(x,z,v)+((epsilon^2*x^2)/2)*(fxx(x,z,v)+gyx(x,z,v)+hzx(x,z,v));
G = @(x,y,z,u,v,w) z^(-1)*u^2+gx(x,z,v)*y+gz(x,z,v)*w-z^(-1)*g(x,z,v)^2+g(x,z,v)*gy(x,z,v)-z^2/x^2*(y-f(x,z,v))*fy(x,z,v)-z^2/v^2*(w-h(x,z,v))*hy(x,z,v)+((epsilon^2*z^2)/2)*(fxy(x,z,v)+gyy(x,z,v)+hzy(x,z,v));
H = @(x,y,z,u,v,w) v^(-1)*w^2+hx(x,z,v)*y+hy(x,z,v)*u-v^(-1)*h(x,z,v)^2+h(x,z,v)*hz(x,z,v)-v^2/x^2*(y-f(x,z,v))*fz(x,z,v)-v^2/z^2*(u-g(x,z,v))*gz(x,z,v)+((epsilon^2*v^2)/2)*(fxz(x,z,v)+gyz(x,z,v)+hzz(x,z,v));

fprintf("读取数据�?...");

AA = load ('LimitCycle_1.csv');
BB = load ('LimitCycle_1_ans.csv');

fprintf("已完成\n");

[NN1,NN2] = size(AA);
x1=zeros(NN2,N+1);
x2=zeros(1,N+1);
x3=zeros(NN2,N+1);
x4=zeros(1,N+1);
x5=zeros(NN2,N+1);
x6=zeros(1,N+1);
num=0;
lst=zeros(1,0);
fprintf("NN1 NN2=%d %d\n",NN1,NN2)
for I=1:NN2
    if (mod(I,100)==0)
        fprintf("%d 计算结束\n",I)
    end
    a1 = AA(1,I);
    a2 = AA(2,I);
    a3 = AA(3,I);
    v1 = BB(1,I);
    v2 = BB(2,I);
    v3 = BB(3,I);
    x1(I,1) = u_L;
    x2(1) = v1;
    x3(I,1) = v_L;
    x4(1) = v2;
    x5(I,1) = w_L;
    x6(1) = v3;
    for i = 1 : N
        x1(I,i+1) = x1(I,i) + x2(i)*dt;
        x2(i+1) = x2(i) + F(x1(I,i),x2(i),x3(I,i),x4(i),x5(I,i),x6(i))*dt;
        x3(I,i+1) = x3(I,i) + x4(i)*dt;
        x4(i+1) = x4(i) + G(x1(I,i),x2(i),x3(I,i),x4(i),x5(I,i),x6(i))*dt;
        x5(I,i+1) = x5(I,i) + x6(i)*dt;
        x6(i+1) = x6(i) + H(x1(I,i),x2(i),x3(I,i),x4(i),x5(I,i),x6(i))*dt;
    end
end

fprintf("保存数据�?...");

csvwrite('x1_trace.csv',x1)
csvwrite('x3_trace.csv',x3)
csvwrite('x5_trace.csv',x5)
fprintf("已完成\n");
