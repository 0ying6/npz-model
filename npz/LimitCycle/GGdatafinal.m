%% carbon cycle sysytem
clear,clc

tic
%% initail parameters
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
    beta=0.7;
    gamma=0.5;
    lamda=0.6;
    mu=0.035;
%% functions 
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


%% main code
% u_L = (b-1)^(-1/gama)*c_p;
% v_L = y0+mu*(theta+nu-theta*c_p^gama/((b-1)*c_x^gama+c_p^gama));
T = 1;

F = @(x,y,z,u,v,w) x^(-1)*y^2+fy(x,z,v)*u+fz(x,z,v)*w-x^(-1)*f(x,z,v)^2+f(x,z,v)*fx(x,z,v)-x^2/z^2*(u-g(x,z,v))*gx(x,z,v)-x^2/v^2*(w-h(x,z,v))*hx(x,z,v)+((epsilon^2*x^2)/2)*(fxx(x,z,v)+gyx(x,z,v)+hzx(x,z,v));
G = @(x,y,z,u,v,w) z^(-1)*u^2+gx(x,z,v)*y+gz(x,z,v)*w-z^(-1)*g(x,z,v)^2+g(x,z,v)*gy(x,z,v)-z^2/x^2*(y-f(x,z,v))*fy(x,z,v)-z^2/v^2*(w-h(x,z,v))*hy(x,z,v)+((epsilon^2*z^2)/2)*(fxy(x,z,v)+gyy(x,z,v)+hzy(x,z,v));
H = @(x,y,z,u,v,w) v^(-1)*w^2+hx(x,z,v)*y+hy(x,z,v)*u-v^(-1)*h(x,z,v)^2+h(x,z,v)*hz(x,z,v)-v^2/x^2*(y-f(x,z,v))*fz(x,z,v)-v^2/z^2*(u-g(x,z,v))*gz(x,z,v)+((epsilon^2*v^2)/2)*(fxz(x,z,v)+gyz(x,z,v)+hzz(x,z,v));
%% main code
u_L =0.2258;
v_L = 0.0469;
w_L = 0.0690;
% u_L =0.9614;
% v_L = 0.0068;
% w_L = 0.0247;
N = 50;
dud = 0 : 2*pi/N : 2*pi;
uu =  0 : 2*pi/N : 2*pi;
vec1 =sin(dud).*cos(uu);
vec2 =sin(dud).*sin(uu);
vec3 =cos(dud);
M = 5;
vel =1e-2;
NN = 10^3;
dt = T/NN;
t = 0 : dt : T;
dy = @(t, y)[y(2);
             F(y(1),y(2),y(3),y(4),y(5),y(6));
             y(4);
             G(y(1),y(2),y(3),y(4),y(5),y(6));
             y(6);
             H(y(1),y(2),y(3),y(4),y(5),y(6))];
tspan = [0,T];


p1=0;
P=zeros(7,M*(N+1)^2);
u = load('LimitCycle_1.txt');
% plot(u(1,:),u(2,:),'r'); hold on
for k = 1 : M
    fprintf("k=%d;\n",k);
    out1 = (vel*k/M)*vec1 ;
    out2 = (vel*k/M)*vec2;
    out3 = (vel*k/M)*vec3;
    for i = 1 : N+1
        for j = 1 : N+1
            for o = 1 : N+1
                y0 = [u_L out1(i) v_L out2(j) w_L out3(o)];
               %options = odeset('RelTol', 1e12, 'AbsTol', 1e-8);
               [x, u2] = ode45(dy, tspan, y0);
                %[x,u2] = ode45(dy,tspan,y0);
                p1=p1+1;
                P(1,p1) =out1(i);
                P(2,p1) = out2(j);
                P(3,p1) = out3(o);
                P(4,p1) = u2(end,1);
                P(5,p1) = u2(end,3);
                P(6,p1) = u2(end,5);
                P(7,p1) = x(end);
            end
        end
    end
end
% for k = M : M
%     out1 = rur(k)*vec1;
%     out2 = rur(k)*vec2;
% for i = 1 : N+1
%     for j = 1 : N+1            
%         x1(1) = u_L;
%         x2(1) = out1(i);
%         x3(1) = v_L;
%         x4(1) = out2(j);
%         for l = 1 : NN
%             x1(l+1) = x1(l) + x2(l)*dt;
%             x2(l+1) = x2(l) + F(x1(l),x2(l),x3(l),x4(l))*dt;
%             x3(l+1) = x3(l) + x4(l)*dt;
%             x4(l+1) = x4(l) + G(x1(l),x2(l),x3(l),x4(l))*dt;
%         end
%         P(1,(k-1)*(N+1)^2+(i-1)*(N+1)+j) = out1(i);
%         P(2,(k-1)*(N+1)^2+(i-1)*(N+1)+j) = out2(j);
%         P(3,(k-1)*(N+1)^2+(i-1)*(N+1)+j) = x1(end);
%         P(4,(k-1)*(N+1)^2+(i-1)*(N+1)+j) = x3(end);
%         plot(x1,x3,'-o');hold on
%     end
% end
% end

save shiyan9 P

toc