%% compute the limit cycle in two dimension

clear,clc;
tic

%% initialization
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
%% carbon system: vector filed (f,g)
f = @(x,y,z)    -(x/(e+x))*(a/(b+c*y))*y+r*y+(beta*lamda*(y^2)/(mu^2+y^2))*z+gamma*d*(z^2)+k*(N_0-x);
g = @(x,y,z)    (x/(e+x))*(a/(b+c*y))*y-r*y-(lamda*(y^2)/(mu^2+y^2))*z-(s+k)*y;
h = @(x,y,z)    (alpha*lamda*(y^2)/(mu^2+y^2))*z-d*(z^2);
%% Compute the trajectory using Euler scheme
T = 100;
N = 1e5;
dt = T/N;
t = 0 : dt: T;
%83.5819254216627;2290.50827400743
%[83.6339868155324;2290.44660430012]
%鍒濆鐐圭殑閫夊彇浼氬奖鍝嶆瀬闄愮幆锛宯u=0.9鐨勬椂鍊欓渶瑕佹敼鍙樹竴涓嬪垵鍊?
 x(1) =0.8;
 y(1) =0.3;
 z(1) =0.31;

%92.6715592143817;2281.91346961487
for i = 1 : N
    x(i+1) = x(i) + f(x(i),y(i),z(i))*dt;
    y(i+1) = y(i) + g(x(i),y(i),z(i))*dt;
    z(i+1) = z(i) + h(x(i),y(i),z(i))*dt;
end
 % plot(x,y,'b');
  hold on
  stable=[x(1:5000);y(1:5000);z(1:5000)];
%  plot(t,x,'r');hold on
%  plot(t,y,'g')
p = 1;
LCV(1,1) = x(end);
LCV(2,1) = y(end);
LCV(3,1) = z(end);
while p > 0 || i==0
    i = i - 1;
    LCV(1,N-i) = x(i+1);
    LCV(2,N-i) = y(i+1);
    LCV(3,N-i) = z(i+1);
    if norm(LCV(:,1)-LCV(:,N-i))<0.5&& t(end)-t(i)>35;
        p=0
    end
end
LCV(:,N-i+1)=[x(end);y(end);z(end)];

plot3(LCV(1,:),LCV(2,:),LCV(3,:),'r')
save('LimitCycle_1.txt','LCV','-ascii','-double');
%save('LimitCycle_11.mat','LCV');
toc
hold on
