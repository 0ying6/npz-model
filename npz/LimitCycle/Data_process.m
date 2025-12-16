%% compute the stable point
clear,clc
load('shiyan9.mat');
SED=P;
u = load('LimitCycle_1.txt');
plot3(u(1,:),u(2,:),u(3,:),'r');
hold on
 %plot(SED(3,:),SED(4,:),'*');hold on

[~,n] = size(SED);
k = 0;
SDE11 = zeros(n,1);
%SDE12 = zeros(n,1);
for j = 1 : n
      if  min(sum((u-SED([4,5,6],j)).^2)) <1e-6
                 k = k + 1;
                 SDE11(k)=j;
      end
end

SDE1 = SED([1,2,3,4,5,6],SDE11(1:k));
%plot(SDE1(3,:),SDE1(4,:),'*');  

l=5; %稀疏度>0，越大越密
SDE21=SDE1([4,5,6],:);
[~, n] = size(SDE21);
xmin=min(SDE21(1,:));xmax=max(SDE21(1,:));
ymin=min(SDE21(2,:));ymax=max(SDE21(2,:));
zmin=min(SDE21(3,:));zmax=max(SDE21(3,:));
SDE21=SDE21-[xmin,ymin,zmin]';
SDE21=floor(SDE21./[(xmax-xmin)/l,(ymax-ymin)/l,(zmax-zmin)/l]');
[~,ia,~]=unique(SDE21',"rows","stable");
SDE2=SDE1(:,ia);
%save sparseshiyan9.mat SDE2
plot3(SDE2(4,:),SDE2(5,:),SDE2(6,:),'*');   



