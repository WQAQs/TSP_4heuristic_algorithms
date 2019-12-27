clc,clear 
load sj.txt
x=sj(:,1:2:8);x=x(:); %x为经度
y=sj(:,2:2:8);y=y(:); %y为纬度
sj=[x y]; %出发点
d1=[70,40]; 
sj=[d1;sj;d1]; %把出发点，终点添加到表中
sj=sj*pi/180; %经纬度换算成距离（单位m）
%距离矩阵d
d=zeros(102); 
for i=1:101 
 for j=i+1:102 
 temp=cos(sj(i,1)-sj(j,1))*cos(sj(i,2))*cos(sj(j,2))+sin(sj(i,2))*sin(sj(j,2)); 

 d(i,j)=6370*acos(temp); 
 end 
end
d=d+d'; % d'表示矩阵d的转置，因为距离矩阵是一个对称阵
S0=[];% S0是访问路线
Sum=inf; % Sum是S0访问路线对应的访问总费用
rand('state',sum(clock)); 
for j=1:1000 % 初始化一个比较好的访问路线S0
    S=[1 1+randperm(100),102]; 
    temp=0;
    for i=1:101 
        temp=temp+d(S(i),S(i+1)); 
    end 
    if temp<Sum 
        S0=S;Sum=temp; 
    end 
end
e=0.1^30;L=30000;at=0.999;T=1; 
sum_iter=[]%保存每一次迭代的代价函数
%退火过程
for k=1:L 
    %产生新解，随机化产生两个数c1，c2（即在S0中要交换的两个访问点）
    c=2+floor(100*rand(1,2));
    c=sort(c);
    c1=c(1);c2=c(2); 
    %计算代价函数
    df=d(S0(c1-1),S0(c2))+d(S0(c1),S0(c2+1))-d(S0(c1-1),S0(c1))-d(S0(c2),S0(c2+1)); %交换S0中第c1和第c2个访问点后的总费用-之前的总费用
    if df<0 %若总费用降低，则直接更新访问路线S0
      S0=[S0(1:c1-1),S0(c2:-1:c1),S0(c2+1:102)]; 
      Sum=Sum+df; 
    elseif exp(-df/T)>rand(1)  %若总费用升高，则以一定的概率exp(-df/T)更新访问路线S0
      S0=[S0(1:c1-1),S0(c2:-1:c1),S0(c2+1:102)]; 
      Sum=Sum+df; 
    end 
    T=T*at; 
    if T<e 
      break; 
    end 
    sum_iter=cat(1,sum_iter,Sum);
end
figure(1)
plot(sum_iter)
%输出巡航路径及路径长度
S0;Sum;
% 绘制巡航路径
res = [];
for i=1:102
    index = S0(i)
    resx = sj(index,1)
    resy = sj(index,2)
    B = [resx resy]
    res = cat(1,res,B)
end
X = res(:,1)*180/pi
Y = res(:,2)*180/pi
figure(2)
plot(X,Y,'-ro','linewidth',1,'markeredgecolor','b','markerfacecolor','0.49,1,0.65','markersize',3);
