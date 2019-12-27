clc,clear 
load sj.txt
x=sj(:,1:2:8);x=x(:); %xΪ����
y=sj(:,2:2:8);y=y(:); %yΪγ��
sj=[x y]; %������
d1=[70,40]; 
sj=[d1;sj;d1]; %�ѳ����㣬�յ���ӵ�����
sj=sj*pi/180; %��γ�Ȼ���ɾ��루��λm��
%�������d
d=zeros(102); 
for i=1:101 
 for j=i+1:102 
 temp=cos(sj(i,1)-sj(j,1))*cos(sj(i,2))*cos(sj(j,2))+sin(sj(i,2))*sin(sj(j,2)); 

 d(i,j)=6370*acos(temp); 
 end 
end
d=d+d'; % d'��ʾ����d��ת�ã���Ϊ���������һ���Գ���
S0=[];% S0�Ƿ���·��
Sum=inf; % Sum��S0����·�߶�Ӧ�ķ����ܷ���
rand('state',sum(clock)); 
for j=1:1000 % ��ʼ��һ���ȽϺõķ���·��S0
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
sum_iter=[]%����ÿһ�ε����Ĵ��ۺ���
%�˻����
for k=1:L 
    %�����½⣬���������������c1��c2������S0��Ҫ�������������ʵ㣩
    c=2+floor(100*rand(1,2));
    c=sort(c);
    c1=c(1);c2=c(2); 
    %������ۺ���
    df=d(S0(c1-1),S0(c2))+d(S0(c1),S0(c2+1))-d(S0(c1-1),S0(c1))-d(S0(c2),S0(c2+1)); %����S0�е�c1�͵�c2�����ʵ����ܷ���-֮ǰ���ܷ���
    if df<0 %���ܷ��ý��ͣ���ֱ�Ӹ��·���·��S0
      S0=[S0(1:c1-1),S0(c2:-1:c1),S0(c2+1:102)]; 
      Sum=Sum+df; 
    elseif exp(-df/T)>rand(1)  %���ܷ������ߣ�����һ���ĸ���exp(-df/T)���·���·��S0
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
%���Ѳ��·����·������
S0;Sum;
% ����Ѳ��·��
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
