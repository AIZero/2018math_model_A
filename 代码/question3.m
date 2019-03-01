clear;clc;
%定义热传导率，密度，比热
v1=0.082;v2=0.37;v3=0.045;v4=0.028;
p1=300;p2=862;p3=74.2;p4=1.18;
c1=1377;c2=2100;c3=1726;c4=1005;
%定义单位时间t，单位长度x
t=0.01;x=0.2*10^-3;
%计算热扩散率与单位时间及单位长度平方倒数的乘积
h1=(v1/(p1*c1))*t/x/x;h2=(v2/(p2*c2))*t/x/x;
h3=(v3/(p3*c3))*t/x/x;h4=(v4/(p4*c4))*t/x/x;
y=0.01;
o=5.669*10^-8;
e=32.9253;
%初始化记录矩阵
accord=[];
%进行穷举，j为第二层的厚度，
for j=0.6:0.2:25
    for four=4:0.1:5.6
    a=ceil((0.6+j+3.6+four)/0.2);       %计算厚度步数
    T=ones(a+1,180000)*37;              %初始化温度矩阵
    T(1,:)=80;                          %初始化边界温度
for  k=1:180000
    for i=2:a
        if(i<=4)
            %T(i,k+1)=h1*(T(i-1,k)+T(i+1,k))+(1-2*h1)*T(i,k);
            T(i,k+1)=h1*(T(i-1,k)+T(i+1,k))+(1-2*h1)*T(i,k)+y*2100/(p1*c1)*exp(-y*i*x);
       else
            if(4<i&i<=(4+j/0.2))
                T(i,k+1)=0.55*h2*(T(i-1,k)+T(i+1,k))+(1-2*0.55*h2)*T(i,k);
            else
                if((4+j/0.2)<i&i<=(4+j/0.2+24))
                    T(i,k+1)=0.55*h3*(T(i-1,k)+T(i+1,k))+(1-2*0.55*h3)*T(i,k);
                else
                    %T(i,k+1)=h4*(T(i-1,k)+T(i+1,k))+(1-2*h4)*T(i,k);
                    T(i,k+1)=-(o*e*(T(floor(4+j/0.2+18),k+1)^4-T((a),k)^4)+5*(T(floor(4+j/0.2+18),k+1)-T((a),k)))*x/v4+T(i-1,k+1);
                    if T(i,k+1)<=37
                        T(i,k+1)=37;
                    end
                end
            end
        end 
    end
end
%获取皮肤外层温度随时间的变化，并选取符合条件的值纪录于accord中，其中第一行为第二层厚度，第二行为第四层厚度，第三行为最大温度数，第四行为温度随时间变化的值
tt=T(floor((0.6+j+3.6+four)/0.2),2:end);
max_num=max(tt);
len=length(find(tt>44));
if max_num<47&&len<30000
    b=[j,four,max_num,tt];
    accord=[accord;b];
end
    end
end
%挑选满足条件的值
ttt=[];
accord1=[];
for i=1:692
    max1=max(accord(i,:));
    ttt=[ttt,max1-accord(i,end)];
end
bala=find(ttt==0);
for i=1:length(bala)
    bal=bala(i);
    accord1=[accord1;accord(bal,:)];
end

last_accord=[];
for i=1:333
    if accord1(i,3)~=37
        last_accord=[last_accord;accord1(i,:)];
    end
end