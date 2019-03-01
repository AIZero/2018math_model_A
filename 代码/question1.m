clear;clc;
%定义厚度,热传导率，密度，比热
d1=0.6*10^-3;d2=6*10^-3;d3=3.6*10^-3;d4=5*10^-3;
v1=0.082;v2=0.37;v3=0.045;v4=0.028;
p1=300;p2=862;p3=74.2;p4=1.18;
c1=1377;c2=2100;c3=1726;c4=1005;
%定义单位时间t，单位长度x
t=0.01;x=0.2*10^-3;
%计算热扩散率与单位时间及单位长度平方倒数的乘积
h1=(v1/(p1*c1))*t/x/x;h2=(v2/(p2*c2))*t/x/x;
h3=(v3/(p3*c3))*t/x/x;h4=(v4/(p4*c4))*t/x/x;
%初始化温度矩阵，行数为位移步数，列数为时间的步数
d=d1+d2+d3+d4;
a=ceil(d/x)+1;  
T=ones(a,540000)*37;
T(1,:)=75;  %初始化边界温度
y=0.01;
o=5.669*10^-8;
e=32.9253;
%计算（厚度，时间）时的温度T
for  k=1:540000
    for i=2:76
        if(i<=4)
            %T(i,k+1)=h1*(T(i-1,k)+T(i+1,k))+(1-2*h1)*T(i,k);
            if i==4
                T(i,k+1)=(h1*T(i-1,k)+h2*T(i+1,k))/(h1+h2);
            else
                T(i,k+1)=h1*(T(i-1,k)+T(i+1,k))+(1-2*h1)*T(i,k)+y*2100/(p1*c1)*exp(-y*i*x);
            end
        else
            if(4<i&i<=34)
                if i==34
                    T(i,k+1)=(h2*T(i-1,k)+h3*T(i+1,k))/(h3+h2);
                else
                    T(i,k+1)=h2*(T(i-1,k)+T(i+1,k))+(1-2*h2)*T(i,k);
                end
            else
                if(34<i&i<=53)
                    if i==52
                        T(i,k+1)=(h3*T(i-1,k)+h4*T(i+1,k))/(h3+h4);
                    else
                        T(i,k+1)=h3*(T(i-1,k)+T(i+1,k))+(1-2*h3)*T(i,k);
                    end
                else
                    %T(i,k+1)=h4*(T(i-1,k)+T(i+1,k))+(1-2*h4)*T(i,k);
                    T(i,k+1)=-(o*e*(T(52,k+1)^4-T(76,k)^4)+5*(T(52,k+1)-T(76,k)))*x/v4+T(i-1,k+1);
                    
                end
            end
        end 
    end
end
A=xlsread('...\2018CUMCM-Problem-ABCD\2018-A-Chinese\CUMCM-2018-Problem-A-Chinese-Appendix.xlsx',2,'A3:B5403');   %读取数据
t=[1:5400;T(76,1:100:540000)]';    %读取皮肤外层的温度及其对应的时间
figure(1)
plot(A(2:end,1),A(2:end,2));       %附件二的数据图
xlabel('时间  t/s');
ylabel('温度  °C');
figure(2);
plot(A(2:end,1),A(2:end,2));       %附件二的数据图
xlabel('时间  t/s');
ylabel('温度  °C');
hold on
plot(t(:,1),t(:,2));               %皮肤外层的温度及其对应的时间图
figure(3);
[xx,yy]=meshgrid(0:0.2:15.2,1:5400);
mesh(xx',yy',T(:,1:100:540000));     %整个过程的温度分布图
xlabel('厚度  mm');
ylabel('时间  t/s');
zlabel('温度  °C');
figure(4)                          %各层交界处温度与时间的关系图
plot(1:5400,T(4,1:100:540000),1:5400,T(34,1:100:540000));
hold on
plot(1:5400,T(52,1:100:540000),1:5400,T(76,1:100:540000));
xlabel('时间  t/s');
ylabel('温度  °C');
TT=T(:,1:100:540000)';
xlswrite('...\problem1.xlsx',TT,1,'A1:BX5400');
