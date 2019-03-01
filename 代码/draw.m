%画图找规律
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
d2=[0.6 10 15 20 25];     %给定第二层d2的厚度
%d4=[3.5 3.9 4.2 5.1];  %给定第四层d3的厚度
%TTT=[55 65 70 80];     %给定外界环境的温度
accord=[];
for z=1:1:4
j=d2(z);     %第二层厚度
%j=6;
%j=6;
kl=5;        %第四层厚度
%k1=d4(z)
%k1=5;
a=ceil((0.6+j+3.6+kl)/0.2);
T=ones(a+1,360000)*37;
T(1,:)=80;       %外界环境温度
%T(1,:)=80;
%T(1,:)=TTT(z);
for  k=1:360000
    for i=2:a
      if(i<=4)
            T(i,k+1)=h1*(T(i-1,k)+T(i+1,k))+(1-2*h1)*T(i,k);
            if i==4
                T(i,k+1)=(h1*T(i-1,k)+h2*T(i+1,k))/(h1+h2);
            else
            T(i,k+1)=h1*(T(i-1,k)+T(i+1,k))+(1-2*h1)*T(i,k)+y*2100/(p1*c1)*exp(-y*i*x);
            end;
       else
            if(4<i&i<=(4+j/0.2))
                if i==(4+j/0.2)
                    T(i,k+1)=(h2*T(i-1,k)+h3*T(i+1,k))/(h3+h2);
                else
                T(i,k+1)=h2*(T(i-1,k)+T(i+1,k))+(1-2*h2)*T(i,k);
                end
            else
                if((4+j/0.2)<i&i<=(4+j/0.2+19))
                   if i==(4+j/0.2+18)
                        T(i,k+1)=(h3*T(i-1,k)+h4*T(i+1,k))/(h3+h4);
                    else
                    T(i,k+1)=h3*(T(i-1,k)+T(i+1,k))+(1-2*h3)*T(i,k);
                    end
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
%zxc=T(a,:);
%qw=length(find(zxc>44));
%bb=[bb,qw];
%gg=T(:,360000);
accord=T(a,:);
hold on
%plot(1:a,gg(1:a,1));
plot(1:length(accord)/100+1,accord(1:100:length(accord)));
end
xlabel('时间  t/s');
ylabel('温度  °C');
