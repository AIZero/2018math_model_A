clear;clc;
%�����ȴ����ʣ��ܶȣ�����
v1=0.082;v2=0.37;v3=0.045;v4=0.028;
p1=300;p2=862;p3=74.2;p4=1.18;
c1=1377;c2=2100;c3=1726;c4=1005;
%���嵥λʱ��t����λ����x
t=0.01;x=0.2*10^-3;
%��������ɢ���뵥λʱ�估��λ����ƽ�������ĳ˻�
h1=(v1/(p1*c1))*t/x/x;h2=(v2/(p2*c2))*t/x/x;
h3=(v3/(p3*c3))*t/x/x;h4=(v4/(p4*c4))*t/x/x;
y=0.01;
o=5.669*10^-8;
e=32.9253;
%��ʼ����¼����
accord=[];
%������٣�jΪ�ڶ���ĺ��
for j=0.6:0.2:25
    a=ceil((0.6+j+3.6+5)/0.2);   %�����Ȳ���
    T=ones(a+1,360000)*37;       %��ʼ���¶Ⱦ���
    T(1,:)=65;                   %��ʼ���߽��¶�
%���㣨��ȣ�ʱ�䣩ʱ���¶�T
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
%��ȡƤ������¶���ʱ��ı仯����ѡȡ����������ֵ��¼��accord�У����е�һ��Ϊ��ȣ��ڶ���Ϊ����¶�����������Ϊ�¶���ʱ��仯��ֵ
tt=T(floor((0.6+j+3.6+5)/0.2),2:end);
max_num=max(tt);
len=length(find(tt>44));
if max_num<47&&len<30000
    b=[j,max_num,tt];
    accord=[accord;b];
end
end
%��ѡ����������ֵ
ttt=[];
last_accord=[];
for i=1:156
    max1=max(accord(i,:));
    ttt=[ttt,max1-accord(i,end)];
end
bala=find(ttt==0);
for i=1:length(bala)
    bal=bala(i);
    last_accord=[last_accord;accord(bal,:)];
end
%plot(last_accord(:,1),last_accord(:,2));
%������Ž�
last_accord(1,1)
