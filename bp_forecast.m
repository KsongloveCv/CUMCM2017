clc;clear all;
disp(sprintf('���������������,�����ĵȴ�...'));
x=xlsread('F:\����\B201727004013_���_�ο�_Ҧ����\����һ.xls','Sheet1');
x=x(1:835,1:3);
p=x';%�������
y=xlsread('F:\����\B201727004013_���_�ο�_Ҧ����\����һ.xls','Sheet1');
y=y(1:835,4);
t=y';
[pn,minp,maxp,tn,mint,maxt]=premnmx(p,t); 
dx=[-1,1;-1,1;-1,1]; %��һ���������СֵΪ-1�����ֵΪ 1
%BP ����ѵ��
net=newff(dx,[3,5,1],{'tansig','tansig','purelin'},'traingdx'); % �� �� ģ�ͣ������ݶ��½���ѵ����
net.trainParam.show=10000; %10000 �ֻ���ʾһ�ν��
net.trainParam.Lr=0.05; %ѧϰ�ٶ�Ϊ 0.05
net.trainParam.epochs=20000; %���ѵ���ֻ�Ϊ 20000 ��
net.trainParam.goal=0.065*10^(-2); %�������
net=train(net,pn,tn); %��ʼѵ�������� pn,tn �ֱ�Ϊ�����������
%����ԭʼ���ݶ� BP �������
an=sim(net,pn); %��ѵ���õ�ģ�ͽ��з���
a=postmnmx(an,mint,maxt); % �ѷ���õ������ݻ�ԭΪԭʼ����������
pnew=xlsread('F:\����\B201727004013_���_�ο�_Ҧ����\����һ.xls','Sheet2');
pnew=pnew(1:835,1:3);
pnew=pnew'
pnewn=tramnmx(pnew,minp,maxp); %����ԭʼ�������ݵĹ�һ�������������ݽ��й�һ����
anewn=sim(net,pnewn); %���ù�һ��������ݽ��з��棻
anew=postmnmx(anewn,mint,maxt) %�ѷ���õ������ݻ�ԭΪԭʼ����������
anew=round(anew)';