clc;clear all;
disp(sprintf('正在载入相关数据,请耐心等待...'));
x=xlsread('F:\备份\B201727004013_杨楷文_宋康_姚晰月\附件一.xls','Sheet1');
x=x(1:835,1:3);
p=x';%输入矩阵
y=xlsread('F:\备份\B201727004013_杨楷文_宋康_姚晰月\附件一.xls','Sheet1');
y=y(1:835,4);
t=y';
[pn,minp,maxp,tn,mint,maxt]=premnmx(p,t); 
dx=[-1,1;-1,1;-1,1]; %归一化处理后最小值为-1，最大值为 1
%BP 网络训练
net=newff(dx,[3,5,1],{'tansig','tansig','purelin'},'traingdx'); % 建 立 模型，并用梯度下降法训练．
net.trainParam.show=10000; %10000 轮回显示一次结果
net.trainParam.Lr=0.05; %学习速度为 0.05
net.trainParam.epochs=20000; %最大训练轮回为 20000 次
net.trainParam.goal=0.065*10^(-2); %均方误差
net=train(net,pn,tn); %开始训练，其中 pn,tn 分别为输入输出样本
%利用原始数据对 BP 网络仿真
an=sim(net,pn); %用训练好的模型进行仿真
a=postmnmx(an,mint,maxt); % 把仿真得到的数据还原为原始的数量级；
pnew=xlsread('F:\备份\B201727004013_杨楷文_宋康_姚晰月\附件一.xls','Sheet2');
pnew=pnew(1:835,1:3);
pnew=pnew'
pnewn=tramnmx(pnew,minp,maxp); %利用原始输入数据的归一化参数对新数据进行归一化；
anewn=sim(net,pnewn); %利用归一化后的数据进行仿真；
anew=postmnmx(anewn,mint,maxt) %把仿真得到的数据还原为原始的数量级；
anew=round(anew)';