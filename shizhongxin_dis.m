function [ret]=shizhongxin_dis()
distribute=xlsread('C:\Users\song\Desktop\CUMCM2017\zhongxin.xlsx');
x=zeros(835,1);
for i=1:835
    data=distribute(i,:)
    x(i,1)=min(data);
end
xlswrite('C:\Users\song\Desktop\CUMCM2017\ÊÐÖÐÐÄ¾àÀë.xlsx',x);