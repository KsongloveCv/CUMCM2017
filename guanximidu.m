function [ret]=guanxi_density()
distribute=xlsread('C:\Users\song\Desktop\CUMCM2017\guanxi.xlsx');
data=distribute';
cnt=zeros(835,1);
for i=1:835
    sum=0;
    for j=1:1877
        %d=sqrt(s^2*exp(-lambda)*k/a*pi);һ��k=0.8,sΪ�������,lambda�˾�GDP,a���˿�-->d=1.5299
        if data(i,j)<1.5299
            sum=sum+1;
        end
    end
    cnt(i,1)=sum;
end
xlswrite('C:\Users\song\Desktop\CUMCM2017\��ϵ�ܶ�.xlsx',round(cnt/10));
