function [ret]=passenger_density()
distribute=xlsread('C:\Users\song\Desktop\CUMCM2017\passenger.xlsx');
sum=0;
data=distribute;
cnt=zeros(835,1);
for i=1:835
    sum=0;
    for j=1:835
        if data(i,j)<1.5299
            sum=sum+1;
        end
    end
    cnt(i,1)=sum;
end
xlswrite('C:\Users\song\Desktop\CUMCM2017\ÈÎÎñÃÜ¶È.xlsx',round(cnt/10));
