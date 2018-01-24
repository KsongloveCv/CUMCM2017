clear all;clc;
p=xlsread('C:\Users\song\Desktop\CUMCM2017\问题一数据.xlsx');
cnt
route=randperm(835);
z=route(:,1:100);
h=z';
for i=1:835
    for j=1:100
        if(i==h(j))
            cnt(j,:)=p(i,:);
        end
    end
end
cnt1=round(cnt(:,1)/10);
cnt2=round(cnt(:,2)/10);
cnt3=cnt(:,3);cnt4=cnt(:,4);
cntt=[cnt1 cnt2 cnt3 cnt4];


