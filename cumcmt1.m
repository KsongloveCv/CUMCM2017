function cumcmt1
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here data process
%子函数    [geoid,msg] = geoidtst(geoid)
%          shortdistance(lat1, lon1, lat2, lon2, ellipsoid)
tic
disp(sprintf('正在载入相关数据,请耐心等待...'));
%%数据准备
distribute1=xlsread('C:\Users\song\Desktop\CUMCM2017\d1.xls');
distribute2=xlsread('C:\Users\song\Desktop\CUMCM2017\d2.xlsx');
number=length(distribute1);pi=3.1415926;
disp(['number is ',num2str(number)]);
x=distribute1(:,1);lati=x*(pi/180);
y=distribute1(:,2);loni=y*(pi/180);
distance=zeros(number,number);
ellipsoid=geoidtst(almanac('earth','ellipsoid'));
for i=1:number
    distance(i,i+1:number) = shortdistance(lati(i)*ones(length(i+1:number),1), loni(i)*ones(length(i+1:number),1),lati(i+1:number), loni(i+1:number),ellipsoid);
    distance(i,1)=distance(i,1)*6371*1000*2*pi/360;  
end
dis=(distance+distance')/10;
distribute1=distribute1(:,1:2);
distribute2=distribute2(:,1:2);
diss=[distribute2;distribute1];
number2=1877+835;
disp(['number2 is ',num2str(number2)]);
xx=diss(:,1);latii=xx*(pi/180);
yy=diss(:,2);longii=yy*(pi/180);
distance2=zeros(number2,number2);
ellipsoid = geoidtst(almanac('earth','ellipsoid'));
for i=1:number2
    distance2(i,i+1:number2) = shortdistance(latii(i)*ones(length(i+1:number2),1), longii(i)*ones(length(i+1:number2),1),latii(i+1:number2), longii(i+1:number2),ellipsoid);
    distance2(i,1)=distance2(i,1)*6371*1000*2*pi/360;  
end
disfarm=(distance2+distance2')/10;
dis2=disfarm(836:2712,1:835);
xlswrite('C:\Users\song\Desktop\CUMCM2017\passenger.xlsx',dis);
xlswrite('C:\Users\song\Desktop\CUMCM2017\guanxi.xlsx',dis2);
toc



