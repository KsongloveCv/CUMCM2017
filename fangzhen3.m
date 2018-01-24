clear all;clc;
tic;
%随机生成500个订单  passenger   随机生成800个会员  huiyuan 随机生成50个中心点center
huiyuan_data=[];data=[];data1=[];
passenger=2066;huiyuan=3000;center=5;
p=12;
passenger_lati=xlsread('C:\Users\song\Desktop\CUMCM2017\附件三.xls','A1:A2066');%订单数据
passenger_loni=xlsread('C:\Users\song\Desktop\CUMCM2017\附件三.xls','B1:B2066');
%转换距离
huiyuan_lati=22+p*rand(huiyuan,1);
huiyuan_loni=112+p*rand(huiyuan,1);
center_lati=22+p*rand(center,1);
center_loni=112+p*rand(center,1);
plot(passenger_lati,passenger_loni,'bo','MarkerSize', 5);
%axis([20 35 111 124]);
hold on
plot(huiyuan_lati,huiyuan_loni,'g*','MarkerSize', 5);
%axis([20 35 111 124]);
hold on
plot(center_lati,center_loni,'mp','MarkerSize', 8);
%axis([20 35 110 125]);
xlabel('纬度')
ylabel('经度')
legend('任务','会员','中心点');
for ii=110:125
    x=20:1:35;y= 0*x+i;
    if mod(ii,5) == 0 
       plot(x,y,'k:');
    end
end
for jj=21:35
    x=110:125;y=0*x+jj; 
    if mod(jj,5) == 0 
       plot(y,x,'k:');
    end
end
passenger_data=[passenger_lati passenger_loni ones(passenger,1)];
huiyuan_data=[huiyuan_lati huiyuan_loni ones(huiyuan,1)];
center_data=[center_lati center_loni ones(center,1)];
data=[passenger_data;center_data;huiyuan_data];
dis=[];dis=data;
nfarms=passenger+center;
disp(['nfarms is ',num2str(nfarms)]);
k=size(dis,1);
dis_x=repmat(dis(:,1),1,k);
dis_y=repmat(dis(:,2),1,k);
dist_x=repmat(dis(:,1)',k,1);
dist_y=repmat(dis(:,2)',k,1);
distance=sqrt((dis_x - dist_x).^2+(dis_y - dist_y).^2);
%提取矩阵
passenger_dis=distance(1:passenger,1:passenger);
center_dis=distance(passenger+1:passenger+center,passenger+1:passenger+center);
%中心点到任务点的距离     
dis_cp=distance(passenger+1:passenger+center,1:passenger); 
%提取人到每个任务的距离矩阵
dis_hr=distance(passenger+center+1:end,1:passenger+center); 
%为了简化 使每个包中点个数相同 每个包放多少点 number
number=10;[m,n]=size(dis_cp);
[Y,I]=sort(dis_cp,2,'ascend');   %按照行进行升序排列
RowCol.value=Y(:,1:number);      %结构体RowCol，存放数值
%RowCol.row=repmat((1:N).',number);  %行坐标值
RowCol.col=I(1,1:number);        %列坐标值
for m=1:center
    aa{m}=I(m,:);
end
for mm=2:center
    for n=1:number
        for nn=mm:center    
        aa{nn}(find(aa{nn}==aa{mm-1}(n)))=[];
        end  
    end
end
center_good=[aa{1}(1:number); aa{2}(1:number); aa{3}(1:number); aa{4}(1:number); aa{5}(1:number)];%5*10
%划分之后 如何去确定去做不做
%确定人去做任务  任务落到意愿之内 优先选择一个离得最近的 然后去做 花费时间为 t走+t做
%假设打包项目和未打包项目对人的吸引力相同  只要落入会员范围之内就会去选择 然后开始做
%假设完成一个任务所需时间为6000*rand();单位为s
speed_do=6000*rand();
r_valid=500*rand();%距离意愿 单位为m
speed=0;
%人走过去的时间  speed_go  假设人的速度为2m/s
%人如果做得是打包的任务 那么要将包内任务做完才能去接下一个任务
all_xuhao=[1:passenger+center]';
all_zuobiao=[passenger_data;center_data];
passenger_xuhao=all_xuhao(1:passenger);
center_xuhao_1=[101 102 103 104 105];
center_xuhao_2=[aa{1}(1:number);aa{2}(1:number);aa{3}(1:number);aa{4}(1:number);aa{5}(1:number)]'; 
center_xuhao=[center_xuhao_1;center_xuhao_2];
%打包成五类
for i=1:length(center_xuhao)
    passenger_xuhao(find(passenger_xuhao==center_xuhao(i)))=[]; 
end
bao_nei=center_xuhao;
bao_wai=passenger_xuhao;
%hui_dis=zeros(huiyuan,passenger+center)
%for i=1:huiyuan
%    hui_dis(i,:)=sort(dis_hr(i,:));%200*85   200个会员距离每个任务的距离
%end
%假设会员编号 按照信誉值高低排序  初始化
c=0;
data1=[passenger_data;center_data];
huiyuan_data(:,3)==0;data1(:,3)==0;renwu_data=data1;
all_K = [];xmax = 111*cos(pi*34/180)*1.4;
ymax = 0.7*111;
%包内的距离用常数去量化 统一为50
for i=1:24*3600 %时间
    r_valid=500*rand();
    %首先更新会员状态
    lc =huiyuan_data(:,3)-0.1;
    lc(lc<=0)=0; %都未完成
    huiyuan_data(:,3)=lc;
    %会员随机一个方向前进0.1
    valid_lines=find(lc==0);
    for m=1:length(valid_lines)
        k = valid_lines(m);
        for e=1:100
             degree = 2*pi*rand();%出行方向
             %下一状态位置变化
             new_x = huiyuan_data(k,1) + 0.1.*cos(degree);
             new_y = huiyuan_data(k,2) + 0.1.*sin(degree);
             if(new_x>=0&& new_x<= r_valid&& new_y>=0 && new_y<=r_valid)
             huiyuan_data(k,1:2) = [new_x,new_y]; %位置变化
                 break
             end
        end
    end
     for j=1:passenger+center
         
         %if renwu_data(j,3)==1;  %已经被完成
          %   continue;
         %end
         ren_num=find(dis_hr(m,j)<r_valid);
         if isempty(ren_num) %意愿内没有任务，下一位会员
             continue;
         else  %第几个货物
             c=j;
             if ismember(c,bao_nei(:,1));
                 for ii=1:11
                     renwu_data(bao_nei(ii,1),3)=0;
                 end
                 speed=speed_do*11+25;
             elseif ismember(c,bao_nei(:,2));
                 for ii=1:11
                     renwu_data(bao_nei(ii,2),3)=0;
                 end
                 speed=speed_do*11+25;
             elseif ismember(c,bao_nei(:,3));
                 for ii=1:11
                     renwu_data(bao_nei(ii,3),3)=0;
                 end
                 speed=speed_do*11+25;
             elseif ismember(c,bao_nei(:,4));
                 for ii=1:11
                     renwu_data(bao_nei(ii,4),3)=0;
                 end
                 speed=speed_do*11+25;
              elseif ismember(c,bao_nei(:,5));
                 for ii=1:11
                     renwu_data(bao_nei(ii,5),3)=0;
                 end
                 speed=speed_do*11+25;
             else
                 renwu_data(j,3)=0;
             end
         end
     end
end   
toc    
            
            
            
