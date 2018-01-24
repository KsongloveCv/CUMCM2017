%�ݻ�����
% drawItems.m
% Copyright 2016
% ��ѭ GPL Э�飨������ί�᲻�ܴ�Э�����ƣ�

function  drawItems ( AB,Q,init,li,lj,k)
Aa = AB(1,1);
Ab = AB(1,2);
Ba = AB(2,1);
Bb = AB(2,2);
for i=0:0.1:li
    for j=0:0.1:lj
        [T,Y]=ode45('fx',init,[i j],[],Q,Aa,Ab,Ba,Bb);
        grid on
        figure(k)
        hold on
        plot(Y(:,1),Y(:,2));
    end
end
