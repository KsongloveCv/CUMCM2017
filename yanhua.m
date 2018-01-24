%Ö÷³ÌÐò
clear
%y-x
for i=0:0.1:1
    for j=0:0.1:1
        [T,Y]=ode45('differential',[5 3.1],[i j]);
        figure(1)
        grid on
        plot(Y(:,1),Y(:,2));
        hold on
    end
end
