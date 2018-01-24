%ÁéÃô¶È·ÖÎö
cc_percent = [];
for i = 50:149
    p1 = floor(i*(1+0.05*rand()));
    p2 = floor(i*(1+0.05*rand()));
    ans1 = abs((G(i)*(0.75*lambda(i)+1)-G(i)*0.25*phi(i))/F(i)-1);
    ans2 = abs((G(p1)*(0.75*lambda(p2)+1)-G(p1)*0.25*phi(p1))/F(p2)-1);
    cc_percent(i-49) = abs((ans1-ans2)/ans1);
end