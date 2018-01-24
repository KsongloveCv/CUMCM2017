clear 
clc 
fitnessfcn=@Fun;
 nvars=4; 
 lb=[0.5,0.5,0.5,0.5]; 
 ub=[3.5,3.5,3.5,3.5];
 A=[];b=[]; 
 Aeq=[];beq=[]; 
options=gaoptimset('paretoFraction',0.3,'populationsize',100,'generations',300,'stallGenLimit',200,'TolFun',...
    1e-10,'PlotFcns',@gaplotpareto);


[x,fval]=gamultiobj(fitnessfcn,nvars,A,b,Aeq,beq,lb,ub,options)
