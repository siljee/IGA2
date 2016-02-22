close all
addpath ../; addpath ../Basis/; addpath ../HelpFunctions/
a = -1; b = 1; m=3;
xb = linspace(knotVec(1),knotVec(end),100)

func = @(x) (x-2).^2;
%func = @(x) sin(pi*x)+1
figure(1)
%plot(x,1-exp(-5*abs(x)),'--k')
plot(x,func(x),'--k')
hold on 
knotVec = [0,0,0,0,1,2,3,4,4,4,4];

maxp = 3;
for p=3:maxp
    %t = [a*ones(1,p),a:1/m:b,b*ones(1,p)]
    %knotVec = [a*ones(1,p),a:1/m:0,0*ones(1,p-2),0:1/m:b,b*ones(1,p)];
    knotVec = makeUniformKnotVector(p,4)
    %ts = (t(2:end-2) +t(3:end-1))/2;
    %t = findGrevillePoints(knotVec,p);

    %t = 0;
    %for i = 1:p
    %    t = t + knotVec(1+i:end-(p+1-i));
    %end
    %t = t/p;
    t = findGrevillePoints(knotVec,p)
    
    
    N = BsplineBasis(knotVec,p,linspace(knotVec(1),knotVec(end),100));
    
    x = xb*(b-a)+a
    t = t *(b-a)+a;
    %%%%%%%%%%%%%%%%%%%%
    % Plot
    %%%%%%%%%%%%%%%%%%%%
    %Qf = sqrt(2)*sin(pi*ts/2)*N;
    
    figure(1)
    Qf= (func(t))*N;
    plot(x,Qf,'r')
    plot(t,func(t),'*-')
    
    figure(2)
    subplot(ceil(sqrt(maxp)),ceil(sqrt(maxp)),p)
    plot(x,func(x)-Qf,'g')
  
   
end

figure(1)
legend('exact','p=1','p=2','p=3','p=4','p=5','p=6')

%plot(x,sin(pi*x/2),'--k', 'linewidth',2)


%%%%%%%%%%%%%%%%%%%%%%%%% GENRAL SPLINE INTERPOLATION %%%%%%%%%%%%%
