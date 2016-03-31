close all; addpath ../; addpath ../Basis/; addpath ../HelpFunctions/
a = -1; b = 1;
el = 5; p = 3;

func = @(x) (x).^1;
% func = @(x) sin(x*2*pi);

x_para = linspace(0,1,100);
x_phys = x_para*(b-a)+a;

knotVec = makeUniformKnotVector(p,el);
greville_para = findGrevillePoints(knotVec,p);
greville_phys = greville_para*(b-a) + a;

N = BsplineBasis(knotVec,p,x_para);
y_phys = func(x_phys);

d = N'\y_phys';


figure(1)
hold on 

% The points 
plot(greville_phys,d,'-k*')

% The exact function
plot(x_phys,func(x_phys),'-b', 'linewidth',2)



% The interpolated solution
plot(x_phys,N'*d,'r-',  'linewidth',2)
% Verification
norm(func(x_phys)'-N'*d)


%%%%%%%%%%%%%%%%%%%%%%%%%%%% TEst av andre punkter %%%%%%%%%%%%%%%%%%%
x_para2 = [0,0.1,0.3,0.35,0.39,0.4,0.7,0.8,0.99,1];
x_phys2 = x_para2*(b-a)+a;

N = BsplineBasis(knotVec,p,x_para2);

norm(func(x_phys2)'-N'*d)


figure(2)
hold on 

% The exact function
plot(x_phys,func(x_phys),'--k')
plot(x_phys2,func(x_phys2),'*m')

% The points 
%plot(greville_phys,d,'-*')

% The interpolated solution
plot(x_phys2,N'*d,'r-')
% Verification
norm(func(x_phys2)'-N'*d)