addpath ../; addpath ../Basis/; addpath ../HelpFunctions/; clear all
p = 2; el = 1; a= 1; b = 2; nx = 10000;
% knotVec = makeUniformKnotVector(p,el);
knotVec = makeRandomNonUniformKnotVector(el,p,2,1);
n = length(knotVec)-p-1;
x = linspace(a,b,nx);
xN = x/(b-a) -a;
Px = findGrevillePoints(knotVec,p)*(b-a) +a;
%Py = rand(n,1);
Py = [1,2,3]
Py = repmat(Py,2,1);
N = BsplineBasis(knotVec,p,xN);

[Py_new, knotVec_new] = h_refinement(knotVec,p,Py)

Px_new = findGrevillePoints(knotVec_new,p)*(b-a) +a;
N_new = BsplineBasis(knotVec_new,p,xN);

close all
figure(1)
hold on
plot(Px,Py,'*-')
plot(x,N'*Py, 'linewidth',2)

axis([a,b,0,1])


hold on
plot(x,N_new'*Py_new')
plot(Px_new,Py_new,'*-')

axis([a,b,0,1])
