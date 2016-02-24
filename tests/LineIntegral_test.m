addpath ../; addpath ../HelpFunctions/; addpath ../Basis/;

p = 2;
el = 3;
a = -1; b = a+el;


knotVec = makeUniformKnotVector(p,el);
[G1,W1] = GaussQuadrature(p+1);
[G2,W2] = GaussQuadrature(p+2);
h = @(x) 2*ones(size(x));

N_eta = BsplineBasis(knotVec,p,G2');
full(N_eta);
Py = findGrevillePoints(knotVec,p)*(b-a) + a
Py*N_eta
G2'*(b-a)+a

% "Element 1"


%N_eta_e = N_eta(1:p+1,:);
%full(N_eta_e);
%Py_e = Py(1:3);

%N_eta_e'*Py_e';
%G*(b/el-a)+a;
%Py_e*N_eta_e(:,1);



%h = @(x) zeros(size(x));

%I = h(Py(1))*N_eta_e(1,:) + h(Py(2))*N_eta_e(2,:) + h(Py(3))*N_eta_e(3,:);
%integral(h,0,1/el);



B_eta1 = BernsteinBasis(p,G1);
B_eta2 = BernsteinBasis(p,G2);
B_eta = BernsteinBasis(p,linspace(0,1,100))

close all
figure(2)
hold on
plot(linspace(0,1,100),B_eta',':b')
plot(G1,B_eta1','*-k')
plot(G2,B_eta2','*r--')

C = BezierExtractionOperator(knotVec, p);

% Element 1

Py_e = Py(1:p+1);

N_eta = C(:,:,1)*B_eta2
Py_e*N_eta
G2'*(b-a)/el+a
I = h(Py(1))*N_eta(:,1)*W2(1) + h(Py(2))*N_eta(:,2)*W2(2) + h(Py(3))*N_eta(:,3)*W2(3)
%I = 2*W2(1) + 2*W2(2) + 2*W2(3)
%I = (N_eta(:,1)*W(1)+N_eta(:,2)*W(2))*(0-a)
int = @(x) 3*x;
int_phys = @(x) int((x+1)/3)
integral(int_phys,a,a+(b-a)/el)
