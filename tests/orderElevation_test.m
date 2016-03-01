addpath ../HelpFunctions/

p_xi = 2; p_eta = 2;
el_xi = 2; el_eta = 1;
knotVec_xi  = makeUniformKnotVector(p_xi, el_xi);
knotVec_eta = makeUniformKnotVector(p_eta,el_eta);

np_xi  = length(knotVec_xi)  - p_xi  - 1;                 % Number of control points in xi direction
np_eta = length(knotVec_eta) - p_eta - 1;                 % Number of control points in eta direction
 

Ax = 0; Bx = 1;
Ay = 0; By = 1;


Px = [0,1,2; 0,1,2; 0, 9/4, 2.5; 3,3,3]
Py = [0,0,0; 3, 3/4, 0.5; 3,2,1; 3,2,1]

Px = (findGrevillePoints(knotVec_xi, p_xi)-knotVec_xi(1))* (Bx-Ax)/(knotVec_xi(end)-knotVec_xi(1)) + Ax
Px = repmat(Px',1,np_eta);
Py = (findGrevillePoints(knotVec_eta, p_eta)-knotVec_eta(1))*(By-Ay)/(knotVec_eta(end)-knotVec_eta(1)) + Ay; 
Py = repmat(Py,np_xi,1);


xN = linspace(0,1,10);

N_xi = BsplineBasis(knotVec_xi, p_xi, xN);
N_eta = BsplineBasis(knotVec_eta, p_eta, xN);

Cx = N_xi'*Px*N_eta;
Cy = N_xi'*Py*N_eta;

close all
figure(1)
hold on
surf(Cx,Cy, zeros(size(Cx)),'linestyle', 'none')
plot(Px,Py,'*-k')
plot(Px',Py','*-k')

% h_refinement(knotVec_xi, knotVec_eta,p_xi,p_eta,Px,Py)
[Px, Py] = p_refinement(knotVec_xi, knotVec_eta, p_xi, p_eta, Px, Py)   

figure(2)
hold on
surf(Cx,Cy, zeros(size(Cx)),'linestyle', 'none')
plot(Px,Py,'*-k')
plot(Px',Py','*-k')

%%%% Use these for p=3
p_xi = 3; p_eta = 3;
el_xi = 2; el_eta = 1;
knotVec_xi  = makeUniformKnotVector(p_xi, el_xi);
knotVec_eta = makeUniformKnotVector(p_eta,el_eta);

N_xi = BsplineBasis(knotVec_xi, p_xi, xN);
N_eta = BsplineBasis(knotVec_eta, p_eta, xN);

Cx = N_xi'*Px*N_eta;
Cy = N_xi'*Py*N_eta;

figure(3)
hold on
surf(Cx,Cy, zeros(size(Cx)),'linestyle', 'none')
plot(Px,Py,'*-k')
plot(Px',Py','*-k')
