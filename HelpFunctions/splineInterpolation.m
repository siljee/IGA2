function g = splineInterpolation(U_exact, knotVec, p, P)

greville_para = findGrevillePoints(knotVec,p);
x_para = greville_para; %linspace(0,1,10);

N = BsplineBasis(knotVec,p,x_para);
x_phys = N'*P;
y_phys = U_exact(x_phys);

g = N'\y_phys;

end