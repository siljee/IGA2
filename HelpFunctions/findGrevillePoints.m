function greville = findGrevillePoints(knotVec,p)
% This function finds the Greville points based on a knot vector and 
% the polynomial order of the spline. 

greville = 0;
for i = 1:p
    greville = greville + knotVec(1+i:end-(p+1-i));
end
greville = greville/p;
end