function knotVec = makeUniformKnotVector(p,el)

knotVec = [zeros(1,p), 0:1/el:1, ones(1,p)];

end