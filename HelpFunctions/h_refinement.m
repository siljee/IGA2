function [knotVec_xi, knotVec_eta, Px, Py] = h_refinement(knotVec_xi, knotVec_eta ,p_xi, p_eta, Px, Py)

uniqueKnotVec_xi = [knotVec_xi(findElements(knotVec_xi)) knotVec_xi(end)];
uniqueKnotVec_eta = [knotVec_eta(findElements(knotVec_eta)) knotVec_eta(end)];
newKnots_xi = uniqueKnotVec_xi(1:end-1) + (uniqueKnotVec_xi(2:end) - uniqueKnotVec_xi(1:end-1))/2;
newKnots_eta = uniqueKnotVec_eta(1:end-1) + (uniqueKnotVec_eta(2:end) - uniqueKnotVec_eta(1:end-1))/2;

[Px, Py, knotVec_xi] = knotInsertion(knotVec_xi, p_xi, newKnots_xi, Px, Py);
[Py, Px, knotVec_eta] = knotInsertion(knotVec_eta, p_eta, newKnots_eta, Py', Px');
Py = Py';
Px = Px';
end