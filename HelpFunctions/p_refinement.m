function [Px, Py] = p_refinement(knotVec_xi, knotVec_eta, p_xi, p_eta, Px, Py)
newKnots_xi = findElements(knotVec_xi);
newKnots_xi = knotVec_xi(newKnots_xi(2:end));
if isempty(newKnots_xi)
    newKnots_xi = (knotVec_xi(end)-knotVec_xi(1))/2;
end

newKnots_eta = findElements(knotVec_eta);
newKnots_eta = knotVec_eta(newKnots_eta(2:end));
if isempty(newKnots_eta)
    newKnots_eta = (knotVec_eta(end)-knotVec_eta(1))/2;
end 

[Px, Py, knotVec_xi] = knotInsertion(knotVec_xi, p_xi, newKnots_xi ,Px,Py);
[Py, Px, knotVec_eta] = knotInsertion(knotVec_eta, p_eta, newKnots_eta,Py',Px');
Px = Px';
Py = Py';

end