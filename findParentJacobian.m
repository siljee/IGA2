function [A, A_xi, A_eta] = findParentJacobian( knotVec_xi, knotVec_eta, j, i)

% Needed when knot vectors are nonuniform.
knot_j = findElements(knotVec_xi);
knot_i = findElements(knotVec_eta);
knot_j = knot_j(j);
knot_i = knot_i(i);

% Jacobian is the area of parameter, since parent domain is unit square.
A_xi = ((knotVec_xi(knot_j+1)- knotVec_xi(knot_j)));
A_eta = ((knotVec_eta(knot_i+1)- knotVec_eta(knot_i)));
A = A_xi*A_eta;
end