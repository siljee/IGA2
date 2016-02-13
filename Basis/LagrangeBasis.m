function L = LagrangeBasis(p, t, x)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Constructs p+1 Lagrange polynomials L, defined by t, evaluated at x. 
%
% Input:
%    p        - The polynomial order of the p+1 basis functions.
%               (p = 0: constant, p=1: linear, p=2: quadratic etc.)
%    t        - Parametric values. 
%    x        - A vector of evaluation points.
%
% Output:
%    L        - Array ofLagrange polynomials. One basis in each row,
%               defined on points x, one in each column. Each basis
%               have value one at one t and zero at all other t's.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

L = ones(p+1,length(x));

for i = 1:p+1           % Bases
    for j = 1:p+1       % Parameters
        if j~= i
            L(i,:) = L(i,:).*(x-t(j))/(t(i)-t(j));
        end
    end
end