function dB = BernsteinDerivative(p, x)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Derivated Bernstein polynomials of order p and return p+1 derivated 
% Bernstein polynomials of order p-1.
%
% Input:
%    p        - The polynomial order of the original p+1 basis functions.
%               (p = 0: constant, p=1: linear, p=2: quadratic etc.)
%    x        - A vector of evaluation points in Bernstein parameter space.
%
% Output:
%    dB       - Array of derivated Bernstein polynomials of order p-1. 
%               One basis in each row, defined on points x, one in each column.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

n = length(x);
dB = zeros(p+1,n);

%Generate the support, such that:
% B_(-1,p-1) = 0, B_(0,p-1), ... B_(p-1,p-1) , B_(p,p-1) = 0
B_supp = zeros(p+2,n);              
B_supp(2:p+1,:) = BernsteinBasis(p-1,x);

for v = 1:p+1
    dB(v,:) = p * (B_supp(v,:) - B_supp(v+1,:)); 
end
end
