function B = BernsteinBasis(p,x)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Constructs p+1 Bernstein polynomials.
%
% Input:
%    p        - The polynomial order of the p+1 basis functions.
%               (p = 0: constant, p=1: linear, p=2: quadratic etc.)
%    x        - A vector of evaluation points in Bernstein parameter space.
%
% Output:
%    B        - Array of Bernstein polynomials. One basis in each row,
%               defined on points x, one in each column.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

n = length(x);
B = zeros(p+1,n);

for v = 0:p 
    B(v+1,:) = nchoosek(p,v) .* x.^v .* (1-x).^(p-v);
end
end
