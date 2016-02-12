function B = BernsteinPolynomial(p,x)
% Generate p+1 Bernstein polynomials B, evaluated at x. Values of the 
% vector x must be in [0,1]. The p+1 Bernstein polynomials can be 
% used to form Bezier curves of degree p.

n = length(x);
B = zeros(p+1,n);

for v = 0:p 
    B(v+1,:) = nchoosek(p,v) .* x.^v .* (1-x).^(p-v);
end
end