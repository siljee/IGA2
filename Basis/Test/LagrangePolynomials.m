function L = LagrangePolynomials(points, t, x)
% Generate Lagrange polynomials for the given p+1 points and parameter
% values t. Each polynomial will be of degree p and is evaluated on the 
% vector x.

p = length(points)-1;
n = length(x);

L = ones(p+1,n);
for i = 1:p+1           % Bases
    for j = 1:p+1       % Parameters
        if j~= i
            L(i,:) = L(i,:).*(x-t(j))/(t(i)-t(j));
        end
    end
end
end