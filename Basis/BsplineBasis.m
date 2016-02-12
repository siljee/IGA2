function N = BsplineBasis(knotVec, p, x)
% Constructs n univariate B-spline basis functions utilizing Cox-de Boor
% recursion formula.
%
% Input:
%    knotVec  - The vector of knots. i = 1,2,...,n+p+1. Knot values can
%               can be repeated.
%    p        - The polynomial order of the n basis functions.
%               (p = 0: constant, p=1: linear, p=2: quadratic etc.)
%    x        - A vector of values in parameter space where the basis
%               functions should be evaluated.
%
% Output:
%    N        - Array of B-spline basis functions. Each tow contain one
%               basis function evaluated on x.
%

n = length(knotVec)-p-1;                    % Number of bases
N = zerothPolynomialOrder(knotVec, n, x);

for order = (1:p)
    for basis = (1: n)
        N_1 = N(basis, :);
        if (basis == n)
            N_2 = 0;
        else
            N_2 = N(basis+1, :);
        end
        
        % Cox-de Boor recursion formula
        N(basis,:) = sparse( ...
            definedFraction(x, knotVec(basis), knotVec(basis+order), knotVec(basis)) .* N_1 ...
            + definedFraction(knotVec(basis+order+1), x, knotVec(basis+order+1), knotVec(basis+1)) .* N_2);
    end
end
end

function N = zerothPolynomialOrder(knotVec, n, x)
% Construct n bases N for the innermost recursion in Cox-de Boor recursion 
% formula, that is when p=0 and the basis function i is 0 everywhere except 
% in [knotVec(i), knotVec(i+1)] where it is constant 1. The basis functions 
% are defined on x in parameter space.

N = zeros(n, length(x));

for i = (1:n-1)
    N(i,:) = 0                    ...
        + (x >= knotVec(i))   .*1 ... 
        - (x >= knotVec(i+1)) .*1;
end
% The last basis should be defined on the last point. 
i = n;
N(n,:) = 0                    ...
    + (x >= knotVec(n))   .*1 ...
    - (x > knotVec(n+1)) .*1;

N = sparse(N);
end

function fraction = definedFraction(num1, num2, den3, den4)
% To make sure all fractions of Cox-de Boor formulae are defined the 
% assumption    0/0 = 0   must be explicitly implemented. 

numerator = num1 - num2;
denominator = den3 - den4;

fraction = numerator/denominator;
if (denominator == 0)
    fraction = 0;
end
end

