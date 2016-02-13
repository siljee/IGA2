function N = BsplineBasis(knotVec, p, x)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Constructs n univariate B-spline basis functions utilizing
% Cox-de Boor recursion formula.
%
% Input:
%    knotVec  - The vector of knots. i = 1,2,...,n+p+1. Knot values 
%               can be repeated.
%    p        - The polynomial order of the n basis functions.
%               (p = 0: constant, p=1: linear, p=2: quadratic etc.)
%    x        - A vector of evaluation points in parameter space.
%
% Output:
%    N        - Array of B-spline basis functions. One basis in each 
%               row, defined on points x, one in each column.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

n = length(knotVec)-p-1;        % Number of bases
N = zeros(n, length(x));        % Initialize basis array

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initializing constant basis functions for p=0 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i = (1:n-1)
    N(i,:) = 0                    ...
        + (x >= knotVec(i))   .*1 ...
        - (x >= knotVec(i+1)) .*1;
end

% The last basis should be defined on the last point.
N(n,:) = 0                    ...
    + (x >= knotVec(n))  .*1  ...
    - (x > knotVec(n+1)) .*1;
N = sparse(N);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Generate basis functions N for p > 0
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
              definedFraction(x, knotVec(basis), ...
                    knotVec(basis+order), knotVec(basis)) .* N_1 ...
            + definedFraction(knotVec(basis+order+1), x, ...
                    knotVec(basis+order+1), knotVec(basis+1)) .* N_2);
    end
end
end


function fraction = definedFraction(num1, num2, den3, den4)
% Make defined fraction based on the condition:  0/0 = 0.

denominator = den3 - den4;
if (denominator == 0)
    fraction = 0; return;
end
fraction = (num1-num2)/denominator;
end
