% This function compute local Bezier extraction operators for all 
% elements. The extraction operator can then be used to map univariate 
% Bernstein polynomials to univariate B-spline basis functions.
%
% Input
%   knotVec     The knotVector used to define the univariate B-spline.
%   p           Polynomial order of B-spline basis function.
%
% Output
%   C           A 3-dimensional matrix with one Bezier extraction 
%               operator matrix for each element.
%
function C =  BezierExtractionOperator(knotVec, p)
m = length(knotVec);
a = p+1;  mult = 1;  el = 1;

C = zeros(p+1,p+1,1);
I = eye(p+1,p+1);
C(:,:,el) = I;

while (a+mult < m)
    % Count multiplicity
    while a+mult<m && (knotVec(a+mult+1) == knotVec(a+mult))
        mult = mult+1;
    end
    
    if mult <= p
        C(:,:,el+1) = I;
        numerator = knotVec(a+mult) - knotVec(a);
        % Generate alphas
        for j = (p:-1:mult+1)
            alphas(j-mult) = numerator/ (knotVec(a+j) - knotVec(a));
        end
        nnK = p - mult;   % Number of new knots.
        
        % Update matrix cofficients for new knots
        for j = 1:nnK
            save = nnK-j+1;
            s = mult + j;
            for k = p+1:-1:s+1
                alpha = alphas(k-s);
                C(:,k,el) = alpha*C(:,k,el) + (1-alpha)*C(:,k-1,el);
            end
            if a+mult < m
                % Update overlapping coefficients for next operator
                C(save:j+save, save, el+1) = C(p-j+1:p+1,p+1, el);
            end
        end
        el = el+1;
    end

    if a + mult < m
        % update idices for next element
        a = a + mult;
        mult = 1;
    end
end
end
