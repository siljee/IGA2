function [P_dir, P_otherDir, knotVec_dir] = knotInsertion(knotVec_dir, p_dir, newKnots_dir, P_dir, P_otherDir)
% Input:
% knotVector    The original knot vector 
% newKnots      The new knots we want to append
% p             The order of the basis functions
%
% Returns:
% C             The Bezier Extraction operator
% knotVector    The new knotVector, updated with newKnots. 
%



n = length(knotVec_dir)-p_dir-1;                 % Number of bases 
m = length(newKnots_dir);                    % Number of new knots
%C_transpose = 1;   


for j = 1:m                                     % index of inserted knot
    k = find(knotVec_dir < newKnots_dir(j),1,'last');   % index of knot left of inserted knot
   
    % Generate alphas
    alpha = zeros(n+j-1,1);
    alpha(1:k-p_dir) = 1;
    for i = (k-p_dir+1:k)
        denominator = knotVec_dir(i+p_dir) - knotVec_dir(i);
        if denominator ~= 0
            alpha(i) = (newKnots_dir(j) - knotVec_dir(i)) / denominator;
        end
    end

    alpha = repmat(alpha(2:end),1,size(P_dir,2));
    
    % Update
    P_dir = [P_dir(1,:); alpha.*P_dir(2:end,:) + (1-alpha).*P_dir(1:end-1,:); P_dir(end,:)];
    P_otherDir = [P_otherDir(1,:); alpha.*P_otherDir(2:end,:) + (1-alpha).*P_otherDir(1:end-1,:); P_otherDir(end,:)];
    
    knotVec_dir = [knotVec_dir(1:k), newKnots_dir(j), knotVec_dir(k+1:end)];
 
end

