function [P_dir, P_otherDir, knotVec_dir] = knotInsertion2D...
    (knotVec_dir, p_dir, newKnots_dir, P_dir, P_otherDir)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% For a two-dimensional object, new knots are inserted in one 
% parametric direction and control points are restructured to preserve 
% the objects parametric and geometric properties. 
%
% Input:
%   knotVect_dir  Original knot vector for the parametric direction
%                 of inserted knots.
%   p_dir         Polynomial order for the parametric direction of
%                 inserted knots.
%   newKnots_dir  Knots to be inserted into the knot vector for the 
%                 parametric direction of inserted knots. 
%   P_dir         Control points for one dimension in physical space.
%   P_otherDir    Control points for the other dimenstion in physical
%                 space.
%
% Otput:
%   P_dir         The new control points in one dimension.
%   P_otherDir    The new control points in the other dimension.
%   knotVec_dir   The new knotVec_dir after insertion of newKnots_dir. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


n = length(knotVec_dir)-p_dir-1;        % Number of bases 
m = length(newKnots_dir);               % Number of new knots

for j = 1:m                % index of inserted knot
    k = find(knotVec_dir <= newKnots_dir(j),1,'last');   
                           % index of knot left of inserted knot
                           
    % Generate alphas
    alpha = zeros(n+j-1,1);
    alpha(1:k-p_dir) = 1;
    for i = (k-p_dir+1:k)
        denominator = knotVec_dir(i+p_dir) - knotVec_dir(i);
        if denominator ~= 0
            alpha(i) = (newKnots_dir(j)-knotVec_dir(i)) / denominator;
        end
    end

    alpha = repmat(alpha(2:end),1,size(P_dir,2));
    
    % Update P and knot vector
    P_dir = [...
        P_dir(1,:); ...
        alpha.*P_dir(2:end,:) + (1-alpha).*P_dir(1:end-1,:); ...
        P_dir(end,:)];
    P_otherDir = [...
        P_otherDir(1,:); ...
        alpha.*P_otherDir(2:end,:)+(1-alpha).*P_otherDir(1:end-1,:);...
        P_otherDir(end,:)];
    
    knotVec_dir = [...
        knotVec_dir(1:k),...
        newKnots_dir(j), ...
        knotVec_dir(k+1:end)];
end

