function elements = countElements(knotVec)
% This function counts how many nontrivial elements a knot vector represents,
% that is all elements with positive measure. This is done by counting
% the number of following knots that are different. 

elements = sum(knotVec(1:end-1)-knotVec(2:end)~= 0);
end