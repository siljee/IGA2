function elements = findElements(knotVec)
% Find startindex of every element. This function is utilized when 
% dealing with a knotVector of repeated knot values.


elements = find(knotVec(1:end-1)-knotVec(2:end)~= 0);

end