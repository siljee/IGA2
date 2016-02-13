function h = findMaxStepSize(knotVec)
    h = max(knotVec(2:end)-knotVec(1:end-1))/2;
end