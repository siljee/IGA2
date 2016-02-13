function knotVec = makeRandomNonUniformKnotVector(elements, p, maxDist, isMult)

% Open knotVector to the left
knotVec = [zeros(1,p+1)];
value = 0;

for i = 1:elements-1
    
    value = value + ceil(rand()*maxDist);
    if isMult
        mult = floor(exp(rand()*log(p+1)));
    else
        mult = 1;
    end
    knotVec = [knotVec, ones(1,mult)*value];
end

% Open knotVector to the right
value = value + ceil(rand()*maxDist);
knotVec = [knotVec, ones(1,p+1)*value]/value;

end
