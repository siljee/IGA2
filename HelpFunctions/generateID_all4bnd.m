function ID = generateID_all4bnd(np_xi, np_eta, leftDirichlet, rightDirichlet, buttomDirichlet, topDirichlet)

Dbnd = 0;                   % Number of dirichlet bnd (left, right or both)

if ~isempty(leftDirichlet)
   Dbnd = Dbnd +1;
end
if ~isempty(rightDirichlet)
   Dbnd = Dbnd +1;
end
start = 0;


ID = [];
for i = 1:np_eta
    if (i==1 && ~isempty(buttomDirichlet))
        % Buttom Dirichlet
        if isscalar(buttomDirichlet)
            ID = buttomDirichlet*ones(1,np_xi);
        else
            ID = buttomDirichlet;
        end
        start = np_xi-Dbnd;
    elseif (i==np_eta && ~isempty(topDirichlet))
        % Top Dirichlet
        if isscalar(topDirichlet)
            ID = [ID, topDirichlet * ones(1,np_xi)];
        else
            ID = [ID, topDirichlet];
        end
    else
        % All other values
        leftValue = leftDirichlet;
        if length(leftDirichlet)>1
            leftValue = leftDirichlet(i);
        end
        
        rightValue = rightDirichlet;
        if length(rightDirichlet)>1
            rightValue = rightDirichlet(i);
        end
        ID = [ID, leftValue,(i-1)*(np_xi-Dbnd)+1-start:i*(np_xi-Dbnd)-start,rightValue];
    end
end