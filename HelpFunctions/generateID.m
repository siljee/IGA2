function ID = generateID(np_xi, np_eta, leftDirichlet, rightDirichlet, buttomDirichlet, topDirichlet)

Dbnd = 0;                   % Number of dirichlet bnd (left, right or both)

if ~isempty(leftDirichlet)
   Dbnd = Dbnd +1;
end
if ~isempty(rightDirichlet)
   Dbnd = Dbnd +1;
end
IDindex = 1;

ID = [];
for i = 1:np_eta
    if (i==1)
        % Buttom Dirichlet
        if isempty(buttomDirichlet)
            % If not a dirichlet bnd
            ID = [ID, IDindex:np_xi];
            IDindex = IDindex + np_xi;
        else            
            if isempty(leftDirichlet)
                % If left neuanann boundary, include the positive equation index
                ID = [ID, IDindex];
                IDindex = IDindex + 1;
                buttomDirichlet = buttomDirichlet(2:end);
            end
            ID = [ID, buttomDirichlet];
            if isempty(rightDirichlet)
                % If right neuanann boundary, include the positive equation index
                ID(end) = IDindex;
                IDindex = IDindex + 1;
            end
        end
    elseif (i==np_eta)
        if isempty(topDirichlet)
            ID = [ID, IDindex:IDindex+np_xi-1];
            IDindex = IDindex+np_xi;
        else
            if isempty(leftDirichlet)
                % If left neuanann boundary, include the positive equation index
                ID = [ID, IDindex];
                IDindex = IDindex + 1;
                topDirichlet = topDirichlet(2:end);
            end
            ID = [ID, topDirichlet];
            if isempty(rightDirichlet)
                % If right neuanann boundary, include the positive equation index
                ID(end) = IDindex;
                IDindex = IDindex + 1;
            end
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
        ID = [ID, leftValue, IDindex:(IDindex-1 + np_xi-Dbnd), rightValue];
        IDindex = IDindex +np_xi-Dbnd;
    end
end