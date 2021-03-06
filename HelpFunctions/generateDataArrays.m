function [ID, IEN, LM, DIR, NEU] = generateDataArrays(knotVec_xi, knotVec_eta, p_xi, p_eta, bnd, u_exact, domain, Px, Py)

DIR = [];

np_xi = length(knotVec_xi) - p_xi -1;                   % Number of control points in xi direction
np_eta = length(knotVec_eta) - p_eta -1;                % Number of control points in eta direction



greville_xi  = findGrevillePoints(knotVec_xi, p_xi);
greville_eta = findGrevillePoints(knotVec_eta, p_eta);




N_xi = BsplineBasis(knotVec_xi, p_xi, greville_xi);
N_eta = BsplineBasis(knotVec_eta, p_eta, greville_eta);


greville_xi_mat  = repmat(Px(:,1)'*N_xi,length(greville_eta),1);
greville_eta_mat = repmat(N_eta*Py(1,:)',1,length(greville_xi));

grevilleValues = N_eta'*u_exact(greville_xi_mat,greville_eta_mat)*N_xi;


dirCount = 1;

% Go through all boundaries (left, right, top, buttom)
% and add LM values <= 0 for all bnd with type 'd'.
%
% When the dirichlet bnd is homogeneous (value = 'h') the
% LM index is 0.
% When the dirichlet bnd is constant (value = 'c') the
% LM index is one negative number.
% When the dirichlet bnd is variable (value = 'v') the
% LM index multiple negative numbers.

bndNames = fieldnames(bnd);
for i = 1:numel(bndNames)
    % Dirichlet Boundary
    % Generate DIR and make vectors necessary for ID array.
    if bnd.(bndNames{i}).type == 'd'
        switch bnd.(bndNames{i}).value
            case 'h'
                switch cell2mat(bndNames(i))
                    case 'left'
                        dirID.(bndNames{i}) = zeros(1,np_eta);
                    case 'right'
                        dirID.(bndNames{i}) = zeros(1,np_eta);
                    case 'buttom'
                        dirID.(bndNames{i}) = zeros(1,np_xi);
                    case 'top'
                        dirID.(bndNames{i}) = zeros(1,np_xi);
                end
            case 'c'
                switch cell2mat(bndNames(i))
                    case 'left'
                        % independent on y
                        DIR = [DIR; u_exact(domain.startX, domain.startY)];
                        dirID.(bndNames{i}) = -dirCount*ones(1,np_eta);
                    case 'right'
                        % independent on y
                        DIR = [DIR; u_exact(domain.endX, domain.startY)];
                        dirID.(bndNames{i}) = -dirCount*ones(1,np_eta);
                    case 'buttom'
                        % independent on x
                        DIR = [DIR; u_exact(domain.startX, domain.startY)];
                        dirID.(bndNames{i}) = -dirCount*ones(1,np_xi);
                    case 'top'
                        % independent on x
                        DIR = [DIR; u_exact(domain.startX, domain.endY)];      
                        dirID.(bndNames{i}) = -dirCount*ones(1,np_xi);
                end
                dirCount = dirCount + 1;
                
            case 'v'
                switch cell2mat(bndNames(i))
                    case 'left'
                        dirID.(bndNames{i}) = -dirCount:-1:-np_eta-dirCount+1;
                        dirCount = dirCount + np_eta;
                        %DIR = [DIR; grevilleValues(:,1)];
                        U_exact_left = @(y) u_exact(Px(1,:)',y);
                        DIR = [DIR; splineInterpolation(U_exact_left, knotVec_eta, p_eta, Py(1,:)')];
                    case 'right'
                        dirID.(bndNames{i}) = -dirCount:-1:-np_eta-dirCount+1;
                        dirCount = dirCount + np_eta;
                        %DIR = [DIR; grevilleValues(:,end)];
                        U_exact_right = @(y) u_exact(Px(end,:)',y);
                        DIR = [DIR; splineInterpolation(U_exact_right, knotVec_eta, p_eta, Py(end,:)')];
                    case 'buttom'
                        %if any(NEU2 == 2)
                            dirID.(bndNames{i}) = -dirCount:-1:-np_xi-dirCount+1;
                        %end
                        dirCount = dirCount + np_xi;
                        %DIR = [DIR; grevilleValues(1,:)'];
                        U_exact_buttom = @(x) u_exact(x,Py(:,1));
                        DIR = [DIR; splineInterpolation(U_exact_buttom, knotVec_xi, p_xi, Px(:,1))];
                    case 'top'
                        dirID.(bndNames{i}) = -dirCount:-1:-np_xi-dirCount+1;
                        dirCount = dirCount + np_xi;
                        %DIR = [DIR; grevilleValues(end,:)'];
                        U_exact_top = @(x) u_exact(x,Py(:,end));
                        DIR = [DIR; splineInterpolation(U_exact_top, knotVec_xi, p_xi, Px(:,end))];
                end
            otherwise
                error('Value of bnd with type d are wrongly specified! Use h, c or v')
        end
    else
        dirID.(bndNames{i}) = [];
    end
end

ID = generateID_all4bnd(np_xi, np_eta, dirID.left ,dirID.right, dirID.buttom, dirID.top);

IEN = generate_IEN(knotVec_xi, knotVec_eta, p_xi, p_eta);
LM = ID(IEN);
ID2 = reshape(ID,np_xi,np_eta)';

NEU = [];
% Neumann Boundary
if bnd.left.type == 'n'
    value = bnd.left.value;
    if value == 'c' || value == 'v';
        NEU = [NEU; ID2(ID2(:,1)>0,1), ones(sum(ID2(:,1)>0),1)];
    end
end
if bnd.right.type == 'n'
    value = bnd.right.value;
    if value == 'c' || value == 'v';
        NEU = [NEU; ID2(ID2(:,end)>0,end), 2*ones(sum(ID2(:,end)>0),1)];
    end
end
if bnd.buttom.type == 'n'
    value = bnd.buttom.value;
    if value == 'c' || value == 'v';
        NEU = [NEU; ID2(1,ID2(1,:)>0)', 3*ones(sum(ID2(1,:)>0),1)];
    end
end
if bnd.top.type == 'n'
    value = bnd.top.value;
    if value == 'c' || value == 'v';
        NEU = [NEU; ID2(end,ID2(end,:)>0)', 4*ones(sum(ID2(end,:)>0),1)];
    end
end

