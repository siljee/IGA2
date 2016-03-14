function U = IGA_Method_Bezier(knotVec_xi, knotVec_eta, p_xi, p_eta, el_xi, el_eta, np_xi, np_eta, Px, Py, f, DIR, ID, IEN, LM, h, NEU, Gp_xi, Gp_eta) 

% Initialize 
dof = sum(ID>0);            % Degrees of freedom
A = sparse(dof, dof);       % Stiffness matrix
b = sparse(dof,1);          % Load vector

% Bezier extraction operator
C_xi  = BezierExtractionOperator(knotVec_xi ,p_xi);
C_eta = BezierExtractionOperator(knotVec_eta,p_eta);

% Berstein basis and derivative need only be evaluated once.
[B_xi, B_eta, dB_xi, dB_eta, W, W_xi, W_eta] = ...
    bernsteinBasisAndDerivative(p_xi, p_eta, Gp_xi, Gp_eta);

for e_eta = 1:el_eta
  for e_xi = 1:el_xi
    e = (e_eta-1)*el_xi + e_xi;
    l = 0;
        
    % Element variables
    LM_e = LM(e,:); IEN_e = IEN(e,:);      
    Px_e = Px(IEN_e)'; Py_e = Py(IEN_e)';
    innerN = LM_e(LM_e > 0);
        
    % Basis functions and Jacobian. 
    [Jpar, Jpar_xi, Jpar_eta] = findParentJacobian(...
        knotVec_xi, knotVec_eta, e_xi, e_eta);
    [N, dNdxi, dNdeta, N_xi, N_eta] = shapeFunctions(B_xi, B_eta,...
        dB_xi, dB_eta, C_xi(:,:,e_xi), C_eta(:,:,e_eta),Jpar_xi, Jpar_eta); 
    [Jphys, dNdx, dNdy, dxdxi, dydeta] = Jacobian(dNdxi, dNdeta,...
        Px_e, Py_e, Gp_xi*Gp_eta);
    J = Jpar*Jphys;
      
    % Assemble load vector and stiffness matrix 
    l = l + localLoadVector(f, N, Px_e, Py_e, LM_e, Gp_xi*Gp_eta,...
        W, J);
    l = l + localNeumannBnd(NEU, h, Px, Py, N_xi, N_eta, dxdxi, dydeta,... 
        Jpar_xi, Jpar_eta, LM_e, IEN_e, Gp_xi, Gp_eta, W_xi, W_eta);
    l = l + localDirichletBnd(DIR, dNdx, dNdy, Gp_xi*Gp_eta, W, ...
        J, LM_e);
    
    b(innerN) = b(innerN) + l;
    A(innerN,innerN) = A(innerN,innerN) + localStiffnessMatrix(...
        dNdx(LM_e>0,:), dNdy(LM_e>0,:), Gp_xi*Gp_eta, W, J);
  end
end

% Solve system
U = A\b;
U = addDirichelBnd(U, ID, DIR);
end


