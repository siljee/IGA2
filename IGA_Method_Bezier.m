function U = IGA_Method_Bezier(knotVec_xi, knotVec_eta, p_xi, p_eta, el_xi, el_eta, np_xi, np_eta, Px, Py, f,  dirichletValues, ID, IEN, LM) %leftVal, rightVal)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%       Initialisere viktige variabler
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%np_xi  = length(knotVec_xi)  - p_xi  - 1;                % Number of control points in xi direction
%np_eta = length(knotVec_eta) - p_eta - 1;                % Number of control points in eta direction

C_xi  = BezierExtractionOperator(knotVec_xi ,p_xi);
C_eta = BezierExtractionOperator(knotVec_eta,p_eta);

%IEN = generate_IEN(knotVec_xi, knotVec_eta, p_xi, p_eta);
%ID  = generateID_all4bnd(np_xi, np_eta, -1, -2,[],[]);
%LM  = ID(IEN);
 
dof = sum(ID>0);                      % Degrees of freedom
Gp_xi = p_xi+1; Gp_eta = p_eta +1;    % Number of Gauss points
Gp = Gp_xi * Gp_eta;

% Initialize
K = sparse(dof, dof);
R = sparse(dof,1);

[G_xi_tilde ,W_xi] = GaussQuadrature(Gp_xi);
[G_eta_tilde ,W_eta] = GaussQuadrature(Gp_eta);
W = kron(W_eta,W_xi);

elementStartKnots_xi = findElements(knotVec_xi);
elementStartKnots_eta = findElements(knotVec_eta);

for i = 1:el_eta
    for j = 1:el_xi
        e = (i-1)*el_xi +j;
        k = 0; l = 0;
        
        knot_j = elementStartKnots_xi(j);
        knot_i = elementStartKnots_eta(i);
       
        A_xi = (knotVec_xi(knot_j+1)- knotVec_xi(knot_j));
        A_eta = (knotVec_eta(knot_i+1)- knotVec_eta(knot_i));
        A = A_xi*A_eta;
        Jparent = A;
        %We = 1;
        %We = ((knotVec_xi(p_xi+j+1)- knotVec_xi(p_xi+j))*(knotVec_eta(p_eta+i+1)- knotVec_eta(p_eta+i)))
        
        %G_xi = G_xi_tilde*A_xi; G_eta = G_eta_tilde*A_eta; % Map gauss points from parent space to parameter space
        
        Ce = kron(C_eta(:,:,i),C_xi(:,:,j));        % (p+1)*(q+1) x (p+1)*(q+1)
        IEN_e = reshape(IEN(e,:),p_xi+1,p_eta+1);   % [IEN(e,:), IEN(e,:)+ncp];
        LM_e = LM(e,:);
        inner = LM_e(LM_e > 0);
        dof = sum(ID>0); 

        
        [N, dNdxi, dNdeta, N_xi, N_eta] = shapeFunctions(p_xi, p_eta, G_xi_tilde, G_eta_tilde, Ce, C_xi(:,:,j), C_eta(:,:,i), A_xi, A_eta); % Px(IEN_e)); %, Py(IEN_e));
        [detJ, dNdx, dNdy] = Jacobian(dNdxi, dNdeta, Px(IEN_e(:)), Py(IEN_e(:)), Gp);
        
        
        % PARAMETER SPACE
        %%%%%%%% Funksjonene er nå i parameter space. De er evaluert for
        %%%%%%%% gausspunkt fra parentspace som er mappet over til
        %%%%%%%% parameter space. Altså trengs ingen mer mapping annet enn
        %%%%%%%% over til physical space nå
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Stiffness matrix
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %gradPhi2g = dNdx(LM_e>0,:)*dNdx(LM_e>0,:)' + dNdy(LM_e>0,:)*dNdy(LM_e>0,:)'
        for g = 1:Gp
            gradPhi2g = dNdx(LM_e>0,g)*dNdx(LM_e>0,g)' + dNdy(LM_e>0,g)*dNdy(LM_e>0,g)';
            k = k + gradPhi2g*W(g)*Jparent*detJ(g);
        end

        K(inner,inner) = K(inner,inner) + k;
        % K er evaluert i parameter space! Nei. De deriverte er med hensyn
        % på x og er vel evaluert i physical space?
       
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Load vector l
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %f = repmat((pi^2*sin(pi*G_xi))',Gp_eta,1); %,1)
%         j
%         A_xi
%         A_xi*G_xi_tilde
%         'gauss punkter i parameter: '
%         A_xi*G_xi_tilde+knotVec_xi(p_xi+j)
%         Px(IEN_e(:,:))'
%         (A_xi*G_xi_tilde+knotVec_xi(p_xi+j))'*Px(IEN_e(:,:))'  % Gauss punkter i parameter space
%         'c'
% %         G_xi_tilde\C_xi(:,:,j)'
% C_xi(:,:,j)'
        %C_xi(:,:,j)'*Px(IEN_e(:,1))
%         C_xi(:,:,j)= Px(IEN_e)
        %knotVec_xi(p_xi+j)
        %knotVec_xi(p_xi+j+1)
        %Px(IEN_e(:,1))
        % N'*Px    daah!    x_physical = 1*(repmat((G_xi_tilde*A_xi + knotVec_xi(p_xi+j)),Gp_eta,1))
        %(G_xi_tilde*A_xi + knotVec_xi(p_xi+j)),Gp_eta,1)'*Px(IEN_e(:,:))
        
        % Kanskje f heller skal evalueres i de (dobbelt)mappede gausspunktene!! 
        % FEILLLLL  fx = N'*f(Px(IEN_e(:)), Py(IEN_e(:))) % Det kan være Py skal transponeres!
        %Px
        %Px(IEN_e)
        %Px(IEN_e(:,1))
%         N'*Px(IEN_e(:));
%         N_xi'*Px(IEN_e(:,1));

        fx = f(N_xi'*Px(IEN_e)*N_eta,N'*Py(IEN_e(:)));
        %fx = f(x_physical,1) % Disse to er de samme

        % Men jeg lurer på om denne N kanskje er evaluert i parameter
        % space?
        for g = 1:Gp
            %N(LM_e>0,g)'* Px(LM_e>0,:)
            l = l + fx(g)*N(LM_e>0,g)*W(g)*Jparent*detJ(g);
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Nonhomo Dirichlet
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % rightside
        
        % Ta en for her
        for dVi = 1:length(dirichletValues)
            %         dVi = 1;   % DirichletValueIndex
            if ismember(-dVi, LM_e) %&& rightVal ~= 0
                nDbnd = find(LM_e == -dVi);
                a = 0;
                
                for g = 1:Gp
                    gradPhi2g = dNdx(:,g)*dNdx(:,g)' + dNdy(:,g)*dNdy(:,g)';
                    a = a + gradPhi2g*W(g)*Jparent*detJ(g);
                end
                l = l - dirichletValues(dVi)*sum(a(LM_e>0, nDbnd),2);
            end
        end
        R(inner) = R(inner) + l;
        %disp('-------------------')
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Sove and add Dirichlet boundary
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        U = K\R;
        U = addDirichelBnd(U, ID, dirichletValues);
end


