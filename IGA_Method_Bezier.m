function U = IGA_Method_Bezier(knotVec_xi, knotVec_eta, p_xi, p_eta, el_xi, el_eta, np_xi, np_eta, Px, Py, f, DIR, ID, IEN, LM, h_neu, NEU) %leftVal, rightVal)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%       Initialisere viktige variabler
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

C_xi  = BezierExtractionOperator(knotVec_xi ,p_xi);
C_eta = BezierExtractionOperator(knotVec_eta,p_eta);

dof = sum(ID>0);                      % Degrees of freedom
Gp_xi = p_xi+1; Gp_eta = p_eta +1;    % Number of Gauss points
Gp = Gp_xi * Gp_eta;

% Initialize
K = sparse(dof, dof);
R = sparse(dof,1);

[G_xi_tilde ,W_xi] = GaussQuadrature(Gp_xi);
[G_eta_tilde ,W_eta] = GaussQuadrature(Gp_eta);
% [G_xi_tilde ,W_xi] = lgwt(Gp_xi,0,1);
% [G_eta_tilde ,W_eta] = lgwt(Gp_eta,0,1);
W = kron(W_eta,W_xi);

elementStartKnots_xi = findElements(knotVec_xi);
elementStartKnots_eta = findElements(knotVec_eta);

I_neu = 0;

for i = 1:el_eta
    for j = 1:el_xi
        e = (i-1)*el_xi +j;
        k = 0; l = 0;
        
        knot_j = elementStartKnots_xi(j);
        knot_i = elementStartKnots_eta(i);
        
        A_xi = ((knotVec_xi(knot_j+1)- knotVec_xi(knot_j)));
        A_eta = ((knotVec_eta(knot_i+1)- knotVec_eta(knot_i)));
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
        
        
        [N, dNdxi, dNdeta, N_xi, N_eta, dNxidxi,dNetadeta] = shapeFunctions(p_xi, p_eta, G_xi_tilde, G_eta_tilde, Ce, C_xi(:,:,j), C_eta(:,:,i), A_xi, A_eta); % Px(IEN_e)); %, Py(IEN_e));
        [detJ, dNdx, dNdy, dxdxi, dydeta] = Jacobian(dNdxi, dNdeta, Px(IEN_e(:)), Py(IEN_e(:)), Gp);

        detJ;
        % PARAMETER SPACE
        %%%%%%%% Funksjonene er nå i parameter space. De er evaluert for
        %%%%%%%% gausspunkt fra parentspace som er mappet over til
        %%%%%%%% parameter space. Altså trengs ingen mer mapping annet enn
        %%%%%%%% over til physical space nå
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Stiffness matrix
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %gradPhi2g = dNdx(LM_e>0,:)*dNdx(LM_e>0,:)' + dNdy(LM_e>0,:)*dNdy(LM_e>0,:)'
        LM_e;
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
        
        
        %fx = f(N_xi'*Px(IEN_e)*N_eta,N'*Py(IEN_e(:)));
        fx = f(N'*Px(IEN_e(:)),N'*Py(IEN_e(:)));
        for g = 1:Gp
            l = l + fx(g)*N(LM_e>0,g)*W(g)*Jparent*detJ(g);
        end
        
        

        
        
        
        
        

        % Neumann boundary
        if ~isempty(NEU)
            for neu = 1:length(NEU(:,1))
                i_neuMat = find(LM_e(LM_e>0) == NEU(neu,1));
                %i_neuBasis = find(LM_e == NEU(neu,1))
                if ~isempty(i_neuMat)
                    
                    h = cell2mat(h_neu(NEU(neu,2)));
                    LM_e_square = reshape(LM_e,Gp_eta,Gp_xi)';
                    IEN_e_square = reshape(IEN_e,Gp_eta,Gp_xi)';
                    NEU(neu,2);
                    N_xi;
                    N_eta;
                    switch NEU(neu,2)
                        case 1 % left boundary
                            h = h(Px(1,1)*ones(1,Gp_eta),Py(IEN_e_square(:,1))'*N_eta);
                            netaBasis = find(LM_e_square(:,1)== NEU(neu,1));
                            integrand = N_eta(netaBasis,:)*(h'.*W_eta)*dydeta(1)*A_eta;
                        case 2 % right bondary
                            h = h(Px(end,1)*ones(1,Gp_eta),Py(IEN_e_square(:,end))'*N_eta);
                            netaBasis = find(LM_e_square(:,end)== NEU(neu,1));
                            integrand = N_eta(netaBasis,:)*(h'.*W_eta)*dydeta(1)*A_eta;
                        case 3 % buttom boundary
                            h = h(Px(IEN_e_square(1,:))*N_xi, Py(1,1)*ones(1,Gp_xi));
                            nxiBasis = find(LM_e_square(1,:)== NEU(neu,1));
                            integrand = N_xi(nxiBasis,:)*(h'.*W_xi)*dxdxi(1)*A_xi;
                        case 4 % top boundary
                            h = h(Px(IEN_e_square(end,:))*N_xi, Py(IEN_e_square(end,:))*N_eta);
                            nxiBasis = find(LM_e_square(end,:)== NEU(neu,1));
                            integrand = N_xi(nxiBasis,:)*(h'.*W_xi)*dxdxi(1)*A_xi;
                    end
                    Px_e = Px(IEN_e');
                    Py_e = Py(IEN_e');
                    %J_eta = dNetadeta'*Py_e(end,:)';
                    
                    %J_eta
                    dydeta;
                    I_neu = I_neu + integrand*A_xi;
                    l(i_neuMat) = l(i_neuMat) + integrand;
                    %
                    %                 h_neu = {@(y)-dudx_exact(domain.startX*ones(1,np_eta),y), @(y) dudx_exact(domain.endX*ones(1,np_eta),y), ...
                    %        %     @(x)-dudy_exact(x,domain.startY*ones(1,np_xi)), @(x) dudy_exact(x,domain.endY*ones(1,np_xi))};
                    %                 'd'
                    %                 W_xi
                    %                 A_xi
                    %                 J11
                    %                 i_neuBasis
                    %
                    %                 LM_e1 = reshape(LM_e,Gp_eta,Gp_xi)'
                    %N_eta
                    %ID2 = reshape(ID,np_xi,np_eta)'
                    %                 l(i_neuMat)
                    %                 h*N_eta(:,:)*W_eta
                    %
                    %                 %netaBasis = find(ID2(:,end)==N EU(neu,1))
                    %                 netaBasis = find(LM_e1(:,end)== NEU(neu,1))
                    %
                    %                 N_eta(netaBasis,:)
                    %                 h.*W_eta
                    
                    %                 l(i_neuMat) = l(i_neuMat) + N_eta(netaBasis,:)*(h.*W_eta')'*A_xi*J11(1);
                    %
                    %
                    %                 N_xi;
                    %                 h = h(N_xi'*Px(IEN_e))
                    %                 for g = 1:3%Gp
                    %                     g;
                    %                     h(g)
                    %                     N(i_neu,g);
                    %                     N_xi(:,g)
                    %                     l(i_neu) = l(i_neu) + h(g)*N_eta(:,g)*W_xi(g)*A_xi*J11(g);
                    %                 end
                    'SILJE! =D';
                end
            end
            
            
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Nonhomo Dirichlet
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        for dVi = 1:length(DIR)
            if ismember(-dVi, LM_e)
                nDbnd = find(LM_e == -dVi);
                a = 0;
                
                for g = 1:Gp
                    gradPhi2g = dNdx(:,g)*dNdx(:,g)' + dNdy(:,g)*dNdy(:,g)';
                    a = a + gradPhi2g*W(g)*Jparent*detJ(g);
                end
                l = l - DIR(dVi)*sum(a(LM_e>0, nDbnd),2);
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
U = addDirichelBnd(U, ID, DIR);
end


