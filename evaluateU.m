function [error0, error1, errorE] = evaluateU(Uw, U_exact, knotVec_xi, knotVec_eta, el_xi, el_eta, p_xi, p_eta, Px, Py, f, dirichletValues, dudx_exact, dudy_exact)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%       Initialisere viktige variabler
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

np_xi = length(knotVec_xi) - p_xi -1;                   % Number of control points in xi direction
np_eta = length(knotVec_eta) - p_eta -1;                % Number of control points in eta direction

C_xi = BezierExtractionOperator(knotVec_xi,p_xi);
C_eta = BezierExtractionOperator(knotVec_eta,p_eta);

IEN = generate_IEN(knotVec_xi, knotVec_eta, p_xi, p_eta);
ID = generateID_left_right(np_xi, np_eta, -1, -2);
LM = ID(IEN);

%gof = sum(ID>0);
Gp_xi = p_xi + 1; Gp_eta = p_eta + 1;    % Number of Gauss points
Gp = Gp_xi * Gp_eta;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Det kan være en feil når p ~= q
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[G_xi,W_xi] = GaussQuadrature(Gp_xi);
[G_eta,W_eta] = GaussQuadrature(Gp_eta);
W = kron(W_eta,W_xi);
%U = zeros(size(Uw));
error0 = 0;
error1 = 0;
errorE = 0;

elementStartKnots_xi = findElements(knotVec_xi);
elementStartKnots_eta = findElements(knotVec_eta);


for i = 1:el_eta
    for j = 1:el_xi
        e = (i-1)*el_xi +j;
        
        knot_j = elementStartKnots_xi(j);
        knot_i = elementStartKnots_eta(i);
        
       
        Ce = kron(C_eta(:,:,i),C_xi(:,:,j));        % (p+1)*(q+1) x (p+1)*(q+1)
      
        IEN_e = IEN(e,:);
%         LM_e = LM(e,:);
        
        Px_e = Px(IEN_e(:));
        Py_e = Py(IEN_e(:));
        
        A_xi = (knotVec_xi(knot_j+1)- knotVec_xi(knot_j));
        A_eta = (knotVec_eta(knot_i+1)- knotVec_eta(knot_i));
        A = A_xi*A_eta;
        J_delta = A;
        
        [N, dNdxi, dNdeta] = shapeFunctions(p_xi, p_eta, G_xi, G_eta, Ce, C_xi(:,:,j), C_eta(:,:,i), A_xi, A_eta); 
        [detJ, dNdx, dNdy] = Jacobian(dNdxi, dNdeta, Px(IEN_e(:)), Py(IEN_e(:)), Gp);
        
        

        % N er evaluert i gausspunkter.
       % disp('----------------- YO! ----------------')
        a = 0; b = 0; e = 0;
        u = U_exact(N'*Px_e, N'*Py_e);
        dudx = dudx_exact(N'*Px_e, N'*Py_e);
        dudy = dudy_exact(N'*Px_e, N'*Py_e);
    
        uh = N'*Uw(IEN_e)';
        duhdx = dNdx'*Uw(IEN_e)';
        duhdy = dNdy'*Uw(IEN_e)';
        du_uhdx = dudx-duhdx;
        du_uhdy = dudy-duhdy;

        for g = 1:Gp
            u_uh = (u(g)-uh(g))^2;
            u_uhdx = (dudx(g)-duhdx(g))^2;
            u_uhdy = (dudy(g) - duhdy(g))^2;
            
            a = a + u_uh * W(g) * detJ(g) * J_delta;
            b = b + (u_uh + u_uhdx + u_uhdy)*W(g)*detJ(g)*J_delta;
          
            gradPhi2g = du_uhdx(g)*du_uhdx(g)' + du_uhdy(g)*du_uhdy(g)';
            e = e + gradPhi2g*W(g)*J_delta*detJ(g);
            
            %a = a + abs(u(g)-uh(g))^2* W(g);
            %b = b + (abs(u(g)-uh(g))^2 + abs(dudx(g)-duhdx(g))^2 + abs(dudy(g) - duhdy(g))^2)*W(g);
        end
        
        error0 = error0 + sqrt(a);
        error1 = error1 + sqrt(b);
        errorE = errorE + sqrt(e);

    end
end
end