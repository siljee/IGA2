function [errorH0, errorE, Uintegral] = evaluateU(Uw, U_exact, knotVec_xi, knotVec_eta, p_xi, p_eta, el_xi, el_eta, np_xi, np_eta, Px, Py, f, DIR, ID, IEN, LM, h_neu, NEU, dudx_exact, dudy_exact, domain, Gp_xi, Gp_eta)%  leftVal, rightVal);)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%       Initialisere viktige variabler
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


C_xi = BezierExtractionOperator(knotVec_xi,p_xi);
C_eta = BezierExtractionOperator(knotVec_eta,p_eta);


%gof = sum(ID>0);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Det kan være en feil når p ~= q
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Gp_xi = p_xi + 1; Gp_eta = p_eta + 1;    % Number of Gauss points
% Gp = Gp_xi * Gp_eta;
[B_xi, B_eta, dB_xi, dB_eta, W, W_xi, W_eta] = bernsteinBasisAndDerivative(p_xi, p_eta, Gp_xi, Gp_eta)

%U = zeros(size(Uw));
errorH0 = 0;
errorH1 = 0;
errorE = 0;
Uintegral = 0;
% 
% elementStartKnots_xi = findElements(knotVec_xi);
% elementStartKnots_eta = findElements(knotVec_eta);

% Ax = domain.startX;
% Bx = domain.endX;
% lx = (Bx-Ax)/el_xi;
% 
% Ay = domain.startY;
% By = domain.endY;
% ly = (By-Ay)/el_eta;

for i = 1:el_eta
    for j = 1:el_xi
        e = (i-1)*el_xi +j;
        
%         
%         knot_j = elementStartKnots_xi(j);
%         knot_i = elementStartKnots_eta(i);
        
%         Ax_e = (j-1)*lx;
%         Bx_e = j*lx;
%         
%         Ay_e = (i-1)*ly;
%         By_e = (i)*ly;
        
       % Ce = kron(C_eta(:,:,i),C_xi(:,:,j));        % (p+1)*(q+1) x (p+1)*(q+1)
     
          
        
        LM_e = LM(e,:); IEN_e = IEN(e,:);
        inner = LM_e(LM_e > 0);      
        Px_e = Px(IEN_e(:)); Py_e = Py(IEN_e(:));
        
%         IEN_e = IEN(e,:);
%         LM_e = LM(e,:);
%         
%         Px_e = Px(IEN_e(:));
%         Py_e = Py(IEN_e(:));
        

        [Jpar, Jpar_xi, Jpar_eta] = findParentJacobian(knotVec_xi, knotVec_eta, j, i);
        [N, dNdxi, dNdeta, N_xi, N_eta] = shapeFunctions(B_xi, B_eta, dB_xi, dB_eta, C_xi(:,:,j), C_eta(:,:,i), Jpar_xi, Jpar_eta); 
        [detJ, dNdx, dNdy, dxdxi, dydeta] = Jacobian(dNdxi, dNdeta, Px_e, Py_e, Gp_xi*Gp_eta);
%         
%         A_xi = (knotVec_xi(knot_j+1)- knotVec_xi(knot_j));
%         A_eta = (knotVec_eta(knot_i+1)- knotVec_eta(knot_i));
%         A = A_xi*A_eta;
%         J_delta = A;
%         
%         [N, dNdxi, dNdeta,N_xi,N_eta] = shapeFunctions(p_xi, p_eta, G_xi, G_eta, Ce, C_xi(:,:,j), C_eta(:,:,i), A_xi, A_eta); 
%         [detJ, dNdx, dNdy] = Jacobian(dNdxi, dNdeta, Px_e, Py_e, Gp);
%         
%         detJ;   

        % N er evaluert i gausspunkter.
       % disp('----------------- YO! ----------------')
        H0 = 0; H1 = 0; e = 0;
        u = U_exact(N'*Px_e, N'*Py_e);

        dudx = dudx_exact(N'*Px_e, N'*Py_e);
        dudy = dudy_exact(N'*Px_e, N'*Py_e);
        
%         G_xi'*lx + Ax_e
%         U_exact(G_xi'*lx + Ax_e,1) 
% 
%         
%         ((((G_xi/el_xi)+(j-1)/el_xi))*10).^4;
%         
%         N'*Px_e;
%      
        
        uh = N'*Uw(IEN_e);
        
        
%         Uh = N_xi'*reshape(Uw(IEN_e),p_eta+1,p_xi+1)  *N_eta;
        %uh = Uh(:);  
%         detJ'*(((uh-u).^2).*W)*J_delta;
        
        
        
        duhdx = dNdx'*Uw(IEN_e);
        duhdy = dNdy'*Uw(IEN_e);
        du_uhdx = dudx-duhdx;
        du_uhdy = dudy-duhdy;
        
%         H0 = 0;
%         uInt = 0;e 
%         J_delta = 1;
        for g = 1:Gp_xi*Gp_eta
            u_uh = abs(u(g)-uh(g))^2;
            u_uhdx = (dudx(g) - duhdx(g))^2;
            u_uhdy = (dudy(g) - duhdy(g))^2;
%             
%             u_uh * W(g)*detJ(g) *J_delta; 
%             
%             detJ(g);
%             J_delta;
            
            detJ(g);
            Jpar;
            
            H0 = H0 + u_uh * W(g) * abs(detJ(g))*abs(Jpar) ;
            %H1 = H1 + (u_uh + u_uhdx + u_uhdy)*W(g)*detJ(g)*J_delta;
          
            gradPhi2g = du_uhdx(g)*du_uhdx(g)' + du_uhdy(g)*du_uhdy(g)';
            %u_uhdx+u_uhdy -gradPhi2g;
            e = e + (u_uhdx + u_uhdy)*W(g)*abs(Jpar)*abs(detJ(g));
            
            
            %uInt = uInt + u(g) * W(g) * detJ(g)*J_delta ;
%             uInt = uInt + (u(g))^2 * W(g) *detJ(g)*J_delta ;
            %a = a + abs(u(g)-uh(g))^2* W(g);
            %b = b + (abs(u(g)-uh(g))^2 + abs(dudx(g)-duhdx(g))^2 + abs(dudy(g) - duhdy(g))^2)*W(g);
            
        end
   
        
        %H0 = ((((uh-u).^2).*W).*detJ)*J_delta;
        
        %fprintf('\n%.2f \t %.2f \t %.2f \t %.2f\n', Ax_e, Bx_e, Ay_e, By_e)
%         uInt;
%         U_ny = @(x,y)  (U_exact(x,y)).^2;
%         uInt_e = dblquad(U_ny,Ax_e ,Bx_e ,Ay_e ,By_e);
%         uInt_e-uInt;
        
        errorH0 = errorH0 + H0;
       % errorH1 = errorH1 + H1;
        errorE = errorE + e;
        
        
%         Uintegral = Uintegral + uInt;

    end
    
end
%'*'
%errorH1-errorH0-errorE

%errorH1 = sqrt(errorH0 + errorE);
errorH0 = sqrt(errorH0);
errorE = sqrt(errorE);
% 
% N_xi = BsplineBasis(knotVec_xi,p_xi,G_xi');
% N_eta = BsplineBasis(knotVec_eta,p_eta,G_eta');
% 
% N = kron(N_eta,N_xi);
% 
% sin(pi*(N'*Px(:)/10));
% 
% u;
% U = U_exact(N'*Px(:),N'*Py(:));
% Uh = N'*Uw;
% 
% 
% errorH02 = sqrt(sum(((U-Uh).^2).*W)*10);
%detJ'*(((uh-u).^2).*W)*J_delta;

end