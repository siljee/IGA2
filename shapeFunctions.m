function [N, dNdxi, dNdeta,N_xi, N_eta]= shapeFunctions(p_xi, p_eta, G_xi, G_eta, Ce, C_xi, C_eta, A_xi, A_eta)
% Make the shape functions for each element. Each row is a basis, and each
% column is the value in that Gauss point.

%%%%%%% PARENT SPACE! %%%%%%%%
% Univariate
% Bernstein is evaluated in parent space
B_xi  = BernsteinBasis(p_xi,G_xi);              % (p+1) x Gpunkt
B_eta = BernsteinBasis(p_eta,G_eta);            % (q+1) x Gpunkt

dB_xi  = BernsteinDerivative(p_xi,G_xi);     
dB_eta = BernsteinDerivative(p_eta,G_eta);       


% Bivariate
B = kron(B_eta,B_xi);                           % (p+1)*(q+1) x Gpunkt

%N - kron(C_eta*B_eta, C_xi*B_xi) %gir 0

dBdxi  = kron(B_eta,dB_xi);
dBdeta = kron(dB_eta,B_xi);


%%%%%%% PARAMETER SPACE %%%%%%%%%
% When multiplied by extraction operator, we are in parameter space!
N = Ce*B;                                       % (p+1)*(q+1) x Gpunkt 
N_xi = C_xi*B_xi;
N_eta = C_eta*B_eta;

% % 
% % Px(:,1);
% % Py(1,:);
% % 
% % kron(N_eta,N_xi)
% % Py;
% % '-';
% 
% 
% N_xi'*Px(:,1);
% N_eta'*Py(1,:)';
% 
% Nx = N_xi'*Px;
% Ny = N_eta'*Py;
% 
% N2 = kron(N_eta'*Py, N_xi'*Px)
A_xi = 1;
A_eta = 1;
%N_eta(:)'*Py(:,1)
%N_xi(:)*Px(1,:)' 
dNdxi  = Ce*dBdxi/A_xi;
dNdeta = Ce*dBdeta/A_eta;
% dNdxi - kron(C_eta*B_eta,C_xi*dB_xi);
% dNdeta - kron(C_eta*dB_eta,C_xi*B_xi);

end