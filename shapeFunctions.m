function [N, dNdxi, dNdeta,N_xi, N_eta]= shapeFunctions(B_xi, B_eta, dB_xi, dB_eta, C_xi, C_eta, Jparent_xi, Jparent_eta)
% Make the shape functions for each element. Each row is a basis, and each
% column is the value in that Gauss point.

% Elemental Bèzier extraction operator
Ce = kron(C_eta,C_xi);

% B-splines
N_xi  = C_xi *B_xi;  
N_eta = C_eta*B_eta; 
N = kron(N_eta, N_xi);

% Derivative
dBdxi  = kron(B_eta,dB_xi);
dBdeta = kron(dB_eta,B_xi);

dNdxi  = Ce*dBdxi /Jparent_xi;   
dNdeta = Ce*dBdeta/Jparent_eta; 

end