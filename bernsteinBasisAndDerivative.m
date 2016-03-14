function [B_xi, B_eta, dB_xi, dB_eta, W, W_xi, W_eta] = bernsteinBasisAndDerivative(p_xi, p_eta, Gp_xi, Gp_eta)
% Generate Bernstein basis and derivative for two dimensions with 
% polynomial order p, evaluated in Gp number of Gauss points in each 
% direction.

[G_xi  ,W_xi]  = GaussQuadrature(Gp_xi);
[G_eta ,W_eta] = GaussQuadrature(Gp_eta);
W = kron(W_eta,W_xi);

B_xi  = BernsteinBasis(p_xi,G_xi);              % (p_xi +1) x Gp
B_eta = BernsteinBasis(p_eta,G_eta);            % (P_eta+1) x Gp

dB_xi  = BernsteinDerivative(p_xi,G_xi);     
dB_eta = BernsteinDerivative(p_eta,G_eta);
end