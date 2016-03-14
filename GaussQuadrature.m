function [xq, w] = GaussQuadrature(Nq)  %a,b,Nq,g)
% Calculates the integral of a one dimensional function using the Gauss
% quadrature
%
% a     Start coordinate
% b     End coordinate
% Nq    Number of integration points: 1,2,3 or 4.
% g     The function to be integrated
% 
% I     The approximate value of the integral

%Nq = Nq + 1;

switch Nq
    case 1 % 
        % 1-point-rule
        xq =  0.5;
        w =  1;
    case 2
        % 2-point-rule
        xq = [0.5 - sqrt(3)/6;
              0.5 + sqrt(3)/6];
        w  = [0.5; 
              0.5];
    case 3
        % 3-point-rule
        xq = [0.5 - sqrt(15)/10;
              0.5;
              0.5 + sqrt(15)/10];
        w  = [5/18;
             4/9;
             5/18];
    case 4
        % 4-point-rule
        xq = [0.5 - sqrt(525 + 70*sqrt(30)) / 70;
              0.5 - sqrt(525 - 70*sqrt(30)) / 70;
              0.5 + sqrt(525 - 70*sqrt(30)) / 70;
              0.5 + sqrt(525 + 70*sqrt(30)) / 70];
        w  = [(18-sqrt(30))/72;
              (18+sqrt(30))/72;
              (18+sqrt(30))/72;
              (18-sqrt(30))/72];
    otherwise
        [xq, w] = lgwt(Nq,0,1);
        return
end
% 
% x = xq*(b-a) + a;
% dx = (b-a);
% 
% I = g(x)'*w * dx; 

end