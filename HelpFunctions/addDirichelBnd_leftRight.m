function U = addDirichelBnd_leftRight(U0, np_xi, np_eta, leftVal, rightVal)
% This function adds Dirichlet boundary values directly to the solution U. 
% The boundary values are only added to the left and right side on the 
% parameter space. 
%
% Input:    
%   U0       - The initial solution U without entries for the Dirichlet boundary.
%   np_xi    - Number of control points in xi-direction.
%   np_eta   - Number of control points in eta-direction.
%   leftVal  - The Dirichlet boundary value on the left side of the parameter space.
%   rightVal - The Dirichlet boundary value on the right side of the parameter space.
%
% Output:
%   U        - The updated solution U, with values of the Dirichlet boundary 
%              on the left and right side of the parameter space.

U = [];
for i = 1:np_eta
    U = [U; leftVal; U0((i-1)*(np_xi-2)+1:i*(np_xi-2)); rightVal];
end
end