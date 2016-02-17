function U = addDirichelBnd(U0, ID, dirVal)
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

U = zeros(length(ID),1);
U(ID>0) = U0;

for i = 1:length(dirVal)
    U(ID==-i) = dirVal(i);
end

