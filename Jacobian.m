function [detJ, dNdx, dNdy,J11] = Jacobian(dNdxi, dNdeta, Px, Py, Gp)

dNdx = zeros(size(dNdxi));
dNdx1 = zeros(size(dNdxi));
dNdy = zeros(size(dNdxi));

J11   = dNdxi' *Px;
J12   = dNdeta'*Px;
J21   = dNdxi' *Py;
J22   = dNdeta'*Py;

detJ = J11.*J22-J21.*J12;

% Dette kan forbedres
for g = 1:Gp
    J = [J11(g), J12(g); J21(g), J22(g)];
    %invJ = inv(J);
    %dNdx1(:,g) = dNdxi(:,g)*invJ(2,2) + dNdeta(:,g)*invJ(2,1);
    
    dNdxtest = J \ [dNdxi(:,g), dNdeta(:,g)]';
    dNdx(:,g) = dNdxtest(1,:);
    dNdy(:,g) = dNdxtest(2,:);
    
end
end
