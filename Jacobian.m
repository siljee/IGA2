function [detJ, dNdx, dNdy, dxdxi, dydeta] = Jacobian(dNdxi, dNdeta, Px, Py, Gp)

dNdx = zeros(size(dNdxi));
dNdy = zeros(size(dNdxi));

dxdxi   = dNdxi' *Px;
dxdeta  = dNdeta'*Px;
dydxi   = dNdxi' *Py;
dydeta  = dNdeta'*Py;

detJ = dxdxi.*dydeta-dydxi.*dxdeta; %OK(-)

% Dette kan forbedres
for g = 1:Gp
    J = [dxdxi(g), dxdeta(g); dydxi(g), dydeta(g)];
    %invJ = inv(J);
    %dNdx1(:,g) = dNdxi(:,g)*invJ(2,2) + dNdeta(:,g)*invJ(2,1);
    
    dNdxtest = J' \ [dNdxi(:,g), dNdeta(:,g)]';
    dNdx(:,g) = dNdxtest(1,:);     %OK
    %dNdxi(:,g)*dxdxi(g)) + dNdeta(:,g)*inv(dxdeta(g))
    dNdy(:,g) = dNdxtest(2,:);     %OK
end
end
