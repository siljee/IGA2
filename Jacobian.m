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
    
    dNdxtest = J' \ [dNdxi(:,g), dNdeta(:,g)]';
    dNdx(:,g) = dNdxtest(1,:);     %OK
    dNdy(:,g) = dNdxtest(2,:);     %OK
end
end
