function k = localStiffnessMatrix(dNdx, dNdy, Gp, W, J)
k = 0;
for g = 1:Gp
    integrand = dNdx(:,g)*dNdx(:,g)' + dNdy(:,g)*dNdy(:,g)';
    k = k + integrand*W(g)*J(g); 
end
end