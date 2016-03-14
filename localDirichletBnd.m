function l = localDirichletBnd(DIR, dNdx, dNdy, Gp, W, J, LM_e)
inner = LM_e(LM_e > 0);
l = zeros(size(inner'));


for dVi = 1:length(DIR)
    if ismember(-dVi, LM_e)
        nDbnd = find(LM_e == -dVi);
        a = 0;
        
        for g = 1:Gp
            gradPhi2g = dNdx(:,g)*dNdx(:,g)' + dNdy(:,g)*dNdy(:,g)';
            a = a + gradPhi2g*W(g)*J(g);
        end
        l = l - DIR(dVi)*sum(a(LM_e>0, LM_e == -dVi),2);
    end
end
end