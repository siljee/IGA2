function l = localLoadVector(f, N, Px, Py, LM_e, Gp, W, J)
l = 0;
f = f(N'*Px,N'*Py);
for g = 1:Gp
    l = l + f(g)*N(LM_e>0,g)*W(g)*J(g);
end
end