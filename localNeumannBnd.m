function l = localNeumannBnd(NEU, h_neu, Px, Py, N_xi, N_eta, dxdxi, dydeta, A_xi, A_eta, LM_e, IEN_e, Gp_xi, Gp_eta, W_xi, W_eta)
inner = LM_e(LM_e > 0);
l = zeros(size(inner'));

% Neumann boundary
if ~isempty(NEU)    
    for neu = 1:length(NEU(:,1));
        i_neuMat = find(inner == NEU(neu,1));
        %i_neuBasis = find(LM_e == NEU(neu,1))
        if ~isempty(i_neuMat)
            h = cell2mat(h_neu(NEU(neu,2)));
            LM_e_square = reshape(LM_e,Gp_eta,Gp_xi)';
            IEN_e_square = reshape(IEN_e,Gp_eta,Gp_xi)';

            switch NEU(neu,2)
                case 1 % left boundary
                    x = Px(IEN_e_square(:,1))' * N_xi;
                    y = Py(IEN_e_square(:,1))' * N_eta;
                    N_index = find(LM_e_square(:,1)== NEU(neu,1));
                    I = lineIntegral(h(x,y), N_eta, N_index, W_eta, dydeta(1)*A_eta);
                    %h = h(Px(IEN_e_square(:,1))' * N_xi,Py(IEN_e_square(:,1))'*N_eta);
                    %netaBasis = find(LM_e_square(:,1)== NEU(neu,1));
                    %integral = N_eta(netaBasis,:)*(h'.*W_eta)*dydeta(1)*A_eta;
                case 2 % right bondary
                    x = Px(IEN_e_square(:,end))' * N_xi; 
                    y = Py(IEN_e_square(:,end))' * N_eta;
                    %h = h(Px(IEN_e_square(:,end))'*N_xi,Py(IEN_e_square(:,end))'*N_eta);
                    N_index = find(LM_e_square(:,end)== NEU(neu,1));
                    I = lineIntegral(h(x,y), N_eta, N_index, W_eta, dydeta(1)*A_eta);
                    %integral = N_eta(netaBasis,:)*(h'.*W_eta)*dydeta(1)*A_eta;
                case 3 % buttom boundary
                    x = Px(IEN_e_square(1,:)) * N_xi;
                    y = Py(IEN_e_square(1,:)) * N_eta;
                    %h = h(Px(IEN_e_square(1,:))*N_xi, Py(1,1)*ones(1,Gp_xi));
                    N_index = find(LM_e_square(1,:)== NEU(neu,1));
                    I = lineIntegral(h(x,y), N_xi, N_index, W_xi, dxdxi(1)*A_xi);
                    %integral = N_xi(nxiBasis,:)*(h'.*W_xi)*dxdxi(1)*A_xi;
                case 4 % top boundary
                    x = Px(IEN_e_square(end,:)) * N_xi;
                    y = Py(IEN_e_square(end,:)) * N_eta;
%                     h = h(Px(IEN_e_square(end,:))*N_xi, );
                    N_index = find(LM_e_square(end,:)== NEU(neu,1));
                    I = lineIntegral(h(x,y), N_xi, N_index, W_xi, dxdxi(1)*A_xi);
                    %integral = N_xi(nxiBasis,:)*(h'.*W_xi)*dxdxi(1)*A_xi;
            end
            l(i_neuMat) = l(i_neuMat) + I;
            
        end
    end
end
end

function I = lineIntegral(h, N, N_index, W, J)  
I = N(N_index,:)*(h'.*W)*J;
end