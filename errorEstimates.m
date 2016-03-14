function [h, errorE, errorH0, U,  p_xi, p_eta, knotVec_xi, knotVec_eta, Px, Py] = errorEstimates(p_xi_start, p_eta_start, knotVec_xi, knotVec_eta, Px, Py, U_exact, f, dudx_exact, dudy_exact, bnd, domain, h_neu, nohr, nopr, shape)
p_xi = p_xi_start;
p_eta = p_eta_start;
el_xi_start = countElements(knotVec_xi);
el_eta_start = countElements(knotVec_eta);

errorH0 = zeros(3, nohr);
errorH1 = zeros(3, nohr);
errorE = zeros(3, nohr);
h = zeros(3,nohr);

isUniformKV = true;

for pr = 1:nopr     % p-refinement 
    pr

    for hr = 1:nohr     % h-refinement
        if hr == 1
            knotVec_xi_last = knotVec_xi;
            knotVec_eta_last = knotVec_eta;
            Px_last = Px;
            Py_last = Py;
        else
             [knotVec_xi, knotVec_eta, Px, Py] = h_refinement(knotVec_xi, knotVec_eta ,p_xi, p_eta, Px, Py);
        end
        Gp_xi = p_xi+1; Gp_eta = p_eta +1;    % Number of Gauss points
        Gp = Gp_xi * Gp_eta;
        
        hr
        disp('------------------------------------')
        %el_eta = el_xi; %el_xi;
        
        %[knotVec_xi, knotVec_eta, Px, Py] = h_refinement(knotVec_xi, knotVec_eta ,p_xi, p_eta, Px, Py);
        
        el_xi = countElements(knotVec_xi);
        el_eta = countElements(knotVec_eta);
        
        np_xi  = length(knotVec_xi)  - p_xi  - 1;                 % Number of control points in xi direction
        np_eta = length(knotVec_eta) - p_eta - 1;                 % Number of control points in eta direction
        
        
        
        h(p_xi, hr) = findMaxRadius(Px,Py); %1/el_xi;%max(findMaxStepSize(knotVec_xi), findMaxStepSize(knotVec_eta));
        
        
        
        
        
        
        
        [ID, IEN, LM, DIR, NEU] = generateDataArrays(knotVec_xi, knotVec_eta, p_xi, p_eta, bnd, U_exact, domain, Px, Py);
        
        %%%%%% Evaluate on x and y grid %%%%%%%%
        nx = 5; ny = 4;
        xi =  linspace(knotVec_xi(1) ,knotVec_xi(end) ,nx);
        eta = linspace(knotVec_eta(1),knotVec_eta(end),ny);
        
        %x = linspace(domain.startX,domain.endX,nx);
        %y = linspace(domain.startY,domain.endY,ny);
        
        
        N_xi = BsplineBasis(knotVec_xi, p_xi, xi);
        N_eta = BsplineBasis(knotVec_eta, p_eta, eta);
        N = kron(N_eta,N_xi);

        
        
        
        %%%%%%% The method %%%%%%
        
        U = IGA_Method_Bezier(knotVec_xi, knotVec_eta, p_xi, p_eta, el_xi, el_eta, np_xi, np_eta, Px, Py, f, DIR, ID, IEN, LM, h_neu, NEU, Gp_xi, Gp_eta);%  leftVal, rightVal);
        Uw = U;
        U = N'*Uw;      % Interpolation
        
        
        
        [errorH0(p_xi,hr), errorE(p_xi,hr), uInt] = evaluateU(Uw, U_exact, knotVec_xi, knotVec_eta, p_xi, p_eta, el_xi, el_eta, np_xi, np_eta, Px, Py, f, DIR, ID, IEN, LM, h_neu, NEU, dudx_exact, dudy_exact, domain, Gp_xi+10, Gp_eta+10);
        U_ny = @(x,y) (U_exact(x,y)).^2;
        %         uInt_exact = dblquad(U_ny,domain.startX,domain.endX,domain.startY, domain.endY)
        %         uInt
        %         uInt - uInt_exact
        
        
        Uw = reshape(Uw,np_xi,np_eta);
        U  = reshape(U,nx,ny);
        
        N_knot_xi = BsplineBasis(knotVec_xi,p_xi, knotVec_xi(1+p_xi:end-p_xi));
        N_knot_eta = BsplineBasis(knotVec_eta,p_eta, knotVec_eta(1+p_eta:end-p_eta));
        
    end
    
    if pr ~= nopr
        
        %[Px, Py] = p_refinement(knotVec_xi_last, knotVec_eta_last, p_xi, p_eta, Px_last, Py_last);
        p_xi = p_xi +1;
        p_eta = p_eta +1;
        
        if (isUniformKV)
            knotVec_xi  = makeUniformKnotVector(p_xi, el_xi_start);
            knotVec_eta = makeUniformKnotVector(p_eta,el_eta_start);
        else
            knotVec_xi  = makeRandomNonUniformKnotVector(el_xi_start , p_xi ,2, true);
            knotVec_eta = makeRandomNonUniformKnotVector(el_eta_start, p_eta, 2, true);
        end
        
        np_xi  = length(knotVec_xi)  - p_xi  - 1;                 % Number of control points in xi direction
        np_eta = length(knotVec_eta) - p_eta - 1;                 % Number of control points in eta direction

        
        Px = (findGrevillePoints(knotVec_xi, p_xi)-knotVec_xi(1))* (domain.endX-domain.startX)/(knotVec_xi(end)-knotVec_xi(1)) + domain.startX; Px = repmat(Px',1,np_eta);
        Py = (findGrevillePoints(knotVec_eta, p_eta)-knotVec_eta(1))*(domain.endY-domain.startY)/(knotVec_eta(end)-knotVec_eta(1)) + domain.startY; Py = repmat(Py,np_xi,1);
        
        switch shape
            case 'parall-x'
                Px = Px+Py;
            case 'parall-y'
                Py = Px+Py;
                
                n_y = Px(end,1) - Px(1,1);
                n_x = Py(1,end) - Py(end,end);
                
                h_neu{3} = @(x,y) -(dudx_exact(x,y)*n_x + dudy_exact(x,y)*n_y);
                h_neu{4} = @(x,y)   dudx_exact(x,y)*n_x + dudy_exact(x,y)*n_y;
        end
        
    end
    
    size(Px)
    size(Py)
end
