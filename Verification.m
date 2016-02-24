close all; clear all;
addpath HelpFunctions/; addpath Basis/; addpath tests/

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Proble specification - Start here!
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Welcome!
%
% Please specify your problem (choose a string from the list below):
problem = 'y2';
%
% Specify your initial domain:
domain.startX = 0; domain.endX = 1;
domain.startY = 0; domain.endY = 1;
%
% In what shape do you want to shape the domain?
% shape = 'rect'     % 'rect'     : Rectangle.
% shape = 'parall-x' % 'parall-x' : Parallelogram in x-direction.
shape = 'parall-y' % 'parall-y' : Parallelogram in y-direction.

%
% Now, to boundary conditions!
% What type are your boundaries? ('d'-Dirichlet, 'n'-Neumann)
bnd.left.type   = 'd';
bnd.right.type  = 'd';
bnd.buttom.type = 'd';
bnd.top.type    = 'n';
%
% Are any of these homogeneous('h'), constant ('c') or variable('v')?
bnd.left.value   = 'v';
bnd.right.value  = 'v';
bnd.buttom.value = 'v';
bnd.top.value    = 'v';
%
% Thanks, that's all! ^_^
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

elementsStart = 2; elementsEnd = 10;
pStart = 1; pEnd = 3;

isUniformKV = true;
isErrorPlot = 1;
diff = 0.80;

% Problems to choose from:
switch (problem)
    case 'x'
        U_exact = @(x,y) x;
        dudx_exact = @(x,y) ones(size(x));
        dudy_exact = @(x,y) zeros(size(x));
        f = @(x,y) zeros(size(x));
    case 'y'
        U_exact = @(x,y) y;
        dudx_exact = @(x,y) zeros(size(y));
        dudy_exact = @(x,y) ones(size(y));
        f = @(x,y) zeros(size(y));
    case 'x+y'
        U_exact = @(x,y) x+y;
        dudx_exact = @(x,y) ones(size(x));
        dudy_exact = @(x,y) ones(size(x));
        f = @(x,y) zeros(size(x));
    case 'x2'
        U_exact = @(x,y) x.^2;
        dudx_exact = @(x,y) 2*x;
        dudy_exact = @(x,y) zeros(size(y));
        f = @(x,y) repmat(-2, size(x));
    case 'y2'
        U_exact = @(x,y) y.^2;
        dudx_exact = @(x,y) zeros(size(x));
        dudy_exact = @(x,y) 2*y;
        f = @(x,y) repmat(-2, size(y));
    case 'x2+y'
        U_exact = @(x,y) x.^2+y;
        dudx_exact = @(x,y) 2*x;
        dudy_exact = @(x,y) ones(size(x));
        f = @(x,y) repmat(-2, size(y));
    case 'x3'
        U_exact = @(x,y) x.^3;
        dudx_exact = @(x,y) 3*x.^2;
        dudy_exact = @(x,y) zeros(size(y));
        f = @(x,y) -6.*x;
    case 'y3'
        U_exact = @(x,y) y.^3;
        dudx_exact = @(x,y) zeros(size(x));
        dudy_exact = @(x,y) 3*y.^2;
        f = @(x,y) -6.*y;
    case 'x3+y2'
        U_exact    = @(x,y) x.^3 + y.^2;
        dudx_exact = @(x,y) 3*x.^2;
        dudy_exact = @(x,y) 2*y;
        f = @(x,y) -6.*x - 2;
    case 'x4'
        U_exact = @(x,y) x.^4;
        dudx_exact = @(x,y) 4*x.^3;
        dudy_exact = @(x,y) zeros(size(y));
        f = @(x,y) -12.*x.^2;
    case 'x3+y4'
        U_exact    = @(x,y) x.^3 + y.^4;
        dudx_exact = @(x,y) 3*x.^2;
        dudy_exact = @(x,y) 4*y.^3;
        f = @(x,y) -6.*x - 12.*y.^2;
    case 'sin(pix)'
        U_exact = @(x,y) sin(pi*x);
        dudx_exact = @(x,y) pi*cos(pi*x);
        dudy_exact = @(x,y) zeros(size(y));
        f = @(x,y) pi^2*sin(pi*x);
    case 'sin(piy)'
        U_exact = @(x,y) sin(pi*y);
        dudx_exact = @(x,y) zeros(size(x));
        dudy_exact = @(x,y) pi*cos(pi*y);
        f = @(x,y) pi^2*sin(pi*y);
    case 'sin(pix/10)'
        U_exact = @(x,y) sin(pi*x/10);
        dudx_exact = @(x,y) pi*cos(pi*x/10)/10;
        dudy_exact = @(x,y) zeros(size(y));
        f = @(x,y) pi^2*sin(pi*x/10)/100;
end

% Works only for squared domaain
h_neu = {@(x,y)-dudx_exact(x,y), dudx_exact, ...
    @(x,y)-dudy_exact(x,y), dudy_exact};

% % Constant Dirichlet values
% leftVal = U_exact(startX,0); rightVal = U_exact(endX,0);
% dirichletValues = [ leftVal,  rightVal];



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initialize plot data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



errorH0 = zeros(3, elementsEnd-elementsStart);
errorH1 = zeros(3, elementsEnd-elementsStart);
errorE = zeros(3, elementsEnd-elementsStart);
h = zeros(3,elementsEnd-elementsStart);

for p_xi = pStart:pEnd
    p_eta = p_xi;
    for el_xi = elementsStart:elementsEnd
        disp('------------------------------------')
        el_eta = el_xi; %el_xi;
        
        if (isUniformKV)
            knotVec_xi  = makeUniformKnotVector(p_xi, el_xi);
            knotVec_eta = makeUniformKnotVector(p_eta,el_eta);
        else
            knotVec_xi  = makeRandomNonUniformKnotVector(el_xi , p_xi , 1, true);
            knotVec_eta = makeRandomNonUniformKnotVector(el_eta, p_eta, 1, true);
        end
        
        
        np_xi  = length(knotVec_xi)  - p_xi  - 1;                  % Number of control points in xi direction
        np_eta = length(knotVec_eta) - p_eta - 1;                % Number of control points in eta direction
        
        % Works only for squared domains
        % h_neu = {@(y)-dudx_exact(domain.startX*ones(1,np_eta),y), @(y) dudx_exact(domain.endX*ones(1,np_eta),y), ...
        %     @(x)-dudy_exact(x,domain.startY*ones(1,np_xi)), @(x) dudy_exact(x,domain.endY*ones(1,np_xi))};
        
        
        
        
        h(p_xi, el_xi - elementsStart+1) = 1/el_xi;%max(findMaxStepSize(knotVec_xi), findMaxStepSize(knotVec_eta));
        
        %%%%%% Make control points %%%%%%%%%%%%%
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
        
        
                
        
        
        [ID, IEN, LM, DIR, NEU] = generateDataArrays(knotVec_xi, knotVec_eta, p_xi, p_eta, bnd, U_exact, domain, Px, Py);
        
        %%%%%% Evaluate on x and y grid %%%%%%%%
        nx = 5; ny = 4;
        xi =  linspace(knotVec_xi(1) ,knotVec_xi(end) ,nx);
        eta = linspace(knotVec_eta(1),knotVec_eta(end),ny);
        
        x = linspace(domain.startX,domain.endX,nx);
        y = linspace(domain.startY,domain.endY,ny);
        
        
        N_xi = BsplineBasis(knotVec_xi, p_xi, xi);
        N_eta = BsplineBasis(knotVec_eta, p_eta, eta);
        N = kron(N_eta,N_xi);
        xx = reshape(N'*Px(:),nx,ny);
        yy = reshape(N'*Py(:),nx,ny);
        
        
        
        %%%%%%% The method %%%%%%
        
        U = IGA_Method_Bezier(knotVec_xi, knotVec_eta, p_xi, p_eta, el_xi, el_eta, np_xi, np_eta, Px, Py, f, DIR, ID, IEN, LM, h_neu, NEU);%  leftVal, rightVal);
        Uw = U;
        U = N'*Uw;      % Interpolation
        
        
        
        [errorH0(p_xi,el_xi-elementsStart+1), errorE(p_xi,el_xi-elementsStart+1), uInt] = evaluateU(Uw, U_exact, knotVec_xi, knotVec_eta, p_xi, p_eta, el_xi, el_eta, np_xi, np_eta, Px, Py, f, DIR, ID, IEN, LM, h_neu, NEU, dudx_exact, dudy_exact, domain);
        U_ny = @(x,y) (U_exact(x,y)).^2;
        %         uInt_exact = dblquad(U_ny,domain.startX,domain.endX,domain.startY, domain.endY)
        %         uInt
        %         uInt - uInt_exact
        
        
        Uw = reshape(Uw,np_xi,np_eta);
        U  = reshape(U,nx,ny);
        
        N_knot_xi = BsplineBasis(knotVec_xi,p_xi, knotVec_xi(1+p_xi:end-p_xi));
        N_knot_eta = BsplineBasis(knotVec_eta,p_eta, knotVec_eta(1+p_eta:end-p_eta));
        
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Plot
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(1)
subplot(3,1,1)
surf(Px,Py,U_exact(Px,Py))% 'linestyle', 'none');  % Det kan være Py skal transponeres!
title(['Exact.   p_\xi=', num2str(p_xi), ', p_\eta=', num2str(p_eta)])
hold on
%plot((N_knot_xi'*Px)', 'k-*')
%plot((N_knot_eta'*Py')','k-*')

subplot(3,1,2)
surf(xx, yy,U)
title(['Approximate.   p_\xi=', num2str(p_xi), ', p_\eta=', num2str(p_eta)])
hold on
plot3(Px,Py,Uw,'*')
xlabel('x')
ylabel('y')

subplot(3,1,3)
surf(xx, yy, U-U_exact(xx,yy))
title(['Error.   p_\xi=', num2str(p_xi), ', p_\eta=', num2str(p_eta)])


if isErrorPlot
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %       H^0 error
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    figure(4)
    loglog(h',errorH0','linewidth',2)
    hold on
    grid on
    
    
    %diff = 0.9;
    factor12 = (errorH0(1,:)/h(1,:).^(2))*diff;
    factor23 = (errorH0(2,:)/h(2,:).^(3))*diff;
    factor34 = (errorH0(3,:)/h(3,:).^(4))*diff;
    
    fig  = gca;
    set(gca,'fontsize', 14);
    set(fig, 'ColorOrderIndex', 1)
    loglog(h(1,:)',factor12*(h(1,:)'.^(2)),'--')
    loglog(h(2,:)',factor23*(h(2,:)'.^(3)),'--')
    loglog(h(3,:)',factor34*(h(3,:)'.^(4)),'--')
    
    legend('u_h, p=1','u_h, p=2', 'u_h, p=3', 'C_1h^2 (reference)', 'C_2h^3 (reference)', 'C_3h^4 (reference)') % 'u=C_3x³')
    % title('H^0,   u = x^3')
    xlabel('h', 'FontWeight','bold')
    ylabel('||u-u_h||_{H^0}', 'FontWeight','bold')
    % ylabel('H^0 error norm:  ||u-u_h||_{H^0}')
    
    %
    %     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %     %       H^1 error
    %     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %     figure(5)
    %     loglog(h',errorH1','linewidth',2)
    %     hold on
    %     grid on
    %
    %     factor11 = errorH1(1,:)*diff/h(1,:).^(1);
    %     factor22 = errorH1(2,:)*diff/h(2,:).^(2);
    %     factor33 = errorH1(3,:)*diff/h(3,:).^(3);
    %
    %     fig = gca;
    %     set(fig, 'ColorOrderIndex', 1)
    %     set(gca,'fontsize', 14);
    %     loglog(h(1,:)',factor11*(h(1,:)'.^(1)),'--')
    %     loglog(h(2,:)',factor22*(h(2,:)'.^(2)),'--')
    %     loglog(h(3,:)',factor33*(h(3,:)'.^(3)),'--')
    %
    %     legend('u_h, p=1','u_h, p=2', 'u_h, p=3', 'C_1h  (reference)', 'C_2h^2 (reference)',  'C_3h^3 (reference)')
    %     % title('H^1,   u = x^3')
    %     xlabel('h', 'FontWeight','bold')
    %     ylabel('||u-u_h||_{H^1}', 'FontWeight','bold')
    %     % ylabel('H^1 error norm:   ||u-u_h||_{H^1}')
    %
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %       Energy error
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    figure(6)
    loglog(h',errorE','linewidth',2)
    hold on
    grid on
    
    factor11 = errorE(1,:)*diff/h(1,:).^(1);
    factor22 = errorE(2,:)*diff/h(2,:).^(2);
    factor33 = errorE(3,:)*diff/h(3,:).^(3);
    fig = gca;
    set(gca,'fontsize', 14);
    set(fig, 'ColorOrderIndex', 1)
    loglog(h(1,:)',factor11*(h(1,:)'.^(1)),'--')
    loglog(h(2,:)',factor22*(h(2,:)'.^(2)),'--')
    loglog(h(3,:)',factor33*(h(3,:)'.^(3)),'--')
    
    legend('u_h, p=1','u_h, p=2', 'u_h, p=3', 'C_1h  (reference)', 'C_2h^2 (reference)', 'C_3h^3 (reference)') % 'u=C_3x³')
    % title('Energy error norm,   u = x^3')
    xlabel('h', 'FontWeight','bold')
    ylabel('||u-u_h||_{E}', 'FontWeight','bold')
    % ylabel('Energy error norm:   ||u-u_h||_{E}')
    %axis([max(min(min(h))-max(max(h))/2,0), max(max(max(h))+max(max(h))/2,0), max(min(min(errorE))-max(max(errorE))/2,0), max(max(max(errorE))+max(max(errorE))/2,0)])
end

