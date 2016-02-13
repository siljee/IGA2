close all; clear all; 
addpath HelpFunctions/; addpath Basis/; addpath Problems/;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Define problem variables
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
problem = 'x3';
switch (problem)
    case 'x2'
        U_exact = @(x,y) x.^2;
        dudx_exact = @(x,y) 2*x;
        dudy_exact = @(x,y) zeros(size(x));
        f = @(x,y) repmat(-2, size(x));
    case 'x3'
        U_exact = @(x,y) x.^3;
        dudx_exact = @(x,y) 3*x.^2;
        dudy_exact = @(x,y) zeros(size(x));
        f = @(x,y) -6.*x;
end

startX = 0; endX = 1;
startY = 0; endY = 1;

leftVal = U_exact(startX,0); rightVal = U_exact(endX,0); 
dirichletValues = [ leftVal,  rightVal];

isUniformKV = false;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initialize plot data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

elementsStart = 2; elementsEnd = 9;

errorH0 = zeros(3, elementsEnd-elementsStart);
errorH1 = zeros(3, elementsEnd-elementsStart);
errorE = zeros(3, elementsEnd-elementsStart);
h = zeros(3,elementsEnd-elementsStart);

for p_xi = 1:3
    p_eta = p_xi;
    for el_xi = elementsStart:elementsEnd
        disp('------------------------------------')
        el_eta = el_xi;
        
        if (isUniformKV)
            knotVec_xi  = makeUniformKnotVector(p_xi, el_xi);
            knotVec_eta = makeUniformKnotVector(p_eta,el_eta);
        else
            knotVec_xi  = makeRandomNonUniformKnotVector(el_xi , p_xi , 1, true);
            knotVec_eta = makeRandomNonUniformKnotVector(el_eta, p_eta, 1, true);
        end
        
        np_xi = length(knotVec_xi) - p_xi -1;                   % Number of control points in xi direction
        np_eta = length(knotVec_eta) - p_eta -1;                % Number of control points in eta direction
        
        h(p_xi, el_xi - elementsStart+1) = 1/el_xi;%max(findMaxStepSize(knotVec_xi), findMaxStepSize(knotVec_eta));
        
        IEN = generate_IEN(knotVec_xi, knotVec_eta, p_xi, p_eta);
        ID = generateID_left_right(np_xi, np_eta, -1,-2);
        LM = ID(IEN);
        
        %%%%%% Make control points %%%%%%%%%%%%%
        Px = (findGrevillePoints(knotVec_xi, p_xi)-knotVec_xi(1))* (endX-startX)/(knotVec_xi(end)-knotVec_xi(1)) + startX; Px = repmat(Px',1,np_eta);
        Py = (findGrevillePoints(knotVec_eta, p_eta)-knotVec_eta(1))*(endY-startY)/(knotVec_eta(end)-knotVec_eta(1)) + startY; Py = repmat(Py,np_xi,1);
        
        
        %%%%%% Evaluate on x and y grid %%%%%%%%
        nx = 250; ny = 30;
        xi =  linspace(knotVec_xi(1) ,knotVec_xi(end) ,nx);
        eta = linspace(knotVec_eta(1),knotVec_eta(end),ny);
        
        x = linspace(startX,endX,nx);
        y = linspace(startY,endY,ny);
        xx = x'*ones(1,ny);
        yy = ones(nx,1)*y;
        
        N_xi = BsplineBasis(knotVec_xi, p_xi, xi);
        N_eta = BsplineBasis(knotVec_eta, p_eta, eta);
        N = kron(N_eta,N_xi);
        
        %%%%%%% The method %%%%%%
        U = IGA_Method_Bezier(knotVec_xi, knotVec_eta, p_xi, p_eta, el_xi, el_eta, np_xi, np_eta, Px, Py, f, dirichletValues);%  leftVal, rightVal);
        Uw = U;
        U = N'*Uw;      % Interpolation
        
        Uw = reshape(Uw,np_xi,np_eta);
        U  = reshape(U,nx,ny);
        
        [errorH0(p_xi,el_xi-elementsStart+1), errorH1(p_xi,el_xi-elementsStart+1), errorE(p_xi,el_xi-elementsStart+1)] = evaluateU(Uw, U_exact, knotVec_xi, knotVec_eta, el_xi, el_eta, p_xi, p_eta, Px, Py, f, dirichletValues, dudx_exact, dudy_exact);
        
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


if 1
   
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %       H^0 error
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    figure(4)
    loglog(h',errorH0','linewidth',2)
    hold on
    grid on
    
    diff = 1/3;
    factor12 = (errorH0(1,:)/h(1,:).^(2))*diff;
    factor23 = (errorH0(2,:)/h(2,:).^(3))*diff;
    factor34 = (errorH0(3,:)/h(3,:).^(4))*diff;
    
    fig  = gca;
    set(gca,'fontsize', 14);
    set(fig, 'ColorOrderIndex', 1)
    loglog(h(1,:)',factor12*(h(1,:)'.^(2)),'--')
    loglog(h(2,:)',factor23*(h(2,:)'.^(3)),'--')
    %loglog(h(3,:)',factor34*(h(3,:)'.^(4)),'--')
    
    legend('u_h, p=1','u_h, p=2', 'u_h, p=3', 'C_1h^2 (reference)', 'C_2h^3 (reference)'); % 'C_3h^4 (reference)') % 'u=C_3x³')
    % title('H^0,   u = x^3')
    xlabel('h', 'FontWeight','bold')
    ylabel('||u-u_h||_{H^0}', 'FontWeight','bold')
    % ylabel('H^0 error norm:  ||u-u_h||_{H^0}')
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %       H^1 error
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    figure(5)
    loglog(h',errorH1','linewidth',2)
    hold on
    grid on
    
    factor11 = errorH1(1,:)*diff/h(1,:).^(1);
    factor22 = errorH1(2,:)*diff/h(2,:).^(2);
    factor33 = errorH1(3,:)*diff/h(3,:).^(3);
    
    fig = gca;
    set(fig, 'ColorOrderIndex', 1)
    set(gca,'fontsize', 14);
    loglog(h(1,:)',factor11*(h(1,:)'.^(1)),'--')
    loglog(h(2,:)',factor22*(h(2,:)'.^(2)),'--')
    %loglog(h(3,:)',factor33*(h(3,:)'.^(3)),'--')
    
    legend('u_h, p=1','u_h, p=2', 'u_h, p=3', 'C_1h  (reference)', 'C_2h^2 (reference)'); %  'C_3h^3 (reference)')
    % title('H^1,   u = x^3')
    xlabel('h', 'FontWeight','bold')
    ylabel('||u-u_h||_{H^1}', 'FontWeight','bold')
    % ylabel('H^1 error norm:   ||u-u_h||_{H^1}')
    
    
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
    %loglog(h(3,:)',factor33*(h(3,:)'.^(3)),'--')
    
    legend('u_h, p=1','u_h, p=2', 'u_h, p=3', 'C_1h  (reference)', 'C_2h^2 (reference)'); % 'C_3h^3 (reference)') % 'u=C_3x³')
    % title('Energy error norm,   u = x^3')
    xlabel('h', 'FontWeight','bold')
    ylabel('||u-u_h||_{E}', 'FontWeight','bold')
    % ylabel('Energy error norm:   ||u-u_h||_{E}')
    %axis([max(min(min(h))-max(max(h))/2,0), max(max(max(h))+max(max(h))/2,0), max(min(min(errorE))-max(max(errorE))/2,0), max(max(max(errorE))+max(max(errorE))/2,0)])
end