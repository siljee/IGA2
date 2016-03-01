close all; clear all;
addpath HelpFunctions/; addpath Basis/; addpath tests/

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Proble specification - Start here!
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Welcome!
%
% Please specify your problem (choose a string from the list below):
problem = 'x2';
%
% Specify your initial domain:
domain.startX = 0; domain.endX = 1;
domain.startY = 0; domain.endY = 1;
%
% In what shape do you want to shape the domain?
% ('rect' - rectangle,
% 'parall-x' parallelogram in x-direction,
% 'parall-y' - parallelogram in y-direction,
% 'plateHole' - square plate with a hole)
% 'L-shape' 
%
% shape = 'rect';
% shape = 'rect-8';
% shape = 'parall-x';     % Only buttom and top can be neumann, for now.
% shape = 'parall-y';     % Only left and right can be neumann, for now.
shape = 'plateHole';
% shape = 'L-shape';

%
% Now, to boundary conditions!
% What type are your boundaries? ('d'-Dirichlet, 'n'-Neumann)
bnd.left.type   = 'd';
bnd.right.type  = 'd';
bnd.buttom.type = 'd';
bnd.top.type    = 'd';
%
% Are any of these homogeneous('h'), constant ('c') or variable('v')?
bnd.left.value   = 'v';
bnd.right.value  = 'v';
bnd.buttom.value = 'v';
bnd.top.value    = 'v';
%
% Thanks, that's all! ^_^
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

el_xi_start = 2; el_eta_start = 2;
nohr = 1;  % number of h-refinement

%elementsStart = 2; elementsEnd = 10;
p_xi_start = 1; p_eta_start = 1;
nopr = 1; % number of p-refinement

isUniformKV = true;
isErrorPlot = 1;
diff = 0.60;

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
    case 'x6'
        U_exact = @(x,y) x.^6;
        dudx_exact = @(x,y) 6*x.^5;
        dudy_exact = @(x,y) zeros(size(y));
        f = @(x,y) -30.*x.^4;
    case 'x8'
        U_exact = @(x,y) x.^8;
        dudx_exact = @(x,y) 8*x.^7;
        dudy_exact = @(x,y) zeros(size(y));
        f = @(x,y) -56.*x.^6;
    case 'x9'
        U_exact = @(x,y) x.^9;
        dudx_exact = @(x,y) 9*x.^8;
        dudy_exact = @(x,y) zeros(size(y));
        f = @(x,y) -72.*x.^7;
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

%%%%%% Make control points %%%%%%%%%%%%%



errorH0 = zeros(3, nohr);
errorH1 = zeros(3, nohr);
errorE = zeros(3, nohr);
h = zeros(3,nohr);

p_xi = p_xi_start; p_eta = p_eta_start;

if (isUniformKV)
    knotVec_xi  = makeUniformKnotVector(p_xi_start, el_xi_start);
    knotVec_eta = makeUniformKnotVector(p_eta_start,el_eta_start);
else
    knotVec_xi  = makeRandomNonUniformKnotVector(el_xi_start , p_xi_start ,2, true);
    knotVec_eta = makeRandomNonUniformKnotVector(el_eta_start, p_eta_start, 2, true);
end

np_xi  = length(knotVec_xi)  - p_xi_start  - 1;                 % Number of control points in xi direction
np_eta = length(knotVec_eta) - p_eta_start - 1;                 % Number of control points in eta direction

el_xi_start;
Px = (findGrevillePoints(knotVec_xi, p_xi_start)-knotVec_xi(1))* (domain.endX-domain.startX)/(knotVec_xi(end)-knotVec_xi(1)) + domain.startX; Px = repmat(Px',1,np_eta);
Py = (findGrevillePoints(knotVec_eta, p_eta_start)-knotVec_eta(1))*(domain.endY-domain.startY)/(knotVec_eta(end)-knotVec_eta(1)) + domain.startY; Py = repmat(Py,np_xi,1);

switch shape
    case 'parall-x'
        Px = Px+Py;
        
    case 'parall-y'
        Py = Px+Py;
        
        n_y = Px(end,1) - Px(1,1);
        n_x = Py(1,end) - Py(end,end);
        
        h_neu{3} = @(x,y) -(dudx_exact(x,y)*n_x + dudy_exact(x,y)*n_y);
        h_neu{4} = @(x,y)   dudx_exact(x,y)*n_x + dudy_exact(x,y)*n_y;
    case 'plateHole'
        
        Px = [0,1,2; 0,1,2; 0, 9/4, 2.5; 3,3,3];
        Py = [0,0,0; 3, 3/4, 0.5; 3,2,1; 3,2,1];
        p_xi_start = 2; p_eta_start = 2; nopr = 1;
        
        knotVec_xi  = makeUniformKnotVector(p_xi_start,2);
        knotVec_eta = makeUniformKnotVector(p_eta_start,1);
    case 'L-shape'
        
         p_xi = 2; p_eta = 2;
         nopr = 1;
         knotVec_xi = [0 0 0 0.5 0.5 1 1 1];
         knotVec_eta= [0 0 0 1 1 1];
         Px = [1,1,1; 0,0,0.5; -1,-0.5,0; -1,-0.55, 0; -1,-0.6,0];
         Py = [-1,-0.6,0; -1, -0.55,0; -1,-0.5,0; 0,0,0.5; 1,1,1];
         el_xi_start = countElements(knotVec_xi); el_eta_start = countElements(knotVec_eta);
    case 'rect-8'
        
        p_xi_start = 1; p_eta_start = 1;
        nopr = 7; % number of p-refinement
        
end

p_xi = p_xi_start;
p_eta = p_eta_start;

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
        
        
        hr
        disp('------------------------------------')
        %el_eta = el_xi; %el_xi;
        
        %[knotVec_xi, knotVec_eta, Px, Py] = h_refinement(knotVec_xi, knotVec_eta ,p_xi, p_eta, Px, Py);
        
        el_xi = countElements(knotVec_xi);
        el_eta = countElements(knotVec_eta);
        
        np_xi  = length(knotVec_xi)  - p_xi  - 1;                 % Number of control points in xi direction
        np_eta = length(knotVec_eta) - p_eta - 1;                 % Number of control points in eta direction
        
        
        
        h(p_xi, hr) = 1/el_xi;%max(findMaxStepSize(knotVec_xi), findMaxStepSize(knotVec_eta));
        
        
        
        
        
        
        
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
        
        
        
        [errorH0(p_xi,hr), errorE(p_xi,hr), uInt] = evaluateU(Uw, U_exact, knotVec_xi, knotVec_eta, p_xi, p_eta, el_xi, el_eta, np_xi, np_eta, Px, Py, f, DIR, ID, IEN, LM, h_neu, NEU, dudx_exact, dudy_exact, domain);
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



nx = 100; ny = 100;
xi =  linspace(knotVec_xi(1) ,knotVec_xi(end) ,nx);
eta = linspace(knotVec_eta(1),knotVec_eta(end),ny);

N_xi = BsplineBasis(knotVec_xi, p_xi, xi);
N_eta = BsplineBasis(knotVec_eta, p_eta, eta);


elements_xi = findElements(knotVec_xi);
elements_eta = findElements(knotVec_eta);

elements_xi = [knotVec_xi(elements_xi), knotVec_xi(end)];
elements_eta = [knotVec_eta(elements_eta), knotVec_eta(end)];

N_xi_k = BsplineBasis(knotVec_xi, p_xi, elements_xi);
N_eta_k = BsplineBasis(knotVec_eta, p_eta, elements_eta);

figure(88) 
colors = get(gca, 'ColorOrder');

Cx = N_xi'*Px*N_eta;
Cy = N_xi'*Py*N_eta;
hold on
surf(Cx,Cy, zeros(size(Cx)), 'linestyle', 'none', 'FaceColor', colors(1,:), 'FaceAlpha',0.5 )
plot(Px,Py,'k')
plot(Px',Py','k')
Cx_k = (N_xi_k'*Px*N_eta_k)
Cy_k = (N_xi_k'*Py*N_eta_k)


% plot(Cx_k ,Cy_k , 'color', colors(1,:), 'linewidth',3)
plot(Cx_k',Cy_k', 'color', colors(1,:), 'linewidth',3)
plot(Px,Py,'o','MarkerEdgeColor', 'none','MarkerFaceColor', 'k', 'MarkerSize', 10)
plot(N_xi_k'*Px*N_eta_k,N_xi_k'*Py*N_eta_k, 'o','MarkerEdgeColor', 'none','MarkerFaceColor', colors(1,:), 'MarkerSize', 6)
set(gca,'fontsize', 18);


axis([min(min(Px)),max(max(Px)),min(min(Py)),max(max(Py))])

if isErrorPlot
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %       H^0 error
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    figure(4)
    
    loglog(h(2,:)',errorH0(2,:)','linewidth',3, 'color', colors(2,:))
    hold on
    grid on
    
    set(gca, 'ColorOrderIndex', 1)
    set(gca,'fontsize', 16);
    
    leg = legend('u_h, p=1','u_h, p=2', 'u_h, p=3', 'C_1h^2 (reference)', 'C_2h^3 (reference)', 'C_3h^4 (reference)'); % 'u=C_3x³')
    set(leg,'FontSize',11);
    xlabel('h', 'FontWeight','bold')
    ylabel('||u-u_h||_{H^0}', 'FontWeight','bold')
    factor = zeros(1,p_xi)
    legendString = cell(1, p_xi);
    for hi = 2:p_xi
        set(gca, 'ColorOrderIndex', hi)
        factor(hi) = (errorH0(hi,:)/h(hi,:).^(hi+1))*diff;
        loglog(h(hi,:)', factor(hi)*h(hi,:)'.^(hi+1), '--')
        legendString{hi} = ['u_h, p=', num2str(hi)];
        t = text(h(hi,end)/1.5,factor(hi)*h(hi,end)'.^(hi+1),['h^', num2str(hi+1)]);
        t.FontSize = 14;
        t.Color = colors(mod(hi-1,p_xi)+1,:);
        
    end
    %legendString{p_xi} = ['u_h, p=', num2str(p_xi)];
    leg = legend(legendString{2:end});
    set(leg,'FontSize',11);
    %t = text(h(end,end)/1.5,errorH0(end,end)/14', 'Machine error');
    %t.FontSize = 14;
    %t.Color = colors(mod(p_xi-1,p_xi)+1,:);

    
    
    
    
    
    
    
    
    
    
    if strcmp(shape,'L-shape')
        close
        diff = 0.3;
        figure(4)
        fig = gca;
        set(fig, 'ColorOrderIndex', 2)
        colors = get(gca, 'ColorOrder');
        loglog(h(2,:)',errorH0(2,:)','linewidth',3,'color',colors(2,:))
        hold on
        fig = gca;
        grid on
        factor22 = (errorH0(2,:)/h(2,:).^(2))*diff;
        factor23 = (errorH0(2,:)/h(2,:).^(3))*diff;
        set(fig, 'ColorOrderIndex', 1)
        loglog(h(2,:)',factor22*(h(2,:)'.^(2)),'--')
        loglog(h(2,:)',factor23*(h(2,:)'.^(3)),'--')
        leg = legend('u_h, p=2', 'C_1h^2 (reference)', 'C_2h^3 (reference)');
        set(leg,'FontSize',11);
        set(gca,'fontsize', 18);
        xlabel('h', 'FontWeight','bold')
        ylabel('||u-u_h||_{H^0}', 'FontWeight','bold')
    end
    

    % title('H^0,   u = x^3')
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
    
    loglog(h(2,:)',errorE(2,:)','linewidth',3, 'color', colors(2,:))
    hold on
    grid on
    
    set(gca,'fontsize', 16);
    set(gca, 'ColorOrderIndex', 1)
    
    factor = zeros(1,p_xi);
    legendString = cell(1,p_xi);
    for hi = 2:p_xi
        set(gca, 'ColorOrderIndex', hi)
        if  sum(errorE(hi,:)) ~=0
            factor(hi) = errorE(hi,:)*diff/h(hi,:).^(hi);
            loglog(h(hi,:)',factor(hi)*(h(hi,:)'.^(hi)),'--')
            legendString{hi} = ['u_h, p=', num2str(hi)];
            
            t = text(h(hi,end)/1.5,factor(hi)*(h(hi,end)'.^(hi)),['h^', num2str(hi),'     ']);
            t.FontSize = 14;
            t.Color = colors(mod(hi-1,p_xi)+1,:);
        end
    end
    %legendString{p_xi} = ['u_h, p=', num2str(p_xi)];
    leg = legend(legendString{2:end});
    set(leg,'FontSize',11);
    
    %t = text(h(end,end)/1.5,errorE(end,end)/14', 'Machine error');
    %t.FontSize = 14;
    %t.Color = colors(mod(p_xi-1,p_xi)+1,:);

    
    xlabel('h', 'FontWeight','bold')
    ylabel('||u-u_h||_{E}', 'FontWeight','bold')

    if strcmp(shape,'L-shape')
        close 
        figure(6)
        colors = get(gca, 'ColorOrder');
        loglog(h(2,:)',errorE(2,:)','linewidth',3, 'color',colors(2,:))
        hold on
        grid on
        set(gca,'fontsize', 18);
        loglog(h(2,:)',factor22*(h(2,:)'.^(2)),'--', 'color',colors(2,:))
        
        legend('u_h, p=2', 'C_1h^2 (reference)')
        xlabel('h', 'FontWeight','bold')
        ylabel('||u-u_h||_{E}', 'FontWeight','bold')
    end
    
    % ylabel('Energy error norm:   ||u-u_h||_{E}')
    % title('Energy error norm,   u = x^3')
    %axis([max(min(min(h))-max(max(h))/2,0), max(max(max(h))+max(max(h))/2,0), max(min(min(errorE))-max(max(errorE))/2,0), max(max(max(errorE))+max(max(errorE))/2,0)])
end

