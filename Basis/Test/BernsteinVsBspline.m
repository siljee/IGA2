addpath ../../; addpath ../../HelpFunctions/;  
close all
p = 2;
elements = 4;
n = elements*100;

%knotVec = makeUniformKnotVector(p,elements);
%knotVec = [0,0,0,1/4,1/4,0.5,0.5,3/4,3/4, 1,1,1]
%knotVec = [0,0,0,1/3,2/3,1,1,1];    
%P = [0,1,0,2,1,3]';
elements = 5;
knotVec = [0,0,0,1,2,3,4,4,5,5,5];
P = [1,0; 1,1; 2,0; 2,3; 0.5,1; 0,2; 1,3; 1,2];

x = linspace(0,knotVec(end),n);
grev = findGrevillePoints(knotVec,p);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  B-spline
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

B = BsplineBasis(knotVec, p,  x);




figure(1)
hold on
plot(x,B, 'linewidth', 2)
%plot(knotVec,zeros(size(knotVec)),'o')
set(gca,'fontsize',13)
plot(x,zeros(size(x)),'w','linewidth',2)
plot(x,zeros(size(x)),'k','linewidth',1)


C = B'*P;
figure(2)
hold on
plot(C(:,1),C(:,2), 'linewidth',2)
plot(P(:,1),P(:,2),'k', 'linewidth', 1)
h = plot(P(:,1),P(:,2), 'o', 'MarkerEdgeColor', 'none', 'MarkerSize',8);
set(h, 'MarkerFaceColor', get(h, 'Color'));
hold off
axis off


% if 0
% figure(2)
% hold on
% plot(x,B'*P,'linewidth',2)
% h = plot(grev,P,'-k');
% h = plot(grev,P,'o','MarkerEdgeColor', 'none', 'MarkerSize',8);
% %set(h, 'MarkerFaceColor', get(h, 'Color'));
% xlabel('x')
% ylabel('y')
% set(get(gca,'YLabel'),'Rotation',0); 
% set(gca,'fontsize',13)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  B-spline
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C = BezierExtractionOperator(knotVec, p);

xb = linspace(0,1,n/elements);
b = BernsteinBasis(p,xb);



for e =1:elements
    e
    xplot = x(n/elements*(e-1)+1:n/elements*e);
    if e ==5
        Pb = P(e+1:e+1+p,:);
    else
        Pb = P(e:e+p,:);
    end
    (C(:,:,e)'*Pb);
    grevP = 0:1/(length(knotVec)-p-1):1;
    grevP = 0:1/elements:1;
    
    figure(3)
    colors = get(gca, 'ColorOrder');
    hold on
    plot(xplot, C(:,:,e)*b, 'Color', colors(e,:))
    
    
    figure(4)
    axis off
    hold on
    colors = get(gca, 'ColorOrder');
    
    grevPb = [grevP(e), (grevP(e)+grevP(e+1))/2, grevP(e+1)];
    Curve = (C(:,:,e)*b)'*Pb
    
    plot((C(:,:,e)'*Pb(:,1)), (C(:,:,e)'*Pb(:,2)), '-k')
    h = plot((C(:,:,e)'*Pb(:,1)), (C(:,:,e)'*Pb(:,2)), 'o','Color', colors(e,:), 'MarkerSize', 8);
    set(h, 'MarkerFaceColor', get(h, 'Color'));

    plot(Curve(:,1),Curve(:,2),'Color', colors(e,:), 'linewidth', 2)
    
    if(e~= 1) 
        colors(e-1,:)
        colors(e,:)

        plot(huskx(p+1), husky(p+1), 'o','Color', colors(e-1,:), 'MarkerSize', 3, 'MarkerFaceColor', colors(e-1,:))
    end
    huskx = (C(:,:,e)'*Pb(:,1));
    husky = (C(:,:,e)'*Pb(:,2));
end
