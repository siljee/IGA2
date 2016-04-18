close all; addpath ../; addpath ../Basis/; addpath ../HelpFunctions/;
fontsize = 13; colorOrder = [5,2,1,3,4];

knotVec = [0,0,0,0,1,2,3,3,4,4,4,5,5,5,5];
el = countElements(knotVec);
p = 3;
x = linspace(knotVec(1), knotVec(end),100);

N = BsplineBasis(knotVec,p,x);
figure(1)
plot(x,N,'linewidth',2);

C = BezierExtractionOperator(knotVec,p);


xb = linspace(0,1,100);

B = BernsteinBasis(p,xb);

figure(2)
ax = gca;
set(ax,'fontsize',fontsize)
colors = get(ax, 'ColorOrder');

hold on
for e = 1:el
    plot(xb+(e-1), B, 'linewidth', 2,'color', colors(colorOrder(e),:));
end


figure(3)
ax = gca;
set(ax,'fontsize',fontsize)
colors = get(ax, 'ColorOrder');

hold on
for e = 1:el
    plot(xb+(e-1), C(:,:,e)*B, 'linewidth', 2,'color', colors(colorOrder(e),:));
end