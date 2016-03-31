close all; addpath ../
knotVector = [0,0,0,0,1,2,2,3,3,3,4,4,4,4];
p = 3;
x = linspace(0,4,1000);

% Constants
WINDOW_WIDTH = 700;
WINDOW_HEIGHT = WINDOW_WIDTH * 2;% / knotVector(end);  
close all


N2 = BsplineBasis(knotVector, p-1,  x);
N4 = BsplineBasis(knotVector, p, x);

[dN1] = BsplineDerivative(p, 1, knotVector, x);
[dN2] = BsplineDerivative(p, 2, knotVector, x);
[dN3] = BsplineDerivative(p, 3, knotVector, x);


fig10 = figure(10);
subplot(4,1,1)
hold on

colors = get(gca, 'ColorOrder');
plot(x(1:125),N4(:,1:125), 'color', colors(1,:),'LineWidth', 2)
for i = 0:2
    start = 125 + 250*i;
    slutt = 125 + 250*(i+1)
    plot(x(start:slutt),N4(:,start:slutt), 'color', colors(i+2,:), 'LineWidth', 2)
end
plot(x(slutt:slutt+125),N4(:,slutt:slutt+125), 'color', colors(5,:), 'LineWidth', 2)

%title('original')
subplot(4,1,2)
hold on
%plot(x, dN1, 'k')
%plot(x, dN1([3,5,7,9,10],:), 'LineWidth', 2)

colors = get(gca, 'ColorOrder');
plot(x(1:125),dN1(:,1:125), 'color', colors(1,:),'LineWidth', 2)
for i = 0:2
    start = 125 + 250*i;
    slutt = 125 + 250*(i+1)
    plot(x(start:slutt),dN1(:,start:slutt), 'color', colors(i+2,:), 'LineWidth', 2)
end
plot(x(slutt:slutt+125),dN1(:,slutt:slutt+125), 'color', colors(5,:), 'LineWidth', 2)


subplot(4,1,3)
hold on
%plot(x,dN2, 'k')
%plot(x,dN2([3,5,7,9,10],:), 'LineWidth', 2)
plot(x(1:125),dN2(:,1:125), 'color', colors(1,:),'LineWidth', 2)
for i = 0:2
    start = 125 + 250*i;
    slutt = 125 + 250*(i+1)
    plot(x(start:slutt),dN2(:,start:slutt), 'color', colors(i+2,:), 'LineWidth', 2)
end
plot(x(slutt:slutt+125),dN2(:,slutt:slutt+125), 'color', colors(5,:), 'LineWidth', 2)



subplot(4,1,4)
plot(x,dN3, 'k')
hold on
%plot(x,dN3([3,5,7,9,10],:), 'LineWidth', 2)
% subplot(5,1,5)
% plot(x,dN4, 'k', 'LineWidth', 2)

plot(x(1:125),dN3(:,1:125), 'color', colors(1,:),'LineWidth', 2)
for i = 0:2
    start = 125 + 250*i;
    slutt = 125 + 250*(i+1)
    plot(x(start:slutt),dN3(:,start:slutt), 'color', colors(i+2,:), 'LineWidth', 2)
end
plot(x(slutt:slutt+125),dN3(:,slutt:slutt+125), 'color', colors(5,:), 'LineWidth', 2)


fig10.Position = [1000,1000, WINDOW_WIDTH, WINDOW_HEIGHT];

% 
% fig1.Position = [1000 1000 WINDOW_WIDTH WINDOW_HEIGHT];          % Determines the size of the window
% fig2.Position = [1000 1000 WINDOW_WIDTH WINDOW_HEIGHT];          % Determines the size of the window
% fig3.Position = [1000 1000 WINDOW_WIDTH WINDOW_HEIGHT];          % Determines the size of the window
% fig4.Position = [1000 1000 WINDOW_WIDTH WINDOW_HEIGHT];          % Determines the size of the window
% fig5.Position = [1000 1000 WINDOW_WIDTH WINDOW_HEIGHT];          % Determines the size of the window

if 0

% Figure 2.8 in Piegl1997tnb
knotVector = [0,0,0,0,1,3,4,6,7,8,8,8,8];
p = 3;  
[dN1, N, x] = derivateBsplineBasisFunctions(p, 1, knotVector, linspace(1,8,100));
N2 = BsplineBasisFunctions(knotVector, p-1, false, x );
N3 = BsplineBasisFunctions(knotVector, p, false, x);

hold off
fig = figure(100);

% Constants
WINDOW_WIDTH = 600;
WINDOW_HEIGHT = WINDOW_WIDTH * 3 / knotVector(end);  

plot(x,N2(5,:),'g--')
hold on
plot(x,N2(6,:),'r--')
plot(x,N3(5,:),'b')
plot(x,dN1(5,:),'m')


    
fig.Position = [1000 1000 WINDOW_WIDTH WINDOW_HEIGHT];          % Determines the size of the window
axis([knotVector(1) knotVector(end) -1 2]);                     % Set the span of the axes
ax = gca;                                                       
ax.FontSize = 12;                                           
ax.XTick = (knotVector(1) :1: knotVector(end));                 % Set the stepsize of axes values
ax.YTick = (-1 :1: 2);
    
title('B-spline Basis Functions');                              % Set title
grid on; 
hold off
end