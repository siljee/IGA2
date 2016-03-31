close all; clear; addpath ../
knotVector = [0,0,0,0,1,2,2,3,3,3,4,4,4,4];
p = 3;
x = linspace(0,4,1000);

N = BsplineBasis(knotVector, p, x);

fig = figure(1);
hold on
ax = gca
colors = get(ax, 'ColorOrder');
plot(x(1:125),N(:,1:125), 'color', colors(1,:),'LineWidth', 2)
for i = 0:2
    start = 125 + 250*i;
    slutt = 125 + 250*(i+1)
    plot(x(start:slutt),N(:,start:slutt), 'color', colors(i+2,:), 'LineWidth', 2)
end
plot(x(slutt:slutt+125),N(:,slutt:slutt+125), 'color', colors(5,:), 'LineWidth', 2)

ax.FontSize = 14;  