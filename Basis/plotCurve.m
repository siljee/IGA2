function plotCurve(C, P)
% A simple help-function to plot curves and control points

figure()
hold on
plot(C(1,:),C(2,:), 'linewidth',2)
plot(P(1,:),P(2,:), 'linewidth', 2)
plot(P(1,:),P(2,:), '*m', 'linewidth',2)
hold off
end