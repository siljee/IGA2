function plotBspline(x, N)
% This function is made to plot B-spline basis functions with a set of 
% plot properties that the developer found pretty. 

fig = figure(33);
plot(x, N, 'LineWidth', 2);                                     % Plot and set thickness of the line
setPlotProperties(x, fig)
end

function setPlotProperties(x, fig)
% Specify the plot properties

% Constants
WINDOW_WIDTH = 600;
WINDOW_HEIGHT = WINDOW_WIDTH * 3 / x(end);

% Properties
fig.Position = [1000 1000 WINDOW_WIDTH WINDOW_HEIGHT];  % Determines the size of the window
axis([x(1) x(end) -0.1 1.1]);                           % Set the span of the axes
ax = gca;
ax.FontSize = 12;
ax.XTick = (x(1) :1: x(end));                           % Set the stepsize of axes values
ax.YTick = (0 :1/4: 1);

% title('B-spline Basis Functions');                    % Set title
grid on;                                                % Set square grid in background
end