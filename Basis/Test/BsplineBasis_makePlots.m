close all
addpath ../
%knotVec = [0,0,0,0,1,2,3,4,5,6,6,6,6]%4,4,4,4]
knotVec = [1,2,3,4,5,6,7,8,9,10];
x = linspace(knotVec(1),knotVec(7), 1000);
n = size(knotVec);
m = 3;

for p = 0:3
    figure(p+1)
    B = BsplineBasis(knotVec, p, x);
    for i = 1:m%4:4+2%(4 : size(B1,1)-3)
    B_i = B(i,:);
   
    subplot(3, 1, i)
    plot(x, B_i ,'LineWidth', 2);                                % Plot and set thickness of the line
    set(gca,'fontsize',16)
    grid on
   % setPlotProperties(knotVector, fig, ['Basis ', num2str(i)]);
    end
end
