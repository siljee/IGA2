n = 40;
xstart = 0;
xstop = 3;
p = 3;
x = linspace(xstart,xstop,n);
t = linspace(xstart,xstop,p+1);
points = [1.5,2,-1,1] %rand(1,p+1)*5-2;
L = ones(p+1,n);


for i = 1:p+1           % Bases
    for j = 1:p+1       % Parameters
        if j~= i
            L(i,:) = L(i,:).*(x-t(j))/(t(i)-t(j));
        end
    end
end

close all
figure(1)
hold on
grid on

colors = get(gca, 'ColorOrder');

plot(x,zeros(size(x)),'k')
plot(x,ones(size(x)),':k')
plot(x,L,'linewidth',2)

for basis = 1:p+1
    plot(t(basis),1,'o', 'MarkerEdgeColor', 'none','MarkerFaceColor', colors(basis,:), 'MarkerSize',8);
    plot(t(basis),0,'ok','MarkerFaceColor', 'auto', 'MarkerSize',7);
end

set(gca,'fontsize',16)

figure(2)
hold on
grid on
colors = get(gca, 'ColorOrder');

plot(x,zeros(size(x)),'k')
plot(x,L.*repmat(points',1,40),'linewidth',2)
plot(x,points*L, 'k','linewidth',3)

for basis = 1:p+1
    plot(t(basis),points(basis),'o', 'MarkerEdgeColor', 'none','MarkerFaceColor', colors(basis,:), 'MarkerSize',8);
    plot(t(basis),0,'ok','MarkerFaceColor', 'auto', 'MarkerSize',7);
end

set(gca,'fontsize',16)



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Samme figur i Bernstein
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



B = bernstein