addpath ../; addpath ../HelpFunctions/;

p_xi = 2; p_eta = 2;
knotVec_xi  = makeUniformKnotVector(p_xi,3);
knotVec_eta = makeUniformKnotVector(p_eta,2);
np_xi = length(knotVec_xi)-p_xi-1;
np_eta = length(knotVec_eta)-p_eta-1;

greville_xi = findGrevillePoints(knotVec_xi,p_xi);
greville_eta = findGrevillePoints(knotVec_eta,p_eta);

elementBndX = repmat(knotVec_xi(p_xi+1:end-p_xi),3,1);
elementBndY = repmat(knotVec_eta(p_eta+1:end-p_eta)',1,4);

pointsX = repmat(greville_xi',1,np_eta);
pointsY = repmat(greville_eta,np_xi,1);

ID = generateID_all4bnd(np_xi,np_eta,-1,-2,-3:-1:-3-np_xi+1,[])
IEN = generate_IEN(knotVec_xi, knotVec_eta, p_xi, p_eta)
LM = ID(IEN)

close all 
figure(1) % ID
axis off
hold on
colors = get(gca, 'ColorOrder');
% <<<<<<< HEAD
% plot(elementBndX,elementBndY,'k','Linewidth',3)
% plot(elementBndX',elementBndY','k','Linewidth',3)
% =======
plot(elementBndX,elementBndY,'color',colors(2,:),'Linewidth',3)
plot(elementBndX',elementBndY','color',colors(2,:),'Linewidth',3)
% >>>>>>> 2e08aeb7c56cc84de5fd300eedaedfa4910cb0bb

plot(pointsX,pointsY,'o','MarkerEdgeColor', 'k','MarkerFaceColor', colors(1,:), 'MarkerSize',9)

for i = 1: np_xi*np_eta
    t = text(pointsX(i),pointsY(i)-0.05,['   ' num2str(ID(i))]);
    t.FontSize = 16;
    t.FontWeight = 'bold';
    t.Color = colors(1,:);
    if ID(i)<0
        plot(pointsX(i),pointsY(i),'o','MarkerEdgeColor', 'k','MarkerFaceColor', colors(5,:)-0.188, 'MarkerSize',9)
        t.Color = colors(5,:)-0.188;
        if ID(i) < -2
            plot(pointsX(i),pointsY(i),'o','MarkerEdgeColor', 'k','MarkerFaceColor', colors(5,:), 'MarkerSize',9)
            t.Color = colors(5,:);
        end
    end
end

set(gca,'fontsize',16)

% while -3:-7 correspond to the interpolated values of g at the boundary. 
%
% On the top where there values are positive there is a natural boundary
% condition, like neumann or robin, the basis functions are included in the
% system of basis fuunctions


figure(2) % IEN
axis off
hold on
colors = get(gca, 'ColorOrder');
plot(elementBndX ,elementBndY ,'color',colors(2,:),'Linewidth',3)
plot(elementBndX',elementBndY','color',colors(2,:),'Linewidth',3)

plot(pointsX,pointsY,'o','MarkerEdgeColor', 'k','MarkerFaceColor', 'k', 'MarkerSize',9)

for i = 1: np_xi*np_eta
    t = text(pointsX(i),pointsY(i)-0.05,['   ' num2str(i)]);
    t.FontSize = 16;
    t.FontWeight = 'bold';
    t.Color = 'k';
%     if ID(i)<0
%         plot(pointsX(i),pointsY(i),'o','MarkerEdgeColor', 'k','MarkerFaceColor', colors(5,:)-0.188, 'MarkerSize',9)
%         t.Color = colors(5,:)-0.188;
%         if ID(i) < -2
%             plot(pointsX(i),pointsY(i),'o','MarkerEdgeColor', 'k','MarkerFaceColor', colors(5,:), 'MarkerSize',9)
%             t.Color = colors(5,:);
%         end
%     end
end

set(gca,'fontsize',16)

picOff = 0.1;
windowPos = [500 100 400 250];
saveFigures = 0;

fig = figure(11); % IEN 1 elemet
set(fig, 'Position', windowPos)
axis([0-picOff,1+picOff,-0-picOff,1+picOff])
axis off
hold on
colors = get(gca, 'ColorOrder');
outside = 0.05;
rectangle('Position',[greville_xi(1)-outside,greville_eta(1)-2*outside,greville_xi(3)-greville_xi(1)+3.5*outside,greville_eta(3)-greville_eta(1)+3*outside],'FaceColor',colors(6,:))

plot(elementBndX,elementBndY,'k','Linewidth',3)
plot(elementBndX',elementBndY','k','Linewidth',3)

plot(pointsX,pointsY,'o','MarkerEdgeColor', 'k','MarkerFaceColor', 'k', 'MarkerSize',15)

for i = 1: np_xi*np_eta
    t = text(pointsX(i),pointsY(i)-0.05,[' ' num2str(i)]);
    t.FontSize = 35;
    %t.FontWeight = 'bold';
    t.Color = 'k';
%     if ID(i)<0
%         plot(pointsX(i),pointsY(i),'o','MarkerEdgeColor', 'k','MarkerFaceColor', colors(5,:)-0.188, 'MarkerSize',9)
%         t.Color = colors(5,:)-0.188;
%         if ID(i) < -2
%             plot(pointsX(i),pointsY(i),'o','MarkerEdgeColor', 'k','MarkerFaceColor', colors(5,:), 'MarkerSize',9)
%             t.Color = colors(5,:);
%         end
%     end
end

set(gca,'fontsize',16)

if saveFigures
    print('Figres/IEN11','-dpng')
end


fig = figure(12); % IEN 2. elemet
set(fig, 'Position', windowPos)
axis([0-picOff,1+picOff,-0-picOff,1+picOff])
axis off
hold on
colors = get(gca, 'ColorOrder');
outside = 0.05;
rectangle('Position',[greville_xi(2)-outside,greville_eta(1)-2*outside,greville_xi(4)-greville_xi(2)+3.5*outside,greville_eta(3)-greville_eta(1)+3*outside],'FaceColor',colors(6,:))

plot(elementBndX,elementBndY,'k','Linewidth',3)
plot(elementBndX',elementBndY','k','Linewidth',3)

plot(pointsX,pointsY,'o','MarkerEdgeColor', 'k','MarkerFaceColor', 'k', 'MarkerSize',15)

for i = 1: np_xi*np_eta
    t = text(pointsX(i),pointsY(i)-0.05,[' ' num2str(i)]);
    t.FontSize = 35;
    %t.FontWeight = 'bold';
    t.Color = 'k';
%     if ID(i)<0
%         plot(pointsX(i),pointsY(i),'o','MarkerEdgeColor', 'k','MarkerFaceColor', colors(5,:)-0.188, 'MarkerSize',9)
%         t.Color = colors(5,:)-0.188;
%         if ID(i) < -2
%             plot(pointsX(i),pointsY(i),'o','MarkerEdgeColor', 'k','MarkerFaceColor', colors(5,:), 'MarkerSize',9)
%             t.Color = colors(5,:);
%         end
%     end
end

set(gca,'fontsize',16)

if saveFigures
    print('Figres/IEN12','-dpng')
end


fig = figure(13); % IEN 2. elemet
set(fig, 'Position', windowPos)
axis([0-picOff,1+picOff*1.5,-0-picOff,1+picOff])
axis off
hold on
colors = get(gca, 'ColorOrder');
outside = 0.05;
rectangle('Position',[greville_xi(3)-outside,greville_eta(1)-2*outside,greville_xi(5)-greville_xi(3)+3.8*outside,greville_eta(3)-greville_eta(1)+3*outside],'FaceColor',colors(6,:))

plot(elementBndX,elementBndY,'k','Linewidth',3)
plot(elementBndX',elementBndY','k','Linewidth',3)

plot(pointsX,pointsY,'o','MarkerEdgeColor', 'k','MarkerFaceColor', 'k', 'MarkerSize',15)

for i = 1: np_xi*np_eta
    t = text(pointsX(i),pointsY(i)-0.05,[' ' num2str(i)]);
    t.FontSize = 35;
    %t.FontWeight = 'bold';
    t.Color = 'k';
%     if ID(i)<0
%         plot(pointsX(i),pointsY(i),'o','MarkerEdgeColor', 'k','MarkerFaceColor', colors(5,:)-0.188, 'MarkerSize',9)
%         t.Color = colors(5,:)-0.188;
%         if ID(i) < -2
%             plot(pointsX(i),pointsY(i),'o','MarkerEdgeColor', 'k','MarkerFaceColor', colors(5,:), 'MarkerSize',9)
%             t.Color = colors(5,:);
%         end
%     end
end

set(gca,'fontsize',16)

if saveFigures
    print('Figres/IEN13','-dpng')
end


fig = figure(21); % IEN 1 elemet
set(fig, 'Position', windowPos)
axis([0-picOff,1+picOff,-0-picOff,1+picOff])
axis off
hold on
colors = get(gca, 'ColorOrder');
outside = 0.05;
rectangle('Position',[greville_xi(1)-outside,greville_eta(2)-2*outside,greville_xi(3)-greville_xi(1)+3.5*outside,greville_eta(4)-greville_eta(2)+3*outside],'FaceColor',colors(6,:))

plot(elementBndX,elementBndY,'k','Linewidth',3)
plot(elementBndX',elementBndY','k','Linewidth',3)

plot(pointsX,pointsY,'o','MarkerEdgeColor', 'k','MarkerFaceColor', 'k', 'MarkerSize',15)

for i = 1: np_xi*np_eta
    t = text(pointsX(i),pointsY(i)-0.05,[' ' num2str(i)]);
    t.FontSize = 35;
    %t.FontWeight = 'bold';
    t.Color = 'k';
%     if ID(i)<0
%         plot(pointsX(i),pointsY(i),'o','MarkerEdgeColor', 'k','MarkerFaceColor', colors(5,:)-0.188, 'MarkerSize',9)
%         t.Color = colors(5,:)-0.188;
%         if ID(i) < -2
%             plot(pointsX(i),pointsY(i),'o','MarkerEdgeColor', 'k','MarkerFaceColor', colors(5,:), 'MarkerSize',9)
%             t.Color = colors(5,:);
%         end
%     end
end

set(gca,'fontsize',16)

if saveFigures
    print('Figres/IEN21','-dpng')
end



fig = figure(22) % IEN 2. elemet
set(fig, 'Position', windowPos)
axis([0-picOff,1+picOff,-0-picOff,1+picOff])
axis off
hold on
colors = get(gca, 'ColorOrder');
outside = 0.05;
rectangle('Position',[greville_xi(2)-outside,greville_eta(2)-2*outside,greville_xi(4)-greville_xi(2)+3.5*outside,greville_eta(4)-greville_eta(2)+3*outside],'FaceColor',colors(6,:))

plot(elementBndX,elementBndY,'k','Linewidth',3)
plot(elementBndX',elementBndY','k','Linewidth',3)

plot(pointsX,pointsY,'o','MarkerEdgeColor', 'k','MarkerFaceColor', 'k', 'MarkerSize',15)

for i = 1: np_xi*np_eta
    t = text(pointsX(i),pointsY(i)-0.05,[' ' num2str(i)]);
    t.FontSize = 35;
    %t.FontWeight = 'bold';
    t.Color = 'k';
%     if ID(i)<0
%         plot(pointsX(i),pointsY(i),'o','MarkerEdgeColor', 'k','MarkerFaceColor', colors(5,:)-0.188, 'MarkerSize',9)
%         t.Color = colors(5,:)-0.188;
%         if ID(i) < -2
%             plot(pointsX(i),pointsY(i),'o','MarkerEdgeColor', 'k','MarkerFaceColor', colors(5,:), 'MarkerSize',9)
%             t.Color = colors(5,:);
%         end
%     end
end

set(gca,'fontsize',16)

if saveFigures
    print('Figres/IEN22','-dpng')
end


fig = figure(23); % IEN 2. elemet
set(fig, 'Position', windowPos)
axis([0-picOff,1+picOff*1.5,-0-picOff,1+picOff])
axis off
hold on
colors = get(gca, 'ColorOrder');
outside = 0.05;
rectangle('Position',[greville_xi(3)-outside,greville_eta(2)-2*outside,greville_xi(5)-greville_xi(3)+3.8*outside,greville_eta(4)-greville_eta(2)+3*outside],'FaceColor',colors(6,:))

plot(elementBndX,elementBndY,'k','Linewidth',3)
plot(elementBndX',elementBndY','k','Linewidth',3)

plot(pointsX,pointsY,'o','MarkerEdgeColor', 'k','MarkerFaceColor', 'k', 'MarkerSize',15)

for i = 1: np_xi*np_eta
    t = text(pointsX(i),pointsY(i)-0.05,[' ' num2str(i)]);
    t.FontSize = 35;
    %t.FontWeight = 'bold';
    t.Color = 'k';
%     if ID(i)<0
%         plot(pointsX(i),pointsY(i),'o','MarkerEdgeColor', 'k','MarkerFaceColor', colors(5,:)-0.188, 'MarkerSize',9)
%         t.Color = colors(5,:)-0.188;
%         if ID(i) < -2
%             plot(pointsX(i),pointsY(i),'o','MarkerEdgeColor', 'k','MarkerFaceColor', colors(5,:), 'MarkerSize',9)
%             t.Color = colors(5,:);
%         end
%     end
end

set(gca,'fontsize',16)

if saveFigures
    print('Figres/IEN23','-dpng')
end