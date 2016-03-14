function hmax = findMaxRadius(Px,Py)
el_x = size(Px,2)-1;
el_y = size(Px,1)-1;
hmax = 0;
for i = 1:el_x
    for j = 1:el_y
        h = longestEdge(Px(j,i), Py(j,i), Px(j+1,i), Py(j+1,i), ...
            Px(j+1,i+1), Py(j+1,i+1), Px(j,i+1), Py(j,i+1));
        if hmax < h
            hmax = h;
        end
    end
end


end

function h = longestEdge(ax, ay, bx, by, cx, cy, dx, dy)
    ab = vectorLength(ax,ay,bx,by);
    bc = vectorLength(bx,by,cx,cy);
    cd = vectorLength(cx,cy,dx,dy);
    da = vectorLength(dx,dy,ax,ay);
    
    h = max([ab,bc,cd,da]);
    
end

function L = vectorLength(ax,ay,bx,by)
    L = sqrt((bx-ax)^2 + (by-ay)^2);
end