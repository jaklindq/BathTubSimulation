function [] = plotTub(xi,p,hmax)
x = p(1,:);
y = p(2,:);
x_edge=[floor(min(x)):hmax:ceil(max(x))];
y_edge=[floor(min(y)):hmax:ceil(max(y))];
[X,Y]=meshgrid(x_edge,y_edge);
Z=griddata(x,y,xi,X,Y);
surf(X,Y,Z)
az = 10;
el = 10;
view(az, el);
axis equal
axis([-0.5 2.5 -0.5 1.5 -0.2 1.2])

end

