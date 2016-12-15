function [] = plotTub(xi,p,hmax,az,elev)

x = p(1,:);
y = p(2,:);
x_edge=[floor(min(x)):hmax:ceil(max(x))];
y_edge=[floor(min(y)):hmax:ceil(max(y))];
[X,Y]=meshgrid(x_edge,y_edge);
Z=griddata(x,y,xi,X,Y);
surf(X,Y,Z)
view(az, elev);
axis equal
axis([-0.1 2.1 -0.05 1.05 -0.7 1])

end

