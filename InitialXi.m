function [xi0] = InitialXi(Minv,S,p,dt)
xi0 = zeros(size(p,2),1);
% Generating the first two vectors of the xi matrix. Which are needed for
% the recursive temporal evolution of xi(t).
% The entries are: xi(1) = xi(t_-1) = u_0(:,x_i) and xi(:,2) = xi_0

f0 = SourceFunction(p,0);

[u0,v0] = InitialConditions(p);

a0 = Minv*(f0 -S*u0);


xi0(:,1) = u0 - dt*v0 + dt^2/2*a0;

xi0(:,2) = u0;

end

function [u0,v0] = InitialConditions(x)

u0 = cos( 5*pi*sqrt( sum(x.^2,1) ) ) ./ ( 1 + 10* sqrt( sum(x.^2,1) ) );
v0 = zeros(size(u0));
u0 = u0';
v0 = v0';
end