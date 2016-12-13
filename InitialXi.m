function [xi0] = InitialXi(u0i,v0i,Minv,S,f0i,dt)

% Generating the first two vectors of the xi matrix. Which are needed for
% the recursive temporal evolution of xi(t).
% The entries are: xi(1) = xi(t_-1) = u_0(:,x_i) and xi(:,2) = xi_0

a0i = Minv(f0i -S*u0i);

xi0(:,1) = u0i - dt*v0i + dt^2/2*a0i;

xi0(:,2) = u0i;

end

