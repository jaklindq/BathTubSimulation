cd ~/Chalmers/TMA690_PDE/
% addpath ~/Script/matlab2tikz/src/
% Questions:
% Enumeration: 1D or 2D?
% How linearise xi
% Use center of gravity approximation for L(v)
% How to use [p,e,t]
% What is the subdomain number alltså t(4,:)? 
% Is D ũ +grad(u) and B the same though in the boundary? Or is B the ũ integral?

%% TESSELATION
hmax = 0.5;
tubGeom =[  3 %Specify rectangle
            4 
            0 % x-coord of vertices
            2
            2
            0
            0 % y-coord of vertices
            0
            1
            1];

[g,bt] = decsg(tubGeom);
[p,e,t] = initmesh(g,'hmax',hmax);
figure(1)
clf
pdemesh(p,e,t)
axis([-.5 2.5 -.5 1.5])
xlabel('x')
ylabel('y')
legend('Discretisation')
axis equal


%%
% MyPoissonSolver1(p,e,t);
%%
% Given the mesh; generate matrices M, S, and vector F. 
% N = size(p,2); %Number of nodes
% E = size(t,2); %Number of elements.
% 
% Write functions for these:
% M = Mmatrix(p,e,t);
% S = Smatrix(p,e,t);
% F = Fvector(); 
% Minv = inv(M); % All matrices are static so it's more efficient to invert
% A = (2*M+dt^2S) %Perhaps overkill for efficiency.
% this once.
% Get xi(-1) and xi(0). xi(0) from initial conditions and xi(-1) from
% Taylor backwards.
% xi(:,1) = ; % xi(-1)
dt = 0.1;
M = MassMatrix(p,t);
Minv = inv(M);
S = StiffnessMatrix(p,t);
f0 = SourceFunction(p,0);
F = LoadVector(f0,p,t);
xi = InitialXi(Minv,S,p,dt);


fig = figure(2)
tid = 0;

while ishandle(fig)
tid = tid + 1
xi =[xi (2*xi(:,end)-xi(:,end-1)+dt^2*Minv\(F-S*xi(:,end)))];
clf

plotTub(xi(:,end-1),p,hmax) 
drawnow
end

% The following line of code is if you use JE's gridfit:
% Z=gridfit(x,y,z,x_edge,y_edge);surf(p(1,:),p(2,:),Z)
%%
% Solve recursively for x(t+1):
dt = 1 % Length of time step.



