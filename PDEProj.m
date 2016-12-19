cd ~/Chalmers/TMA690_PDE/
% addpath ~/Script/matlab2tikz/src/
clear all
clc
% TODO:
% Fix randi to rand in load vector
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%% PART I %%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% TESSELATION
hmax = 0.05;
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

% tubGeom =[  4 %Specify rectangle
%     1
%     0.5 % x-coord of vertices
%     1
%     0.5
%     0
%     0
%     ];

[g,bt] = decsg(tubGeom);
[p,e,t] = initmesh(g,'hmax',hmax);
% figure(1)
% clf
% pdemesh(p,e,t)
% axis([-.5 2.5 -.5 1.5])
% xlabel('x')
% ylabel('y')
% legend('Discretisation')
% axis equal

%
% Given the mesh; generate matrices M, S, and vector F.
m = size(p,2);
dt = 0.001;
M = MassMatrix(p,t);
Minv = inv(M);
S = StiffnessMatrix(p,t);
F = zeros(m,1);
xi = InitialXi(Minv,S,p,dt);

A = dt^2*Minv*F; %Ok onlu do this once.
B = dt^2*Minv*S;


fig = figure(2);
tid = 0;
speedParam = 8;
while ishandle(fig)
    tid = tid + 1;
    %     xi =[xi (2*xi(:,end)-xi(:,end-1)+dt^2*Minv*(F-S*xi(:,end)))];
    xi =[xi (2*xi(:,end)-xi(:,end-1)+A-B*xi(:,end))];
    
    xi = xi(:,end-2:end); % Crop xi to reduce memory load; Comment out
    % if you wish to obtain whole solution.
    % Plot each speedParam time step.
    if mod(tid,speedParam) == 0
        
        pdeplot(p, e, t, 'zdata', xi(:,end), 'xydata', xi(:, end),...
            'mesh','on','colorbar','off');
        set(gca,'ZLim',[-1 1])
        drawnow
        
    end
end
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%% PART II %%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% TESSELATION
clear all
hmax = 0.025;
tubGeom =[  3 %Specify rectangle
    4
    0 % x-coord of vertices
    1
    1
    0
    0 % y-coord of vertices
    0
    1
    1];

% tubGeom =[  4 %Specify ellipse
%     1
%     0.5
%     1
%     0.5
%     0
%     0
%     ];

[g,bt] = decsg(tubGeom);
[p,e,t] = initmesh(g,'hmax',hmax);
% figure(1)
% clf
% pdemesh(p,e,t)
% axis([-.5 2.5 -.5 1.5])
% xlabel('x')
% ylabel('y')
% legend('Discretisation')
% axis equal

%
% Given the mesh; generate matrices M, S, and vector F.

dt = 0.001;
% hmax/dt
M = MassMatrix(p,t);
Minv = inv(M);
S = StiffnessMatrix(p,t);
f0 = SourceFunction(p,0);
%F = LoadVector(p,t,0,pi);
xi = InitialXi(Minv,S,p,dt);
%
B = dt^2*Minv*S;

fig = figure(2);
tid = 0;

%global az; global elev;
az = 85;
elev = 16.2;
sld1 = uicontrol('Style', 'slider',...
    'Min',0,'Max',180,'Value',az,...
    'Position', [400 20 120 20],...
    'Callback', @surfazimuth);

sld2 = uicontrol('Style', 'slider',...
    'Min',0,'Max',180,'Value',elev,...
    'Position', [400 40 120 20],...
    'Callback', @surfelev);

speedParam = 1;
pause(2)

while ishandle(fig)
    
    F = LoadVector(p,t,tid);
    %     xi =[xi (2*xi(:,end)-xi(:,end-1)+dt^2*Minv*(F-S*xi(:,end)))];
    xi =[xi (2*xi(:,end)-xi(:,end-1)+dt^2*Minv*F-B*xi(:,end))];
    
    xi = xi(:,end-2:end); % Crop xi to reduce memory load; Comment out
    % if you wish to obtain whole solution.
    
    % Plot each speedParam time step.
    % if mod(plotTicker,speedParam) == 0
    pdeplot(p, e, t, 'zdata', xi(:,end), 'xydata', xi(:, end),...
            'mesh','on','colorbar','off');
        set(gca,'ZLim',[-1 1])
        drawnow
    % end
    tid = tid+1
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% NON TRIVIAL SOURCE TERM f
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cd ~/Chalmers/TMA690_PDE/

clear all
hmax = 0.05;
tubGeom =[  3 %Specify rectangle
    4
    0 % x-coord of vertices
    1
    1
    0
    0 % y-coord of vertices
    0
    1
    1];

% tubGeom =[  4 %Specify ellipse
%     1
%     0.5
%     1
%     0.5
%     0
%     0
%     ];

% Generate geometry for the rectangle (pre-built functions)
[g,bt] = decsg(tubGeom);
[p,e,t] = initmesh(g,'hmax',hmax);


% Parameters and constants
dt = 0.001;
T = 2e3;
m = size(p, 2);
% Place drop in centre
% dropInd = randi(m)
[~,dropInd] = min(sum(([0.5; 0.5] - p).^2,1));
dropTime = 150;
f = -20 % Value of source term (Diracesque)

% Matrices and and preallocation
xi = zeros(m, T);
F = LoadVector(p,t,hmax,f,dropInd,dropTime,T);
S = StiffnessMatrix(p, t);       
M = MassMatrix(p, t);


% Iterate for xi(t+1)
for tIter = 3:T
    b = dt^2 * F(:, tIter) + M*(2*xi(:, tIter - 1) - xi(:, tIter-2)) - dt^2 * S * xi(:, tIter-1);
    xi(:, tIter) = M\b;
end

% Plotting
fig = figure(2);
clf
    for tIter = 1:2:T
        pdeplot(p, e, t, 'zdata', xi(:,tIter), 'xydata', p,...
            'mesh','on','colorbar','off');
        set(gca,'ZLim',[-1 1])
        str = strcat('Time until drop: ', num2str((dropTime-tIter)*(dropTime>tIter)));
        text(0.9,0.9,0.9,str)
        drawnow
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% DIRICHLET BOUNDARY CONDITIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cd ~/Chalmers/TMA690_PDE/

clear all
hmax = 0.05;
tubGeom =[  3 %Specify rectangle
    4
    0 % x-coord of vertices
    1
    1
    0
    0 % y-coord of vertices
    0
    1
    1];

% tubGeom =[  4 %Specify ellipse
%     1
%     0.5clear 
%     1
%     0.5
%     0
%     0
%     ];

% Generate geometry for the rectangle (pre-built functions)
[g,bt] = decsg(tubGeom);
[p,e,t] = initmesh(g,'hmax',hmax);


% Parameters and constants
dt = 0.001;
T = 2e3;
m = size(p, 2);


% Matrices and and preallocation
xi = zeros(m, T);
F = zeros(m, T);
S = StiffnessMatrix(p, t);       
M = MassMatrix(p, t);

% Remove boundary vertices from the problem. 
[F,S,M,BoundNodesCompl] = Dirichlet(F,S,M,e);

% Generate load vector corresponding to source term (Diracesque)
% dropInd = randi(m)
dropInd = 280
F(dropInd, 20) = -40;


% Iterate for xi(t+1)
for tIter = 3:T
    b = dt^2 * F(:, tIter) + M*(2*xi(BoundNodesCompl, tIter - 1) - xi(BoundNodesCompl, tIter-2)) - dt^2 * S * xi(BoundNodesCompl, tIter-1);
    xi(BoundNodesCompl, tIter) = M\b;
end

% xi = [zeros(length(BoundNodes),T); xi];
fig = figure(2);
clf
    for tIter = 1:T
        pdeplot(p, e, t, 'zdata', xi(:,tIter), 'xydata', p,...
            'mesh','on','colorbar','off');
        set(gca,'ZLim',[-1 1])
        drawnow
    end

%%
