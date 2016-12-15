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

%%
% Given the mesh; generate matrices M, S, and vector F.

dt = 0.002;
hmax/dt
M = MassMatrix(p,t);
Minv = inv(M);
S = StiffnessMatrix(p,t);
F = LoadVector(p,t);
xi = InitialXi(Minv,S,p,dt);

A = dt^2*Minv*F;
B = dt^2*Minv*S;


fig = figure(2);
tid = 0;

global az; global elev;
az = 45;
elev = 45;
speedParam = 6;
sld1 = uicontrol('Style', 'slider',...
        'Min',0,'Max',180,'Value',az,...
        'Position', [400 20 120 20],...
        'Callback', @surfazimuth);
    
sld2 = uicontrol('Style', 'slider',...
        'Min',0,'Max',180,'Value',elev,...
        'Position', [400 40 120 20],...
        'Callback', @surfelev);
    
while ishandle(fig)
    tid = tid + 1;
%     xi =[xi (2*xi(:,end)-xi(:,end-1)+dt^2*Minv*(F-S*xi(:,end)))];
    xi =[xi (2*xi(:,end)-xi(:,end-1)+A-B*xi(:,end))];

    xi = xi(:,end-2:end); % Crop xi to reduce memory load; Comment out
    % if you wish to obtain whole solution.
    % Plot each speedParam time step.
    if mod(tid,speedParam) == 0
        
        plotTub(xi(:,end),p,hmax,az,elev)
        
        drawnow
    end
end
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%% PART II %%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% TESSELATION
hmax = 0.04;
tubGeom =[  3 %Specify rectangle
    4
    0 % x-coord of vertices
    3
    3
    0
    0 % y-coord of vertices
    0
    2
    2];

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

dt = 0.0001;
% hmax/dt
M = MassMatrix(p,t);
Minv = inv(M);
S = StiffnessMatrix(p,t);
f0 = SourceFunction(p,0,pi,1);
% F = LoadVector(p,t,0,pi);
xi = InitialXi(Minv,S,p,dt);
%
B = dt^2*Minv*S;

fig = figure(2);
tid = 0;

global az; global elev;
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
    
plotTicker = 0;
speedParam = 1;

while ishandle(fig)

    tid = tid + dt;
    plotTicker = plotTicker +1 ;
    F = LoadVector(p,t,tid,dt);
%     xi =[xi (2*xi(:,end)-xi(:,end-1)+dt^2*Minv*(F-S*xi(:,end)))];
    xi =[xi (2*xi(:,end)-xi(:,end-1)+Minv*F-B*xi(:,end))];

    xi = xi(:,end-2:end); % Crop xi to reduce memory load; Comment out
    % if you wish to obtain whole solution.
    % Plot each speedParam time step.
%     if mod(plotTicker,speedParam) == 0
        plotTub(xi(:,end),p,hmax,az,elev)        
        drawnow
%     end
%     pause(0.5)
    
end

%%
tid = 0;
dt = 0.005
for i=1:100
    tid = [tid tid(end)+dt];
end
tid = tid';
tid./dt

