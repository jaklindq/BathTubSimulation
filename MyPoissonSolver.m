function [U] = MyPoissonSolver(p, t, e)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%% Purpose:  Solves  - div (d grad u) = f
%%           in a domain  D,  
%%           with boundary conditions
%%           d n grad u = b ( g - u), where  n  is the outer
%%           normal to the boundary.
%%   
%%        p  2 x n  vector of node coordinates, that is  p(1, nn)
%%           and  p(2, nn)  are the  x1  and  x2  coordinates,
%%           respectively, of node number  nn. 
%%
%%        t  4 x m  matrix of element node numbers  t(1, el) = nn1,
%%           t(2, el) = nn2,  and  t(3, el) = nn3, 1 <= el <= m,
%%           pointing to corresponding node coordinates  p(:, nn1),
%%           p(:, nn2)  and  p(:, nn3),  and subdomain reference
%%           tags  t(4, el),  useful if data are given by different
%%           expressions in different parts of the domain.
%%
%%        e  7 x q  matrix where  e(:, bel)  holds the following
%%           information about bdry element  bel:  e(1:2, bel)  are
%%           the node numbers (pointers to  p) of the two (endpoint)
%%           nodes of the bdry element,  e(3:4, bel)  are the relative
%%           positions of these two nodes on their bdry segment, with
%%           reference number  e(5, bel),  and  e(6:7, bel) are the
%%           reference tags of the domains/triangle elements tags to
%%           the left and right, respectively, of the bdry element.
%%           The exterior region has reference number 0.
%%
%%        U  n x 1  matrix , the solution values
%%           at node points  p . 
%%           
%%
%% Example:  To solve  - div(d grad u) = f  in
%%           D = [0 1] x [0,1]  with  d = 1,  and
%%           f = sin(10 * x1), with
%%           boundary condition  u = 0  (corresponding to  b = oo and
%%           g = 0)  we give  a  and  f  as node functions in
%%           the code, and similarly  b  and  g.  To define the domain
%%           and mesh we may use the matlab function initmesh, by first
%%           defining the matrix  g = [c1 c2 c3 c4], with columns  ci
%%           defining the four edges/sides of  D. The first number in  ci
%%           should be  2  indicating a linear edge/side, followed by the
%%           two  x1  coordinates and then the two  x2  coordinates for
%%           the two endpoints of the edge, and finally domain reference
%%           numbers to the domains to the left and right, respectively,
%%           of the edge, where the outer domain has reference number 0.
%%           For example  c1 = [2 0 1 0 0 1 0]  defines the edge  (x1, 0),
%%           0 <= x2 <= 1  with  x1-end-coordinates  (0, 1)  and  x2-end-
%%           coordinates  (0, 0),  and giving the domain  D  reference
%%           number  1.  Once  g  is defined the call (note the order)
%%           [p, e, t] = initmesh(g)  gives  p,  t,  and  e.
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('Du måste lägga till en massa.')
disp('ElemBounMatrix, ElemLoadVector och ElemBounVector saknas!')
error('Programmet är inte färdigt!')
% initializing

  n = size(p, 2);
  x1 = p(1,:)';                     % list of x1 coordinates of all nodes
  x2 = p(2,:)';
  d = ones(n, 1);                   % conductivity
  D = DiffMatrix(p, t, d);
  b = 10000*ones(n, 1);             % bdry conductivity
  B = BounMatrix(p, e, b);
  f = sin(10*x1);                   % heat production, EXAMPLE


  F = LoadVector(p, t, f);
  g = x1.*x2.*(1-x2);               % exterior temperature, EXAMPLE
  G = BounVector(p, e, b.*g);
  
    
  A = ( D + B);
  
  
  b = F + G;
  U = A\b;
  
% visualization

  set(gcf, 'DoubleBuffer','on')
  shg
  pdeplot(p, e, t, 'zdata', U, 'xydata', U,...
                              'mesh','on','colorbar','off')
  set(gca,'ZLim',[-.5 .5])
  drawnow
  
  
  
% subroutines
  
% Assembling routines.  
  
function D = DiffMatrix(p, t, d)
  n = size(p, 2);
  D = zeros(n, n);
  for el = 1 : size(t, 2)
    dD = ElemDiffMatrix(p, t, d, el);
    nn = t(1:3, el);
    D(nn, nn) = D(nn, nn) + dD;
  end
  
  
function F = LoadVector(p, t, f)
  n = size(p, 2);
  F = zeros(n, 1);
  for el = 1 : size(t, 2)
    dF = ElemLoadVector(p, t, f, el);
    nn = t(1:3, el);
    F(nn) = F(nn) + dF;
  end
  
  
function B = BounMatrix(p, e, b)
  n = size(p, 2);
  B = zeros(n, n);
  for bel = 1 : size(e, 2)
     dB = ElemBounMatrix(p, e, b, bel);
     nn = e(1:2, bel);
     B(nn, nn) = B(nn, nn) + dB;
  end
  
  
function G = BounVector(p, e, g)
  n = size(p, 2);
  G = zeros(n, 1);
  for bel = 1 : size(e, 2)
     dG = ElemBounVector(p, e, g, bel);
     nn = e(1:2, bel);
     G(nn) = G(nn) + dG;
  end

%------------------------------------------------------------------------------------
%
% INCOMPLETE: ElemBounMatrix, ElemLoadVector and ElemBounVector
%
%------------------------------------------------------------------------------------

% Elementmatrix/vector function
%phi*phi = {0,1/4}. gradphi*gradphi = {C1*C2,0} (I think)
% Don't the integral for every base, but do it for every element instead and superposition the contributions.
function dD = ElemDiffMatrix(p, t, d, el)

    % returns the element diffusion matrix dD for element number
    % el defined by t(1:4, el) and p(t(1:3, el).

      NodeCoords = p(:, t(1:3, el));
      d = sum(d(t(1:3, el)))/3;
      Dphi = ElementPhiGradients(NodeCoords);
      dx = area(NodeCoords);               % element area
      dD = d * Dphi' * Dphi * dx;

      
function dF = ElemLoadVector(p, t, f,el)

    %returns the element load vector dF for element number el
    % Approximated with the value at triangle center.
    NodeCoords = p(:, t(1:3, el));
    TriCenter = [NodeCoords(:,1)+NodeCoords(:,2), NodeCoords(:,1)+NodeCoords(:,3), NodeCoords(:,2)+NodeCoords(:,3)] / 2;
    elArea = area(NodeCoords);
    dF = elArea * f(TriCenter);
    
function dB = ElemBounMatrix(p, t, b, el)



function dG = ElemBounVector(p, e, g, el)



function Dphi = ElementPhiGradients(Nodes)

% returns the gradients of the three element basis functions phi on a
% triangle with nodes Nodes

  v = Nodes(:,3)-Nodes(:,1);
  w = Nodes(:,2)-Nodes(:,1);
  gr = v - dot(v,w)*w/norm(w)^2;
  Dphi3 = gr/norm(gr)^2;
  gr = w - dot(w,v)*v/norm(v)^2;
  Dphi2 = gr/norm(gr)^2;
  Dphi1 = - (Dphi2 + Dphi3);
  Dphi = [Dphi1 Dphi2 Dphi3]; 


function ar=area(NodeCoords)

% returns the area of a triangular element with given node coordinates

  x1 = NodeCoords(1,:);
  x2 = NodeCoords(2,:);
  ar = (x1(2)-x1(1))*(x2(3)-x2(1))-(x1(3)-x1(1))*(x2(2)-x2(1));
  ar = ar/2;