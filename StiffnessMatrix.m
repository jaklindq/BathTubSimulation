function [S] = StiffnessMatrix(p,t)

% Generating the stiffness matrix (corresponding to the laplace operator)
% for a given mesh and base functions-
% returns the element diffusion matrix dD for element number
% el defined by t(1:4, el) and p(t(1:3, el).
n = size(p, 2);
S = zeros(n, n);
    for el = 1 : size(t, 2)
        dS = ElemStiffnessMatrix(p, t, el);
        nn = t(1:3, el);
        S(nn, nn) = S(nn, nn) + dS;
    end

end

function dS = ElemStiffnessMatrix(p,t,el)
    NodeCoords = p(:, t(1:3, el));
    % d = sum(d(t(1:3, el)))/3;
    dPhi = ElementPhiGradients(NodeCoords); % REPLACE
    dx = ElemArea(NodeCoords);               % element area
    dS = dPhi' * dPhi * dx;
end

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
end