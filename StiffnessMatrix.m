function [S] = StiffnessMatrix(p,e,t)

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
dPhi = 1;% ElementPhiGradients(NodeCoords); % REPLACE
dx = ElemArea(NodeCoords);               % element area
dD = dPhi' * dPhi * dx;
end