function [M] = MassMatrix(p,t)

% Generating the mass matrix for a given mesh and base functions-
n = size(p, 2);
M = zeros(n, n);
for el = 1 : size(t, 2)
    dM = ElemMassMatrix(p, t, el);
    nn = t(1:3, el);
    M(nn, nn) = M(nn, nn) + dM;
end

end

function dM = ElemMassMatrix(p,t,el)
NodeCoords = p(:, t(1:3, el));
dx = ElemArea(NodeCoords);               % element area
dM = 1/4 * dx;
end

