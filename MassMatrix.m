function [M] = MassMatrix(p,t)

% Generating the mass matrix for a given mesh and base functions-
m = size(p, 2);
M = zeros(m, m);
for el = 1 : size(t, 2)
    dM = ElemMassMatrix(p, t, el);
    nn = t(1:3, el);
    M(nn, nn) = M(nn, nn) + dM;
end

end

function dM = ElemMassMatrix(p, t, el)
    % returns the element diffusion matrix dD for element number
    % el defined by t(1:4, el) and p(t(1:3, el).
    NodeCoords = p(:, t(1:3, el));
    dx = ElemArea(NodeCoords);
    dM = dx/12 * [2,1,1;1,2,1;1,1,2];
end
