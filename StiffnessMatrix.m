function S = StiffnessMatrix(p, t)
    n = size(p, 2);
    S = zeros(n, n);
    for el = 1 : size(t, 2)
        dD = ElemStiffnessMatrix(p, t, el);
        nn = t(1:3, el);
        S(nn, nn) = S(nn, nn) + dD;
    end
end

function dS = ElemStiffnessMatrix(p, t, el)
    NodeCoords = p(:, t(1:3, el));
    Dphi = ElementPhiGradients(NodeCoords);
    dx = ElemArea(NodeCoords);               
    dS = Dphi' * Dphi * dx;
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
   