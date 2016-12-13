function [F] = LoadVector(f,p,t)
% Generating the load vector for a given function f(t,x), mesh and
% base functions. Will have to be updated each time for a time dependent f

n = size(p, 2);
F = zeros(n, 1);
for el = 1 : size(t, 2) % el = triangel
    dF = ElemLoadVector(p, t, f, el);
    nn = t(1:3, el);
    F(nn) = F(nn) + dF;
    
end

function dF =  ElemLoadVector(p,t,f,el)
        nodes = p(:,t(1:3,el) );
        center = [sum(nodes(1,:)), sum(nodes(2,:))]/3;
        dF = ElemArea(nodes) * 0;%f(center);
end
end        
