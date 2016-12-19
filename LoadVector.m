function [F] = LoadVector(p,t,hmax,sourceValue,dropInd,dropTime,endTime)
    % Generating the load vector for a given function f(t,x), mesh and
    % base functions. Will have to be updated each time for a time dependent f
    
    
    
    m = size(p, 2);
    F = zeros(m, endTime);
        for el = 1 : size(t, 2) % el = triangel
            dF = ElemLoadVector(p, t,el,hmax,sourceValue,dropInd);
            F(dropInd, dropTime) = F(dropInd, dropTime) + dF;            
        end
end

function dF =  ElemLoadVector(p,t,el,hmax,sourceValue,dropInd)
    f = sourceValue/hmax^2; % Scale source with geometry.
    nodes = p(:,t(1:3,el) );
    
    % Element containing drop index
    isDropInd = sum(t(1:3,el)==dropInd);    
    dF = isDropInd*f*ElemArea(nodes)/2;
end
