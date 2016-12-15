function [F] = LoadVector(p,t)
% Generating the load vector for a given function f(t,x), mesh and
% base functions. Will have to be updated each time for a time dependent f
n = size(p, 2);
F = zeros(n, 1);


for el = 1 : size(t, 2) % el = triangel
    dF = ElemLoadVector(p, t, el);
    nn = t(1:3, el);
    F(nn) = F(nn) + dF;
    
end

function dF =  ElemLoadVector(p,t,el)
        nodes = p(:,t(1:3,el) );
%         ind0 = (n-1)*round(rand(1,1),2)+1

        center = [sum(nodes(1,:)), sum(nodes(2,:))]'/3;
        dF = ElemArea(nodes)/2 * SourceFunction(center);
end

end        

%%%%%%%%
%PART II
%%%%%%%%
% n = size(p, 2);
% F = zeros(n, 1);
% % mod(time/dt,10)
% if (mod(time/dt,100)) == 1
%     ind0 = randi(n);
%    F(ind0) = -0.000001;
% else
%    F = zeros(n,1);
% end