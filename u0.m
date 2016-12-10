function [ xi_0] = u0(p)

xi_0 = cos( 5*pi*sqrt( sum(p.^2) ) ) ./ ( 1 + 10* sqrt( sum(p.^2) ) );

end

