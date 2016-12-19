function [F,S,M,BoundNodesCompl] = Dirichlet(F,S,M,e)
% Truncates all the rows (and columns for S,M) of the vector (and matrices)
% corresponding to vertices on the boundary dOmega. Also returns indices
% for all elements not on the boundary to be used to single out which xi:s
% to be updated in the iteration.

m = size(M,1);
BoundNodes = [e(1,:) e(2,:)];
BoundNodes = unique(BoundNodes);

BoundNodesCompl = comple(BoundNodes,m);

S(:,BoundNodes) = [];
S(BoundNodes,:) = [];
M(:,BoundNodes) = [];
M(BoundNodes,:) = [];
F(BoundNodes,:)= [];


end

