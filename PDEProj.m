model = createpde(1);
geometryFromEdges(model,@lshapeg);
mesh = generateMesh(model);
[p, e, t] = meshToPet(mesh);
pdeplot(model)


% model = createpde(1);
% [p,e,t] = initmesh(model)
% [p,e,t] = initmesh(model)
