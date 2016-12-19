function ar=ElemArea(NodeCoords)

% returns the area of a triangular element with given node coordinates

  x1 = NodeCoords(1,:);
  x2 = NodeCoords(2,:);
  ar = (x1(2)-x1(1))*(x2(3)-x2(1))-(x1(3)-x1(1))*(x2(2)-x2(1));
  ar = ar/2;
end