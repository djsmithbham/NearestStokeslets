function [x,X]=GenerateSpherePoints(m,q,a)

  % generates points on surface of sphere using cube projection
  % Input:   m,    number of subdivisions per face for collocation/traction points
  %          q,    number of subdivisions per face for stokeslets (quadrature) points
  %          a,    sphere radius
  % Output:  x,    array of collocation/traction points, with all x1 coordinates first, then all x2 coordinates...
  %          X,    array of stokeslet points

  [s,xc,yc,zc]=UniformSphereGen(1,m);
  [s,xq,yq,zq]=UniformSphereGen(1,q);

  M=6*m*m;
  Q=6*q*q;

  x=a*[reshape(permute(xc,[2,3,1]),M,1);reshape(permute(yc,[2,3,1]),M,1);reshape(permute(zc,[2,3,1]),M,1)]; 
  X=a*[reshape(permute(xq,[2,3,1]),Q,1);reshape(permute(yq,[2,3,1]),Q,1);reshape(permute(zq,[2,3,1]),Q,1)];
