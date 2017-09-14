function MM=CalcMomentFromForceDistribution(f,X,X0,NN)

  % modification 2017_SEP_04 to translate points before
  % multiplying by nearest-neighbour matrix

  [X1,X2,X3]=ExtractComponents(X);
  X1=X1-X0(1);
  X2=X2-X0(2);
  X3=X3-X0(3);
  
  x=[X1' X2' X3']*NN;
  [x1,x2,x3]=ExtractComponents(x);Ze=0*x1;
  
  AM=[Ze -x3 x2; x3 Ze -x1; -x2 x1 Ze];
  MM=AM*f;
