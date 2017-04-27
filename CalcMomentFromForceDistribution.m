function MM=CalculateMomentFromForceDistribution(f,X,X0,NN)

  [X1 X2 X3]=ExtractComponents(X'*NN);Ze=0*X1;
  X1=X1-X0(1);
  X2=X2-X0(2);
  X3=X3-X0(3);
  AM=[Ze -X3 X2; X3 Ze -X1; -X2 X1 Ze];
  MM=AM*f;
