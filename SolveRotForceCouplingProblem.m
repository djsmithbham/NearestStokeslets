function FF=SolveRotForceCouplingProblem(x,X,X0,Om,ep,domain,blockSize)

  % lhs assembly
  [AS, NN]=AssembleStokesletMatrix(x,X,x,ep,domain,blockSize);

  % rhs assembly  
  [x1, x2, x3]=ExtractComponents(x);ze=0*x1;
  x1=x1-X0(1);
  x2=x2-X0(2);
  x3=x3-X0(3);
  rhs=[Om(2)*x3-Om(3)*x2;Om(3)*x1-Om(1)*x3;Om(1)*x2-Om(2)*x1];

  % solve system
  f=AS\rhs;  
  
  % calculate moment
  FF=CalcForceFromForceDistribution(f,NN);
  
