function MM=SolveTransMomentCouplingProblem(x,X,X0,U,ep,domain,blockSize)

% calculates the moment MM about X0 due to translation with velocity U

  % lhs assembly
  [AS NN]=AssembleStokesletMatrix(x,X,x,ep,domain,blockSize);
  
  % rhs assembly
  [x1 x2 x3]=ExtractComponents(x);ze=0*x1;
  rhs=[ze+U(1);ze+U(2);ze+U(3)];

  % solve system
  f=AS\rhs;  
  
  % calculate moment
  MM=CalcMomentFromForceDistribution(f,X,X0,NN);
