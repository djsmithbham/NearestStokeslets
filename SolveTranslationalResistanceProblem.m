function [FF]=SolveTranslationalResistanceProblem(x,X,U,ep,domain,blockSize)

  % Solves the problem of finding the resistance on a body for the 3 translational modes
  % Input:  x          collocation/force grid
  %         X          stokeslet grid
  %         U          imposed body velocity vector
  %         ep         regularisation parameter
  %         domain     'i' infinite fluid (regularised stokeslet) or 'h' halfspace x3 > 0 (regularised blakelet)
  %         blockSize  maximum memory storage for stokeslet matrix in GB
  % Output: FF         force due to translation

  % lhs assembly
  [AS NN]=AssembleStokesletMatrix(x,X,x,ep,domain,blockSize);
  
  % rhs assembly
  [x1 x2 x3]=ExtractComponents(x);ze=0*x1;
  rhs=[ze+U(1);ze+U(2);ze+U(3)];

  % solve system
  f=AS\rhs;  
  
  % calculate force
  FF=CalcForceFromForceDistribution(f,NN);
  
