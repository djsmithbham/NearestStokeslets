function B=CalcTransMomentCouplingMatrix(x,X,X0,eps,domain,blockSize)

% calculates the matrix giving moments due to translational modes about centre X0

MM1=SolveTransMomentCouplingProblem(x,X,X0,[1 0 0],eps,domain,blockSize);
MM2=SolveTransMomentCouplingProblem(x,X,X0,[0 1 0],eps,domain,blockSize);
MM3=SolveTransMomentCouplingProblem(x,X,X0,[0 0 1],eps,domain,blockSize);

B=[MM1 MM2 MM3];
