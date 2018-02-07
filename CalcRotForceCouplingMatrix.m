function C=CalcRotForceCouplingMatrix(x,X,X0,ep,domain,blockSize)

F1=SolveRotForceCouplingProblem(x,X,X0,[1 0 0],ep,domain,blockSize);
F2=SolveRotForceCouplingProblem(x,X,X0,[0 1 0],ep,domain,blockSize);
F3=SolveRotForceCouplingProblem(x,X,X0,[0 0 1],ep,domain,blockSize);

C=[F1 F2 F3];
