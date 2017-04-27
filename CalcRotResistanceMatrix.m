function C=CalcRotResistanceMatrix(x,X,X0,ep,domain,blockSize)

M1=SolveRotationalResistanceProblem(x,X,X0,[1 0 0],ep,domain,blockSize);
M2=SolveRotationalResistanceProblem(x,X,X0,[0 1 0],ep,domain,blockSize);
M3=SolveRotationalResistanceProblem(x,X,X0,[0 0 1],ep,domain,blockSize);

C=[M1 M2 M3];
