function A=CalcTransResistanceMatrix(x,X,eps,domain,blockSize)

% Calculates the translational resistance matrix for an object from its body
% frame points
% input:  theta, phi spherical polar angles
%         x          body frame collocation/force points
%         X          body frame stokeslet points
%         eps        numerical regularisation length
%         domain     'i' for infinite fluid, 'h' for half space above no-slip boundary (x3>0)
%         blockSize  controls matrix assembly blocking size for numerical efficiency, 0.2 is a safe choice
% output: A          3x3 matrix of translational resistance coefficients 
%                    (scaled wrt 1 / mu)

F1=SolveTranslationalResistanceProblem(x,X,[1 0 0],eps,domain,blockSize);
F2=SolveTranslationalResistanceProblem(x,X,[0 1 0],eps,domain,blockSize);
F3=SolveTranslationalResistanceProblem(x,X,[0 0 1],eps,domain,blockSize);

A=[F1 F2 F3];