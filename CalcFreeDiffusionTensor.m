function D=CalcFreeDiffusionTensor(x,X,X0,kB,mu,T,lenscale,ep,domain,blockSize)

% calculates dimensional free diffusion tensor about X0
%
% procedure from Wegener, Biopolymers, 20:303–326, 1981: 
%        take D(4:6,4:6) to get point-independent rotational diffusion matrix

        A=CalcTransResistanceMatrix(x,X,ep,domain,blockSize);
        B=CalcTransMomentCouplingMatrix(x,X,X0,ep,domain,blockSize);
        C=CalcRotResistanceMatrix(x,X,X0,ep,domain,blockSize);
        D=kB*T*inv([A B';B C])/mu/lenscale^3;
