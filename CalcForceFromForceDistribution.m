function FF=CalculateForceFromForceDistribution(f,NN)

    % f is force distribution from nearest neighbour regularised stokeslet code
    % NN is nearest neighbour matrix
    fq=NN*f;
    [f1 f2 f3]=ExtractComponents(fq);
    FF=[sum(f1); sum(f2); sum(f3)];
    
