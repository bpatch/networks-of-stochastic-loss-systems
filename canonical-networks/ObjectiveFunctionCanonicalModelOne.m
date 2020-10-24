function value = ObjectiveFunctionCanonicalModelOne(Weights, BlockingSamples, c, Mu,  Lambda, CapacityCost, BatchTime, BurnInTime, epsilonTilde, deltaTilde,MaxTime,m)
    lambda = [Lambda(1,1)*ones(6,1); Lambda(1,2)*ones(6,1)];       
    BlockingEstimate = zeros(12,3);
    for i = 1:BlockingSamples
        Blocking = CanonicalEstimatorModelOne(m, c, deltaTilde, epsilonTilde, BurnInTime, BatchTime, MaxTime, Mu, Lambda);
        BlockingEstimate = BlockingEstimate + Blocking;
    end
    BlockingEstimate = BlockingEstimate./BlockingSamples;
    value = sum(lambda.*prod(1-BlockingEstimate, 2).*Weights) - dot(c, CapacityCost);
end