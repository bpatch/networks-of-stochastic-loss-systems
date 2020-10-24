function value = ObjectiveFunctionCanonicalModelTwo(Weights, BlockingSamples, c, Mu,  Lambda, CapacityCost, BatchTime, BurnInTime, epsilonTilde, deltaTilde,MaxTime,m)
    lambda = [Lambda(1,1)*ones(6,1); Lambda(1,2)*ones(6,1)];       
    BlockingEstimate = zeros(12,3);
    for i = 1:BlockingSamples
        Blocking = CanonicalEstimatorModelTwo(m, c, deltaTilde, epsilonTilde, BurnInTime, BatchTime, MaxTime, Mu, Lambda);
        BlockingEstimate = BlockingEstimate + Blocking;
    end
    BlockingEstimate = BlockingEstimate./BlockingSamples;
    value = sum(lambda.*(1-BlockingEstimate(:,1)).*Weights(:,1)) ...
                   +sum(lambda.*BlockingEstimate(:,1).*(1-BlockingEstimate(:,2)).*Weights(:,2)) ... 
                   +sum(lambda.*BlockingEstimate(:,1).*BlockingEstimate(:,2).*(1-BlockingEstimate(:,3)).*Weights(:,3)) ...
                   -dot(c, CapacityCost);
end