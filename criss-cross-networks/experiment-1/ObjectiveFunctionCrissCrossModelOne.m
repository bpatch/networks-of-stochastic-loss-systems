function value = ObjectiveFunctionCrissCrossModelOne(Weights, CapacityCost, BlockingSamples, c, Lambda, Mu, alpha, NumberOfStations, BatchTime, DiscardedBatches, epsilonTilde, deltaTilde,MaxTime)
    lambda = 0.5*mean(Lambda);        
    BlockingEstimate1 = zeros(2,NumberOfStations);
    for i = 1:BlockingSamples
        [Blocking1, ~] = CrissCrossEstimatorModelOne(c, Lambda, Mu, alpha, NumberOfStations, BatchTime, DiscardedBatches, deltaTilde, epsilonTilde,MaxTime);
        BlockingEstimate1 = BlockingEstimate1 + Blocking1;
    end
    BlockingEstimate = BlockingEstimate1./BlockingSamples;
    value = lambda*prod(1-BlockingEstimate(1,:)).*Weights(1,1) ...
                      + lambda*prod(1-BlockingEstimate(2,:)).*Weights(2,1) ...
                      - dot(c,CapacityCost);
end