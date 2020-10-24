function value = ObjectiveFunctionCrissCrossModelTwo(Weights, CapacityCost, BlockingSamples, c, Lambda, Mu, alpha, NumberOfStations, BatchTime, DiscardedBatches, epsilonTilde, deltaTilde,MaxTime)
    lambda = 0.5*mean(Lambda);        
    BlockingEstimate1 = zeros(2,NumberOfStations);
    for i = 1:BlockingSamples
        [Blocking1, ~] = CrissCrossEstimatorModelTwo(c, Lambda, Mu, alpha, NumberOfStations, BatchTime, DiscardedBatches, deltaTilde, epsilonTilde,MaxTime);
        BlockingEstimate1 = BlockingEstimate1 + Blocking1;
    end
    BlockingEstimate = BlockingEstimate1./BlockingSamples;
    value = lambda*(1-BlockingEstimate(1,1))*Weights(1,1) ...
                   + lambda*sum(cumprod(BlockingEstimate(1,1:(end-1))).*Weights(1,2:end).*(1-BlockingEstimate(1,2:end))) ...
                   + lambda*(1-BlockingEstimate(2,end))*Weights(2,end) ...
                   + lambda*sum(cumprod(fliplr(BlockingEstimate(2,2:end))).*fliplr(Weights(2,1:(end-1))).*fliplr(1-BlockingEstimate(2,1:(end-1)))) ...
                   - dot(c,CapacityCost);
end