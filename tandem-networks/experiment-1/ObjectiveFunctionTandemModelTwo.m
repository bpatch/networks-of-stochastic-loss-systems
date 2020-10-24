function value = ObjectiveFunctionTandemModelTwo(Weights, CapacityCost, BlockingSamples, c, lambda, Mu, BatchTime, DiscardedBatches, epsilonTilde, deltaTilde)
    BlockingEstimate = zeros(2,1);
    for i = 1:BlockingSamples
        Blocking1 = TandemEstimatorModelTwo(c, lambda, Mu, BatchTime, DiscardedBatches, deltaTilde, epsilonTilde);
        BlockingEstimate = BlockingEstimate + Blocking1;
    end
    BlockingEstimate = BlockingEstimate./BlockingSamples;
    value = -c(1,1)*CapacityCost(1,1)-c(2,1)*CapacityCost(2,1)+lambda*Weights(1,1)*(1-BlockingEstimate(1,1))+lambda*Weights(2,1)*(1-BlockingEstimate(2,1))*BlockingEstimate(1,1);
end