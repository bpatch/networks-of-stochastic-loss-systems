function jacobian = JacobianEstimatorModelTwo(c, EstSize, lambda, Mu, CapacityCost, Weights,BlockingSamples, BatchTime, DiscardedBatches, CIgapSize, epsilonTilde)
   objective = @(C,Blocking1) -C(1,1)*CapacityCost(1,1)-C(2,1)*CapacityCost(2,1)+lambda*Weights(1,1)*(1-Blocking1(1,1)) + lambda*Weights(1,1)*(1-Blocking1(2,1))*Blocking1(1,1);
    jacobian = zeros(1,2);
    for station = 1:2
        %seed = randi(100000, 1,1);
        CapacityTempUp = c;
        CapacityTempUp(station,1) = CapacityTempUp(station,1)+EstSize;
        CapacityTempDown = c;
        CapacityTempDown(station,1) = max(0, CapacityTempDown(station,1)-EstSize);
        BlockingEstimateDown = zeros(2,1);
        BlockingEstimateUp = zeros(2,1);
        for i = 1:BlockingSamples
            %rng(seed)
            Blocking1down = TandemEstimatorModelTwo(CapacityTempDown , lambda, Mu, BatchTime, DiscardedBatches, CIgapSize, epsilonTilde);
            %rng(seed)
            Blocking1up = TandemEstimatorModelTwo(CapacityTempUp, lambda, Mu, BatchTime, DiscardedBatches, CIgapSize, epsilonTilde);
            BlockingEstimateDown = BlockingEstimateDown + Blocking1down;
            BlockingEstimateUp = BlockingEstimateUp + Blocking1up;
        end
        BlockingEstimateDown = BlockingEstimateDown/BlockingSamples;
        BlockingEstimateUp = BlockingEstimateUp/BlockingSamples;
        PerformanceUp = objective(CapacityTempUp, BlockingEstimateUp);
        PerformanceDown = objective(CapacityTempDown, BlockingEstimateDown);
        jacobian(1,station) = (PerformanceUp-PerformanceDown)/(2*EstSize);
    end
end