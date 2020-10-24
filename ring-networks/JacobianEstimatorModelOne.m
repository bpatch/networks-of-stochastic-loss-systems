function jacobien = JacobianEstimatorModelOne(c, StepSize, Lambda, Mu, Alpha, CapacityCost, Weights, NumberOfStations,BlockingSamples, BatchTime, DiscardedBatches, CIgap, epsilonTilde,MaxTime)
   lambda = 0.5*mean(Lambda); % Divide by 2 since arrivals are split between the two routes
   objective = @(c,Blocking1) lambda*prod(1-Blocking1(1,:)).*Weights(1,1) ...
                      + lambda*prod(1-Blocking1(2,:)).*Weights(2,1) ...
                      - dot(c,CapacityCost);
    jacobien = zeros(1,NumberOfStations);
    for station = 1:NumberOfStations
        %seed = randi(100000, 1,1);
        CapacityTempUp = c;
        CapacityTempUp(1,station) = CapacityTempUp(1,station)+StepSize;
        CapacityTempDown = c;
        CapacityTempDown(1,station) = max(0, CapacityTempDown(1,station)-StepSize);
        BlockingEstimateDown = zeros(2,NumberOfStations);
        BlockingEstimateUp = zeros(2,NumberOfStations);
        for i = 1:BlockingSamples
            %rng(seed)
            [Blocking1down, ~] = CrissCrossEstimatorModelOne(CapacityTempDown, Lambda, Mu, Alpha, NumberOfStations, BatchTime, DiscardedBatches, CIgap, epsilonTilde,MaxTime);
            %rng(seed)
            [Blocking1up, ~] = CrissCrossEstimatorModelOne(CapacityTempUp, Lambda, Mu, Alpha, NumberOfStations, BatchTime, DiscardedBatches, CIgap, epsilonTilde,MaxTime);
            BlockingEstimateDown = BlockingEstimateDown + Blocking1down;
            BlockingEstimateUp = BlockingEstimateUp + Blocking1up;
        end
        BlockingEstimateDown = BlockingEstimateDown/BlockingSamples;
        BlockingEstimateUp = BlockingEstimateUp/BlockingSamples;
        PerformanceUp = objective(CapacityTempUp, BlockingEstimateUp);
        PerformanceDown = objective(CapacityTempDown, BlockingEstimateDown);
        jacobien(1,station) = (PerformanceUp-PerformanceDown)/(2*StepSize);
    end
end