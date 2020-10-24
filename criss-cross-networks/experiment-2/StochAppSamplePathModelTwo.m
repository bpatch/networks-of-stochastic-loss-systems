function [values, solutions, iterationTimes] = StochAppSamplePathModelTwo(Weights, CapacityCost, MaxIterations, BlockingSamples, c, Lambda, Mu, Alpha, NumberOfStations, BatchTime, DiscardedBatches, deltaTilde, epsilonTilde, Alpha2, Alpha3, b2, b3,MaxTime) 
    
    
    
%clear all; clc; close all;
% NumberOfStations = 2;
% k = 2;
% Weights = ones(1,NumberOfStations);
% Weights = [Weights; Weights];
% Lambda = [10,20,30];
% Mu = 0.2+0.8.*ones(2,NumberOfStations);
% alpha = 1;
% BatchTime = 10;
% DiscardedBatches = 5;
% deltaTilde = 0.01;
% epsilonTilde = 0.0000000001;    
% BlockingSamples = 20;
% MaxIterations = 3;
% CapacityCost = 0.1+0.2.*rand(1,NumberOfStations);
% C = randi(100, 1, NumberOfStations);
% [values1, solutions1, iterationTimes1] = ModWeibullSamplePath(k, Weights, CapacityCost, MaxIterations, BlockingSamples, C, Lambda, Mu, alpha, NumberOfStations, BatchTime, DiscardedBatches, deltaTilde, epsilonTilde);
% 
% save('SystemParameters.mat')
    CIfunction = @(j) deltaTilde;%/sqrt(j);
    values = NaN(MaxIterations+1, 1);
    iterationTimes = NaN(MaxIterations+1, 1);
    solutions = NaN(MaxIterations+1, NumberOfStations);
    solutions(1,:) = c;
    iterationTimes(1,1) = 0;
    lambda = 0.5*mean(Lambda);
    
    rho = 0.8;
    cParam = 0.5;
    StepSize = @(t) b2*t.^(-Alpha2/3); 
    EstimateSize = @(t) b3*t.^(-Alpha3/6);
    
    BlockingEstimate = zeros(2,NumberOfStations);
    for i = 1:BlockingSamples
        [Blocking1, ~] = CrissCrossEstimatorModelTwo(c, Lambda, Mu, Alpha, NumberOfStations, BatchTime, DiscardedBatches, CIfunction(MaxIterations), epsilonTilde,MaxTime);
        BlockingEstimate = BlockingEstimate + Blocking1;
    end
    BlockingEstimate = BlockingEstimate./BlockingSamples;

    NewPerformance = lambda*(1-BlockingEstimate(1,1))*Weights(1,1) ...
                   + lambda*sum(cumprod(BlockingEstimate(1,1:(end-1))).*Weights(1,2:end).*(1-BlockingEstimate(1,2:end))) ...
                   + lambda*(1-BlockingEstimate(2,end))*Weights(2,end) ...
                   + lambda*sum(cumprod(fliplr(BlockingEstimate(2,2:end))).*fliplr(Weights(2,1:(end-1))).*fliplr(1-BlockingEstimate(2,1:(end-1)))) ...
                   - dot(c,CapacityCost);
    
    values(1,1) = NewPerformance;

    for j = 1:MaxIterations
        tic;
        Direction = JacobianEstimatorModelTwo(c, EstimateSize(j), Lambda, Mu, Alpha, CapacityCost, Weights, NumberOfStations, BlockingSamples, BatchTime, DiscardedBatches, CIfunction(j), epsilonTilde,MaxTime);        
        counter = 1;
        StoppingCondition2 = 0;
        while StoppingCondition2 == 0 
            CapacityTest = max(0,c+StepSize(j)*rho^counter.*Direction);
            for i = 1:BlockingSamples
                [Blocking1, ~] = CrissCrossEstimatorModelTwo(CapacityTest, Lambda, Mu, Alpha, NumberOfStations, BatchTime, DiscardedBatches, CIfunction(j), epsilonTilde,MaxTime);
                BlockingEstimate = BlockingEstimate + Blocking1;
            end
            BlockingEstimate = BlockingEstimate./BlockingSamples;
            PerformanceTest = lambda*(1-BlockingEstimate(1,1))*Weights(1,1) ...
                   + lambda*sum(cumprod(BlockingEstimate(1,1:(end-1))).*Weights(1,2:end).*(1-BlockingEstimate(1,2:end))) ...
                   + lambda*(1-BlockingEstimate(2,end))*Weights(2,end) ...
                   + lambda*sum(cumprod(fliplr(BlockingEstimate(2,2:end))).*fliplr(Weights(2,1:(end-1))).*fliplr(1-BlockingEstimate(2,1:(end-1)))) ...
                   - dot(CapacityTest,CapacityCost);
            if PerformanceTest >= (NewPerformance+cParam*StepSize(j)*rho^counter*(Direction*Direction'))
                StoppingCondition2 = 1;
                c = CapacityTest;
                NewPerformance = PerformanceTest;
            end
            counter = counter + 1;
            if counter == 20
                StoppingCondition2 = 1;
            end
        end
        solutions(j+1,:) = max(0, c);
        iterationTime = toc;
        iterationTimes(j+1,1) = iterationTime;
        BlockingEstimate = zeros(2,NumberOfStations);
        for i = 1:BlockingSamples
            [Blocking1, ~] = CrissCrossEstimatorModelTwo(c, Lambda, Mu, Alpha, NumberOfStations, BatchTime, DiscardedBatches, CIfunction(MaxIterations), epsilonTilde,MaxTime);
            BlockingEstimate = BlockingEstimate + Blocking1;
        end
        BlockingEstimate = BlockingEstimate./BlockingSamples;

        values(j+1,1) = lambda*(1-BlockingEstimate(1,1))*Weights(1,1) ...
                   + lambda*sum(cumprod(BlockingEstimate(1,1:(end-1))).*Weights(1,2:end).*(1-BlockingEstimate(1,2:end))) ...
                   + lambda*(1-BlockingEstimate(2,end))*Weights(2,end) ...
                   + lambda*sum(cumprod(fliplr(BlockingEstimate(2,2:end))).*fliplr(Weights(2,1:(end-1))).*fliplr(1-BlockingEstimate(2,1:(end-1)))) ...
                   - dot(c,CapacityCost);
    end
    BlockingEstimate = zeros(2,NumberOfStations);
    for i = 1:BlockingSamples
        [Blocking1, ~] = CrissCrossEstimatorModelTwo(c, Lambda, Mu, Alpha, NumberOfStations, BatchTime, DiscardedBatches, CIfunction(MaxIterations), epsilonTilde,MaxTime);
        BlockingEstimate = BlockingEstimate + Blocking1;
    end
    BlockingEstimate = BlockingEstimate./BlockingSamples;

    values(MaxIterations+1,1) = lambda*(1-BlockingEstimate(1,1))*Weights(1,1) ...
                   + lambda*sum(cumprod(BlockingEstimate(1,1:(end-1))).*Weights(1,2:end).*(1-BlockingEstimate(1,2:end))) ...
                   + lambda*(1-BlockingEstimate(2,end))*Weights(2,end) ...
                   + lambda*sum(cumprod(fliplr(BlockingEstimate(2,2:end))).*fliplr(Weights(2,1:(end-1))).*fliplr(1-BlockingEstimate(2,1:(end-1)))) ...
                   - dot(c,CapacityCost);
%     hold on
%          plot(1:(MaxIterations+1), values(:,1))
end
            