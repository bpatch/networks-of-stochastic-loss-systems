function [values, solutions, iterationTimes] = StochAppSamplePathModelOne(Weights, CapacityCost, MaxIterations, BlockingSamples, c, lambda, Mu, BatchTime, DiscardedBatches, deltaTilde, epsilonTilde, Alpha2, Alpha3, b2, b3) 
%     
%     c = randi(100, 2, 1)
%     b2 = 150;
    
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
    solutions = NaN(MaxIterations+1, 2);
    solutions(1,:) = c;
    iterationTimes(1,1) = 0;
    
    rho = 0.8;
    cParam = 0.5;
    StepSize = @(t) b2*t.^(-Alpha2/3); 
    EstimateSize = @(t) b3*t.^(-Alpha3/6);
    


    if c(1,1) <= 100 && solutions(2,1) <= 100
        NewPerformance = ExactEvaluationModelOne(round(c(1,1)),round(c(2,1)),Weights, Mu, lambda, CapacityCost);
    else
        BlockingEstimate = zeros(2,1);
        for i = 1:BlockingSamples
            Blocking1 = TandemEstimatorModelOne(c , lambda, Mu, BatchTime, DiscardedBatches, deltaTilde, epsilonTilde);
            BlockingEstimate = BlockingEstimate + Blocking1;
        end
        BlockingEstimate = BlockingEstimate./BlockingSamples;
        NewPerformance  = -c(1,1)*CapacityCost(1,1)-c(2,1)*CapacityCost(2,1)+lambda*Weights(1,1)*(1-BlockingEstimate(1,1))*(1-BlockingEstimate(2,1));
    end    
    
    values(1,1) = NewPerformance;

    for j = 1:MaxIterations
        tic;
        Direction = JacobianEstimatorModelOne(c, EstimateSize(j), lambda, Mu, CapacityCost, Weights, BlockingSamples, BatchTime, DiscardedBatches, CIfunction(j), epsilonTilde);        
        counter = 1;
        StoppingCondition2 = 0;
        while StoppingCondition2 == 0 
            CapacityTest = max(0,c+StepSize(j)*rho^counter.*Direction');
            for i = 1:BlockingSamples
                Blocking1 =  TandemEstimatorModelOne(CapacityTest , lambda, Mu, BatchTime, DiscardedBatches, deltaTilde, epsilonTilde);
                BlockingEstimate = BlockingEstimate + Blocking1;
            end
            BlockingEstimate = BlockingEstimate./BlockingSamples;
            PerformanceTest = -CapacityTest(1,1)*CapacityCost(1,1)-CapacityTest(2,1)*CapacityCost(2,1)+lambda*Weights(1,1)*(1-BlockingEstimate(1,1))*(1-BlockingEstimate(2,1));
            if PerformanceTest >= (NewPerformance+cParam*StepSize(j)*rho^counter*(Direction*Direction'))
                StoppingCondition2 = 1;
                c = CapacityTest;
                NewPerformance = PerformanceTest;
            end
            counter = counter + 1;
            if counter == 20 
                StoppingCondition2 = 1;
%                 c = CapacityTest;
%                 NewPerformance = PerformanceTest;
            end
        end
        solutions(j+1,:) = max(0, round(c));
        iterationTime = toc;
        iterationTimes(j+1,1) = iterationTime;
        
%         c %% comment
        
        if solutions(j+1,1) <= 100 && solutions(j+1,2) <= 20
            values(j+1,1) = ExactEvaluationModelOne(solutions(j+1,1),solutions(j+1,2),Weights, Mu, lambda, CapacityCost);
        else
            BlockingEstimate = zeros(2,1);
            for i = 1:(10*BlockingSamples)
                BlockingEstimate = BlockingEstimate + TandemEstimatorModelOne(solutions(j+1,:)' , lambda, Mu, BatchTime, DiscardedBatches, deltaTilde, epsilonTilde);
            end
            BlockingEstimate = BlockingEstimate/(10*BlockingSamples);
            values(j+1,1) = -c(1,1)*CapacityCost(1,1)-c(2,1)*CapacityCost(2,1)+lambda*Weights(1,1)*(1-BlockingEstimate(1,1))*(1-BlockingEstimate(2,1));
        end
    end
        
%     hold on
%          plot(1:(MaxIterations+1), values(:,1))
end
            