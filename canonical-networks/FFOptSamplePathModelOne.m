function [values, solutions, iterationTimes] = FFOptSamplePathModelOne(k, c, Weights, Mu, Lambda, CapacityCost, MaxIterations, BlockingSamples, deltaTilde, epsilonTilde, BatchTime, BurnInTime,MaxTime, m)
    
    
    
%     clear all; clc; close all;
%     NumberOfStations = 10;
%     k = 2;
%     Weights = ones(1,NumberOfStations);
%     Weights = [Weights; Weights];
%     Lambda = [10,20,30];
%     Mu = 0.2+0.8.*rand(2,NumberOfStations);
%     alpha = 1;
%     BatchTime = 10;
%     DiscardedBatches = 5;
%     deltaTilde = 0.1;
%     epsilonTilde = 0.0000000001;    
%     BlockingSamples = 5;
%     MaxIterations = 20;
%     CapacityCost = 0.1+0.2.*rand(1,NumberOfStations);
    NumberOfStations = 6;
    options = optimoptions(@fmincon, 'Algorithm', 'interior-point', 'Display', 'off');
    values = NaN(MaxIterations+1, 1);
    iterationTimes = NaN(MaxIterations+1, 1);
    solutions = NaN(MaxIterations+1, NumberOfStations);
    solutions(1,:) = c;
    iterationTimes(1,1) = 0;
    
    lambda = [Lambda(1,1)*ones(6,1); Lambda(1,2)*ones(6,1)];
    for j = 1:MaxIterations
        tic;
        BlockingEstimate = zeros(12,3);
        for i = 1:BlockingSamples
            Blocking = CanonicalEstimatorModelOne(m, c, deltaTilde, epsilonTilde, BurnInTime, BatchTime, MaxTime, Mu, Lambda);
            BlockingEstimate = BlockingEstimate + Blocking;
        end
        BlockingEstimate = BlockingEstimate./BlockingSamples;
        
        values(j,1) = sum(lambda.*prod(1-BlockingEstimate, 2).*Weights) - dot(c, CapacityCost);
        
        tau = (1./[c(1), 0, 0;
              c(1), c(3), 0;
              c(1), c(3), c(5);
              c(1), c(4), 0;
              c(1), c(4), c(5);
              c(1), c(4), c(6);
              c(2), 0, 0;
              c(2), c(4), 0;
              c(2), c(4), c(6);
              c(2), c(3), 0;
              c(2), c(3), c(6);
              c(2), c(3), c(5)]).*(-log(BlockingEstimate)).^(1/k);
        A = [];
        b = []; 
        Aeq = [];
        beq = [];
        lb = zeros(1,6);    
        ub = [];
        fun = @(C) -sum(lambda.*prod(1-exp(-([C(1), inf, inf;
              C(1), C(3), inf;
              C(1), C(3), C(5);
              C(1), C(4), inf;
              C(1), C(4), C(5);
              C(1), C(4), C(6);
              C(2), inf, inf;
              C(2), C(4), inf;
              C(2), C(4), C(6);
              C(2), C(3), inf;
              C(2), C(3), C(6);
              C(2), C(3), C(5)].*tau).^k), 2).*Weights) + dot(C, CapacityCost);
        c = fmincon(fun,c,A,b,Aeq,beq,lb,ub,@mycon,options);
        solutions(j+1,:) = round(c);
        iterationTime = toc;
        iterationTimes(j+1,1) = iterationTime;
    end
    BlockingEstimate = zeros(12,3);
    for i = 1:BlockingSamples
        Blocking = CanonicalEstimatorModelOne(m, c, deltaTilde, epsilonTilde, BurnInTime, BatchTime, MaxTime, Mu, Lambda);
        BlockingEstimate = BlockingEstimate + Blocking;
    end
    BlockingEstimate = BlockingEstimate./BlockingSamples;

    values(MaxIterations+1,1) = sum(lambda.*prod(1-BlockingEstimate, 2).*Weights) - dot(c, CapacityCost);
        
end
            