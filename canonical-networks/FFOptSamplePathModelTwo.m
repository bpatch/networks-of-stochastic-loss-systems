function [values, solutions, iterationTimes] = FFOptSamplePathModelTwo(k, c, Weights, Mu, Lambda, CapacityCost, MaxIterations, BlockingSamples, deltaTilde, epsilonTilde, BatchTime, BurnInTime,MaxTime, m)
    
    
    
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

    Routes = [1, 6, 6, 6;
              1, 3, 6, 6;
              1, 3, 5, 6;
              1, 4, 6, 6;
              1, 4, 5, 6;
              1, 4, 6, 6;
              2, 6, 6, 6;
              2, 4, 6, 6;
              2, 4, 6, 6;
              2, 3, 6, 6;
              2, 3, 6, 6;
              2, 3, 5, 6];

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
            Blocking = CanonicalEstimatorModelTwo(m, c, deltaTilde, epsilonTilde, BurnInTime, BatchTime, MaxTime, Mu, Lambda);
            BlockingEstimate = BlockingEstimate + Blocking;
        end
        BlockingEstimate = BlockingEstimate./BlockingSamples;
        
        values(j,1) = sum(lambda.*(1-BlockingEstimate(:,1)).*Weights(:,1)) ...
                   +sum(lambda.*BlockingEstimate(:,1).*(1-BlockingEstimate(:,2)).*Weights(:,2)) ... 
                   +sum(lambda.*BlockingEstimate(:,1).*BlockingEstimate(:,2).*(1-BlockingEstimate(:,3)).*Weights(:,3)) ...
                   -dot(c, CapacityCost);
        
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
        fun = @(C) -sum(lambda.*(1-exp(-(C(Routes(:,1))'.*tau(:,1)).^k)).*Weights(:,1)) ...
                   -sum(lambda.*exp(-(C(Routes(:,1))'.*tau(:,1)).^k).*(1-exp(-(C(Routes(:,2))'.*tau(:,2)).^k)).*Weights(:,2)) ... 
                   -sum(lambda.*exp(-(C(Routes(:,1))'.*tau(:,1)).^k).*exp(-(C(Routes(:,2))'.*tau(:,2)).^k).*(1-exp(-(C(Routes(:,3))'.*tau(:,3)).^k)).*Weights(:,3)) ...
                   +dot(C, CapacityCost);
        c = fmincon(fun,c,A,b,Aeq,beq,lb,ub,@mycon,options);
        solutions(j+1,:) = round(c);
        iterationTime = toc;
        iterationTimes(j+1,1) = iterationTime;
    end
    BlockingEstimate = zeros(12,3);
    for i = 1:BlockingSamples
        Blocking = CanonicalEstimatorModelTwo(m, c, deltaTilde, epsilonTilde, BurnInTime, BatchTime, MaxTime, Mu, Lambda);
        BlockingEstimate = BlockingEstimate + Blocking;
    end
    BlockingEstimate = BlockingEstimate./BlockingSamples;

    values(MaxIterations+1,1) = sum(lambda.*(1-BlockingEstimate(:,1)).*Weights(:,1)) ...
                   +sum(lambda.*BlockingEstimate(:,1).*(1-BlockingEstimate(:,2)).*Weights(:,2)) ... 
                   +sum(lambda.*BlockingEstimate(:,1).*BlockingEstimate(:,2).*(1-BlockingEstimate(:,3)).*Weights(:,3)) ...
                   -dot(c, CapacityCost);
end
            