function [values, solutions, iterationTimes] = FFOptSamplePathModelTwo(k, c, Weights, Mu, Lambda, Alpha, CapacityCost, MaxIterations, BlockingSamples, deltaTilde, epsilonTilde, BatchTime, DiscardedBatches, MaxTime)
    
    
    
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
    NumberOfStations = length(CapacityCost);
    CIfunction = @(j) deltaTilde;%/sqrt(j);
    options = optimoptions(@fmincon, 'Algorithm', 'interior-point', 'Display', 'off');
    values = NaN(MaxIterations+1, 1);
    iterationTimes = NaN(MaxIterations+1, 1);
    solutions = NaN(MaxIterations+1, NumberOfStations);
    solutions(1,:) = c;
    iterationTimes(1,1) = 0;
    
    lambda = 0.5*mean(Lambda);
    for j = 1:MaxIterations
        tic;
        BlockingEstimate = zeros(2,NumberOfStations);
        for i = 1:BlockingSamples
            [Blocking1, ~] = CrissCrossEstimatorModelTwo(c, Lambda, Mu, Alpha, NumberOfStations, BatchTime, DiscardedBatches, CIfunction(1), epsilonTilde, MaxTime);
            BlockingEstimate = BlockingEstimate + Blocking1;
        end
        BlockingEstimate = BlockingEstimate./BlockingSamples;
        
        values(j,1) =  lambda*(1-BlockingEstimate(1,1))*Weights(1,1) ...
                   + lambda*sum(cumprod(BlockingEstimate(1,1:(end-1))).*Weights(1,2:end).*(1-BlockingEstimate(1,2:end))) ...
                   + lambda*(1-BlockingEstimate(2,end))*Weights(2,end) ...
                   + lambda*sum(cumprod(fliplr(BlockingEstimate(2,2:end))).*fliplr(Weights(2,1:(end-1))).*fliplr(1-BlockingEstimate(2,1:(end-1)))) ...
                   - dot(c,CapacityCost);
        
        tau = [(1./c).*(-log(BlockingEstimate(1,:))).^(1/k); (1./c).*(-log(BlockingEstimate(2,:))).^(1/k)];
        A = [];
        b = []; 
        Aeq = [];
        beq = [];
        lb = zeros(1,NumberOfStations);    
        ub = [];
        fun = @(c) - lambda*(1-exp(-(c(1,1)*tau(1,1))^k))*Weights(1,1) ...
                   - lambda*sum(cumprod(exp(-(c(1:(end-1)).*tau(1,1:(end-1))).^k)).*Weights(1,2:end).*(1-exp(-(c(2:end).*tau(1,2:end).^k)))) ...
                   - lambda*(1-exp(-(c(1,end)*tau(2,end))^k))*Weights(2,end) ...
                   - lambda*sum(cumprod(fliplr(exp(-(c(2:end).*tau(2,2:end)).^k))).*fliplr(Weights(2,1:(end-1))).*fliplr(1-exp(-(c(1:(end-1)).*tau(2,1:(end-1))).^k))) ...
                   + dot(c,CapacityCost);
        c = fmincon(fun,c,A,b,Aeq,beq,lb,ub,@mycon,options);
        solutions(j+1,:) = round(c);
        iterationTime = toc;
        iterationTimes(j+1,1) = iterationTime;
    end
    BlockingEstimate = zeros(2,NumberOfStations);
    for i = 1:BlockingSamples
        [Blocking1, ~] = CrissCrossEstimatorModelTwo(c, Lambda, Mu, Alpha, NumberOfStations, BatchTime, DiscardedBatches, CIfunction(j), epsilonTilde, MaxTime);
        BlockingEstimate = BlockingEstimate + Blocking1;
    end
    BlockingEstimate = BlockingEstimate./BlockingSamples;

    values(MaxIterations+1,1) = lambda*(1-BlockingEstimate(1,1))*Weights(1,1) ...
                   + lambda*sum(cumprod(BlockingEstimate(1,1:(end-1))).*Weights(1,2:end).*(1-BlockingEstimate(1,2:end))) ...
                   + lambda*(1-BlockingEstimate(2,end))*Weights(2,end) ...
                   + lambda*sum(cumprod(fliplr(BlockingEstimate(2,2:end))).*fliplr(Weights(2,1:(end-1))).*fliplr(1-BlockingEstimate(2,1:(end-1)))) ...
                   - dot(c,CapacityCost);
end
            