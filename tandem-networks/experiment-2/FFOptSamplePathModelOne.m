function [values, solutions, iterationTimes] = FFOptSamplePathModelOne(k, c, Weights, Mu, lambda, CapacityCost, MaxIterations, BlockingSamples, deltaTilde, epsilonTilde, BatchTime, DiscardedBatches)
    options = optimoptions(@fmincon, 'Algorithm', 'interior-point', 'Display', 'off');
    values = NaN(MaxIterations+1, 1);
    iterationTimes = NaN(MaxIterations, 1);
    solutions = NaN(MaxIterations+1, 2);
    solutions(1,:) = [max(1,round(c(1,1))),  max(1,round(c(2,1)))];
    for j =  1: MaxIterations
        tic;
        BlockingEstimate = zeros(2,1);
        for i = 1:BlockingSamples
            BlockingEstimate = BlockingEstimate + TandemEstimatorModelOne(c, lambda, Mu, BatchTime, DiscardedBatches, deltaTilde, epsilonTilde);
        end
        BlockingEstimate = BlockingEstimate/BlockingSamples;
        tau = (1./c).*(-log(BlockingEstimate)).^(1/k);
        A = [];
        b = []; 
        Aeq = [];
        beq = [];
        lb = zeros(2,1);    
        ub = [];
        fun = @(C) -lambda*Weights(1,1)*(1-exp(-C(1,1)*tau(1,1)))*(1-exp(-C(2,1)*tau(2,1))) + dot(C,CapacityCost);
        c = fmincon(fun,c,A,b,Aeq,beq,lb,ub,@mycon,options);
        solutions(j+1,:) = [max(0,round(c(1,1))),  max(0,round(c(2,1)))];
        iterationTime = toc;
        iterationTimes(j,1) = iterationTime;
    end
    for j = 0:MaxIterations
        if solutions(j+1,1) <= 100 && solutions(j+1,2) <= 100
            values(j+1,1) = ExactEvaluationModelOne(solutions(j+1,1),solutions(j+1,2),Weights, Mu, lambda, CapacityCost);
        else
            BlockingEstimate = zeros(2,1);
            for i = 1:(10*BlockingSamples)
                BlockingEstimate = BlockingEstimate + TandemEstimatorModelOne(solutions(j+1,:)' , lambda, Mu, BatchTime, DiscardedBatches, deltaTilde, epsilonTilde);
            end
            BlockingEstimate = BlockingEstimate/(10*BlockingSamples);
            values(j+1,1) = -solutions(j+1,1)*CapacityCost(1,1)-solutions(j+1,2)*CapacityCost(2,1)+lambda*Weights(1,1)*(1-BlockingEstimate(1,1))*(1-BlockingEstimate(2,1));
        end
    end
end