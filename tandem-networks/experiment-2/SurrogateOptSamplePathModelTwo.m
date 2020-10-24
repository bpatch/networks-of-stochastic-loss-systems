function [values, solutions, iterationTimes] = SurrogateOptSamplePathModelTwo(Weights, CapacityCost, MaxIterations, BlockingSamples, c, lambda, Mu, BatchTime, DiscardedBatches, deltaTilde, epsilonTilde)
%     clear all; clc; close all;
% rng(2020)
% NumberOfPaths = 2; 
% NumberOfStations = 10;
% k = 2;
% Weights = rand(2,NumberOfStations);
% Lambda = [10,20,30];
% Mu = 1.2.*ones(2,NumberOfStations);
% Alpha = 1;
% BatchTime = 10;
% DiscardedBatches = 1;
% deltaTilde = 0.03;
% Alpha2 = 1;
% Alpha3 = 1;
% b2 = 150;
% b3 = 5;
% epsilonTilde = 0.0000000001;    
% BlockingSamples = 10;
% MaxIterations = 10;
% CapacityCost = 0.1+0.2.*rand(1,NumberOfStations);
% 
% C = randi(100,1,NumberOfStations);

    CIfunction = @(j) deltaTilde;%/sqrt(j);


    values = NaN(MaxIterations+1, 1);
    values(1,1) = ObjectiveFunctionTandemModelTwo(Weights, CapacityCost, BlockingSamples, c, lambda, Mu, BatchTime, DiscardedBatches, epsilonTilde, CIfunction(1));
    iterationTimes = NaN(MaxIterations+1, 1);
    solutions = NaN(MaxIterations+1, 2);
    solutions(1,:) = c;
    iterationTimes(1,1) = 0;

    lb = zeros(1,2);
    ub = 100*ones(1,2);
    intcon =1:2;
    seed = randi(10^5);
    for j = 1:MaxIterations
        rng(seed);
        if j == 1
            objconstr = @(x) -ObjectiveFunctionTandemModelTwo(Weights, CapacityCost, BlockingSamples, x', lambda, Mu, BatchTime, DiscardedBatches, epsilonTilde, CIfunction(1));
            opts = optimoptions('surrogateopt','MaxFunctionEvaluations',1,'Display','off','InitialPoints',c','PlotFcn',[]);
            [x,fval,~,output] = surrogateopt(objconstr, lb, ub, intcon, opts);
            solutions(j+1,:) = x;
            values(j+1,1) = -fval;
            previousTime = output.elapsedtime;
            iterationTimes(2,1) = previousTime;
        else
    %         objconstr = @(x) -ObjectiveFunctionTandemModelOne(Weights, CapacityCost, BlockingSamples, x', lambda, Mu, BatchTime, DiscardedBatches, epsilonTilde, CIfunction(j));
    %         opts = optimoptions('surrogateopt','MaxFunctionEvaluations', j-1, 'Display', 'off','InitialPoints',c', 'PlotFcn',[]);
    %         [~, ~, ~, output1] = surrogateopt(objconstr, lb, ub, intcon, opts);

            opts = optimoptions('surrogateopt','MaxFunctionEvaluations',j,'Display','off','InitialPoints',c', 'PlotFcn',[]);
            [x, fval, ~, output2] = surrogateopt(objconstr, lb, ub, intcon, opts);
            solutions(j+1,:) = x;
            values(j+1, 1) = -fval;
            iterationTimes(j+1, 1) = output2.elapsedtime-previousTime;
            previousTime = output2.elapsedtime;
        end
    %     j/MaxIterations
    end

end


