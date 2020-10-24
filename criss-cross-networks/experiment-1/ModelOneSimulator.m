clear all; clc; close all;
% %% The commented code generates data that is stored in ModelOne.mat
rng(2020)
NumberSamples = 10; 
NumberOfStations = 10;
k = 2;
Weights = 10*rand(2,1);
Lambda = [10,20,30];
Mu = 1.2.*ones(2,NumberOfStations);
Alpha = 1;
CapacityCost = 0.1+0.2.*rand(1,NumberOfStations);
BatchTime = 5; %10;
DiscardedBatches = 2;
Alpha2 = 1;
Alpha3 = 1;
b2 = 150;
b3 = 5;
deltaTilde = 0.1;
epsilonTilde = 1/10^6;
BlockingSamples = 5;
MaxIterations = 10; %10;
SuggIterations = 300;% 100;
MaxTime = 100;

lb = zeros(1,10);
ub = 150*ones(1,10);
intcon =1:NumberOfStations;
objconstr = @(x) -ObjectiveFunctionCrissCrossModelOne(Weights, CapacityCost, BlockingSamples, x, Lambda, Mu, Alpha, NumberOfStations, BatchTime, DiscardedBatches, epsilonTilde, deltaTilde,MaxTime);
objconstr_BayesOpt = @(x) -ObjectiveFunctionCrissCrossModelOne(Weights, CapacityCost, BlockingSamples, [x.capacity1, x.capacity2, x.capacity3, x.capacity4, x.capacity5, x.capacity6, x.capacity7, x.capacity8, x.capacity9, x.capacity10], Lambda, Mu, Alpha, NumberOfStations, BatchTime, DiscardedBatches, epsilonTilde, deltaTilde,MaxTime);

FFoptObjData = NaN(NumberSamples, MaxIterations+1);
FFoptCPUData = NaN(NumberSamples, MaxIterations+1);
FFoptSolData = NaN(MaxIterations+1,NumberOfStations,NumberSamples);

SAoptObjData = NaN(NumberSamples, MaxIterations+1);
SAoptCPUData = NaN(NumberSamples, MaxIterations+1);
SAoptSolData = NaN(MaxIterations+1,NumberOfStations,NumberSamples);

SurrOptObjData = NaN(NumberSamples, SuggIterations);
SurrOptCPUData = NaN(NumberSamples, SuggIterations);
SurrOptSolData = NaN(NumberSamples, NumberOfStations);

BayesOptObjData = NaN(NumberSamples, SuggIterations);
BayesOptCPUData = NaN(NumberSamples, SuggIterations);
BayesOptSolData = NaN(NumberSamples, NumberOfStations);

for n = 1:NumberSamples
    c = randi(100, 1, 10);
    [values1, solutions1, iterationTimes1] = FFOptSamplePathModelOne(k, c, Weights, Mu, Lambda, Alpha, CapacityCost, MaxIterations, BlockingSamples, deltaTilde, epsilonTilde, BatchTime, DiscardedBatches, MaxTime);
    FFoptObjData(n, :) = values1';
    FFoptCPUData(n, :) = cumsum(iterationTimes1)';
    FFoptSolData(:,:,n) = solutions1;

    [values2, solutions2, iterationTimes2] = StochAppSamplePathModelOne(Weights, CapacityCost, MaxIterations, BlockingSamples, c, Lambda, Mu, Alpha, NumberOfStations, BatchTime, DiscardedBatches, deltaTilde, epsilonTilde, Alpha2, Alpha3, b2, b3, MaxTime);
    SAoptObjData(n, :) = values2';
    SAoptCPUData(n, :) = cumsum(iterationTimes2)';
    SAoptSolData(:,:,n) = solutions2;
%     
    diary(strcat('SurrOutput_Model1_Sample',num2str(n),'.txt'))

    opts = optimoptions('surrogateopt','MaxFunctionEvaluations',SuggIterations,'Display','iter','InitialPoints',c, 'PlotFcn',[]);
    [x, fval, ~, output2] = surrogateopt(objconstr, lb, ub, intcon, opts) 

    diary off
    
    SurrOptSolData(:,n) = x';

    temp = textread(strcat('SurrOutput_Model1_Sample',num2str(n),'.txt'),'%s','delimiter','\n'); %#ok<*DTXTRD>
    Index = find(contains(temp,'Surrogate'), 1);
    temp = temp(8:(Index-1), 1);

    Index = find(contains(temp,'F-'));

    for i = 1:length(Index)
       tempIndex = find(contains(temp,'F-'), 1);
       temp = [temp(1:(tempIndex-2),1); temp((tempIndex+2):end,1)];
    end

    for i = 1:length(temp)
       foo = strsplit(cell2mat(temp(i)));
       SurrOptCPUData(n,i) = str2double(foo(2));
       SurrOptObjData(n,i) = -str2double(foo(3));
    end

    delete(strcat('SurrOutput_Model1_Sample',num2str(n),'.txt'))
    
    cap1 = optimizableVariable('capacity1', [0,150],'Type','integer');
    cap2 = optimizableVariable('capacity2', [0,150],'Type','integer');
    cap3 = optimizableVariable('capacity3', [0,150],'Type','integer');
    cap4 = optimizableVariable('capacity4', [0,150],'Type','integer');
    cap5 = optimizableVariable('capacity5', [0,150],'Type','integer');
    cap6 = optimizableVariable('capacity6', [0,150],'Type','integer');
    cap7 = optimizableVariable('capacity7', [0,150],'Type','integer');
    cap8 = optimizableVariable('capacity8', [0,150],'Type','integer');
    cap9 = optimizableVariable('capacity9', [0,150],'Type','integer');
    cap10 = optimizableVariable('capacity10', [0,150],'Type','integer');

    capacity1 = c(1,1);
    capacity2 = c(1,2);
    capacity3 = c(1,3);
    capacity4 = c(1,4);
    capacity5 = c(1,5);
    capacity6 = c(1,6);
    capacity7 = c(1,7);
    capacity8 = c(1,8);
    capacity9 = c(1,9);
    capacity10 = c(1,10);

    
    C0_BayesOpt = table(capacity1, capacity2,capacity3, capacity4,capacity5, capacity6,capacity7, capacity8, capacity9, capacity10);
    
    results = bayesopt(objconstr_BayesOpt,[cap1, cap2,cap3, cap4,cap5, cap6,cap7, cap8,cap9, cap10],'Verbose',0,...
    'AcquisitionFunctionName','expected-improvement-plus','InitialX',C0_BayesOpt,'MaxObjectiveEvaluations',SuggIterations);

    BayesOptCPUData(n,:) = cumsum(results.IterationTimeTrace');
    BayesOptObjData(n,:) = -results.ObjectiveMinimumTrace';
    BayesOptSolData(n,:) = [results.XAtMinObjective.capacity1, results.XAtMinObjective.capacity2, ...
                            results.XAtMinObjective.capacity3, results.XAtMinObjective.capacity4, ...
                            results.XAtMinObjective.capacity5, results.XAtMinObjective.capacity6, ...
                            results.XAtMinObjective.capacity7, results.XAtMinObjective.capacity8, ...
                            results.XAtMinObjective.capacity9, results.XAtMinObjective.capacity10, ...
                            ];
    close all;    
    
    disp('Progress:')
    n/NumberSamples
end


save('ModelOneSimulationData.mat')

matrix2latex(CapacityCost, 'out1.tex')