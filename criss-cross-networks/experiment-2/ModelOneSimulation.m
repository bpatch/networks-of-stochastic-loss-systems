clear all; clc; close all;
rng(2020)
NumberSamples = 10; 
NumberScenarios = 20;
NumberOfStations = 10;
k = 2;
BatchTime = 10;
DiscardedBatches = 2;
Alpha2 = 1;
Alpha3 = 1;
b2 = 150;
b3 = 5;
deltaTilde = 0.01;
epsilonTilde = 1/10^6;
BlockingSamples = 5;
MaxIterations = 20;
SuggIterations = 200;
MaxTime = 10^2;
upper_bound = 100;

FFoptObjData = NaN(MaxIterations+1, NumberSamples, NumberScenarios);
FFoptCPUData = NaN(MaxIterations+1, NumberSamples, NumberScenarios);
FFoptSolData = NaN(MaxIterations+1, NumberOfStations, NumberSamples, NumberScenarios);

SAoptObjData = NaN(MaxIterations+1, NumberSamples, NumberScenarios);
SAoptCPUData = NaN(MaxIterations+1, NumberSamples, NumberScenarios);
SAoptSolData = NaN(MaxIterations+1, NumberOfStations, NumberSamples,NumberScenarios);

SurrOptObjData = NaN(SuggIterations, NumberSamples, NumberScenarios);
SurrOptCumulativeCPUTimes = NaN(SuggIterations, NumberSamples, NumberScenarios);
SurrOptSolData = NaN(NumberOfStations, NumberSamples,NumberScenarios);

BayesOptObjData = NaN(SuggIterations, NumberSamples, NumberScenarios);
BayesOptCumulativeCPUTimes = NaN(SuggIterations, NumberSamples, NumberScenarios);
BayesOptSolData = NaN(NumberOfStations, NumberSamples,NumberScenarios);

lambdaValues = NaN(NumberScenarios,3);
alphaValues = NaN(NumberScenarios,1);
capacityCostValues = NaN(NumberScenarios,NumberOfStations);
weightValues1 = NaN(NumberScenarios,1);
weightValues2 = NaN(NumberScenarios,1);

muValues = NaN(NumberScenarios,NumberOfStations);

for m = 1:NumberScenarios
    Weights = 10*rand(2,1);
    Lambda = 2*rand*[10,20,30];
    MuTemp = 2.4*rand(1,NumberOfStations);
    Mu = [MuTemp; MuTemp];
    Alpha = 2*rand;
    CapacityCost = 0.1+0.2*rand(1,NumberOfStations);
    
    alphaValues(m,1) = Alpha;
    lambdaValues(m,:) = Lambda;
    capacityCostValues(m,:) = CapacityCost;
    weightValues1(m,:) = Weights(1,:);
    weightValues2(m,:) = Weights(2,:);

    muValues(m,:) = MuTemp;
    
    for n = 1: NumberSamples
        c = randi(100, 1, 10);
        [values1, solutions1, iterationTimes1] = FFOptSamplePathModelOne(k, c, Weights, Mu, Lambda, Alpha, CapacityCost, MaxIterations, BlockingSamples, deltaTilde, epsilonTilde, BatchTime, DiscardedBatches, MaxTime);
        FFoptCPUData(:,n,m) = iterationTimes1;
        FFoptSolData(:,:,n,m) = solutions1;
        FFoptObjData(:,n,m) = values1;

        [values2, solutions2, iterationTimes2] = StochAppSamplePathModelOne(Weights, CapacityCost, MaxIterations, BlockingSamples, c, Lambda, Mu, Alpha, NumberOfStations, BatchTime, DiscardedBatches, deltaTilde, epsilonTilde, Alpha2, Alpha3, b2, b3, MaxTime);
        SAoptCPUData(:,n,m) = iterationTimes2;
        SAoptSolData(:,:,n,m) = solutions2;
        SAoptObjData(:,n,m) = values2;
    
        lb = zeros(1,10);
        ub = upper_bound*ones(1,10);
        intcon =1:NumberOfStations;
        objconstr = @(x) -ObjectiveFunctionCrissCrossModelOne(Weights, CapacityCost, BlockingSamples, x, Lambda, Mu, Alpha, NumberOfStations, BatchTime, DiscardedBatches, epsilonTilde, deltaTilde, MaxTime);
        objconstr_BayesOpt = @(x) -ObjectiveFunctionCrissCrossModelOne(Weights, CapacityCost, BlockingSamples, [x.capacity1, x.capacity2, x.capacity3, x.capacity4, x.capacity5, x.capacity6, x.capacity7, x.capacity8, x.capacity9, x.capacity10], Lambda, Mu, Alpha, NumberOfStations, BatchTime, DiscardedBatches, epsilonTilde, deltaTilde,MaxTime);
        
        diary(strcat('SurrOutput_Model1_Scenario_',num2str(m),'Sample_',num2str(n),'.txt'))

        opts = optimoptions('surrogateopt', 'MaxFunctionEvaluations', SuggIterations, 'Display', 'iter', 'InitialPoints', c, 'PlotFcn', []);
        [x, fval, ~, output2] = surrogateopt(objconstr, lb, ub, intcon, opts)

        diary off

        SurrOptSolData(:,n,m) = x';

        temp = textread(strcat('SurrOutput_Model1_Scenario_', num2str(m),'Sample_',num2str(n),'.txt'),'%s','delimiter','\n'); %#ok<*DTXTRD>
        Index = find(contains(temp,'Surrogate'),1);
        temp = temp(8:(Index-1),1);

        Index = find(contains(temp,'F-'));

        for i = 1:length(Index)
           tempIndex = find(contains(temp,'F-'),1);
           temp = [temp(1:(tempIndex-2),1); temp((tempIndex+2):end,1)];
        end

        for i = 1:length(temp)
           foo = strsplit(cell2mat(temp(i)));
           SurrOptCumulativeCPUTimes(i,n,m) = str2double(foo(2));
           SurrOptObjData(i,n,m) = -str2double(foo(3));
        end        
                
        delete(strcat('SurrOutput_Model1_Scenario_', num2str(m),'Sample_',num2str(n),'.txt'))

        cap1 = optimizableVariable('capacity1', [0,upper_bound],'Type','integer');
        cap2 = optimizableVariable('capacity2', [0,upper_bound],'Type','integer');
        cap3 = optimizableVariable('capacity3', [0,upper_bound],'Type','integer');
        cap4 = optimizableVariable('capacity4', [0,upper_bound],'Type','integer');
        cap5 = optimizableVariable('capacity5', [0,upper_bound],'Type','integer');
        cap6 = optimizableVariable('capacity6', [0,upper_bound],'Type','integer');
        cap7 = optimizableVariable('capacity7', [0,upper_bound],'Type','integer');
        cap8 = optimizableVariable('capacity8', [0,upper_bound],'Type','integer');
        cap9 = optimizableVariable('capacity9', [0,upper_bound],'Type','integer');
        cap10 = optimizableVariable('capacity10', [0,upper_bound],'Type','integer');

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

        BayesOptCumulativeCPUTimes(:,n,m) = cumsum(results.IterationTimeTrace);
        BayesOptObjData(:,n,m) = -results.ObjectiveMinimumTrace';
        BayesOptSolData(:,n,m) = [results.XAtMinObjective.capacity1; results.XAtMinObjective.capacity2; ...
                                results.XAtMinObjective.capacity3; results.XAtMinObjective.capacity4; ...
                                results.XAtMinObjective.capacity5; results.XAtMinObjective.capacity6; ...
                                results.XAtMinObjective.capacity7; results.XAtMinObjective.capacity8; ...
                                results.XAtMinObjective.capacity9; results.XAtMinObjective.capacity10 ...
                                ];
        close all;    
        
        
        disp('Progress:')
        (NumberSamples*(m-1)+n)/(NumberSamples*NumberScenarios)
    end
end

save('ModelOneSimulationData.mat')


matrix2latex([CapacityCost; Weights], 'out2.tex')