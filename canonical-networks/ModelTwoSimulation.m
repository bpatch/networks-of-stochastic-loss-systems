clear all; clc; close all;
rng(2020)
NumberSamples = 10; 
NumberScenarios = 20;
NumberOfStations = 6;
k = 2;
BatchTime = 10;
BurnInTime = 10;
Alpha2 = 1;
Alpha3 = 1;
b2 = 150;
b3 = 5;
deltaTilde = 0.01;
epsilonTilde = 1/10^6;
BlockingSamples = 10;
MaxIterations = 30;
SuggIterations = 200;
MaxTime = 10^2;
upper_bound = 100;

FFoptObjData = NaN(MaxIterations+1, NumberSamples, NumberScenarios);
FFoptCPUData = NaN(MaxIterations+1, NumberSamples, NumberScenarios);
FFoptSolData = NaN(MaxIterations+1, NumberOfStations, NumberSamples, NumberScenarios);

SurrOptObjData = NaN(SuggIterations, NumberSamples, NumberScenarios);
SurrOptCumulativeCPUTimes = NaN(SuggIterations, NumberSamples, NumberScenarios);
SurrOptSolData = NaN(NumberOfStations, NumberSamples,NumberScenarios);

BayesOptObjData = NaN(SuggIterations, NumberSamples, NumberScenarios);
BayesOptCumulativeCPUTimes = NaN(SuggIterations, NumberSamples, NumberScenarios);
BayesOptSolData = NaN(NumberOfStations, NumberSamples,NumberScenarios);

lambdaValues = NaN(NumberScenarios,2);
weightValues = NaN(12,3,NumberScenarios);
muValues = NaN(NumberScenarios,6);
capacityCostValues = NaN(NumberScenarios,6);


for m = 1:NumberScenarios
    Weights = 50*rand(12,3);
    Lambda = 0.2*rand(1,2);
    Mu = 0.2*rand(1,6);
    CapacityCost = 0.1+0.2*rand(1,6);
    
    lambdaValues(m,:) = Lambda;
    weightValues(:,:,m) = Weights;
    muValues(m,:) = Mu;
    capacityCostValues(m,:) = CapacityCost;
    
    for n = 1: NumberSamples
        c = randi(100, 1, 6);
        [values1, solutions1, iterationTimes1] = FFOptSamplePathModelTwo(k, c, Weights, Mu, Lambda, CapacityCost, MaxIterations, BlockingSamples, deltaTilde, epsilonTilde, BatchTime, BurnInTime,MaxTime,m);
        FFoptCPUData(:,n,m) = iterationTimes1;
        FFoptSolData(:,:,n,m) = solutions1;
        FFoptObjData(:,n,m) = values1;
        
        
        lb = zeros(1,6);
        ub = upper_bound*ones(1,6);
        intcon =1:6; 
        objconstr = @(x) -ObjectiveFunctionCanonicalModelTwo(Weights, BlockingSamples, x, Mu,  Lambda, CapacityCost, BatchTime, BurnInTime, epsilonTilde, deltaTilde,MaxTime,m);
        objconstr_BayesOpt = @(x) -ObjectiveFunctionCanonicalModelTwo(Weights, BlockingSamples, [x.capacity1, x.capacity2, x.capacity3, x.capacity4, x.capacity5, x.capacity6], Mu,  Lambda, CapacityCost, BatchTime, BurnInTime, epsilonTilde, deltaTilde,MaxTime,m);

        diary(strcat('SurrOutput_Model2_Scenario_',num2str(m),'Sample_',num2str(n),'.txt'))

        opts = optimoptions('surrogateopt', 'MaxFunctionEvaluations', SuggIterations, 'Display', 'iter', 'InitialPoints', c, 'PlotFcn', []);
        [x, fval, ~, output2] = surrogateopt(objconstr, lb, ub, intcon, opts)

        diary off
       SurrOptSolData(:,n,m) = x';

        temp = textread(strcat('SurrOutput_Model2_Scenario_', num2str(m),'Sample_',num2str(n),'.txt'),'%s','delimiter','\n'); %#ok<*DTXTRD>
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
                
        delete(strcat('SurrOutput_Model2_Scenario_', num2str(m),'Sample_',num2str(n),'.txt'))
        
        cap1 = optimizableVariable('capacity1', [0,upper_bound],'Type','integer');
        cap2 = optimizableVariable('capacity2', [0,upper_bound],'Type','integer');
        cap3 = optimizableVariable('capacity3', [0,upper_bound],'Type','integer');
        cap4 = optimizableVariable('capacity4', [0,upper_bound],'Type','integer');
        cap5 = optimizableVariable('capacity5', [0,upper_bound],'Type','integer');
        cap6 = optimizableVariable('capacity6', [0,upper_bound],'Type','integer');

        capacity1 = c(1,1);
        capacity2 = c(1,2);
        capacity3 = c(1,3);
        capacity4 = c(1,4);
        capacity5 = c(1,5);
        capacity6 = c(1,6);

        C0_BayesOpt = table(capacity1, capacity2,capacity3, capacity4,capacity5, capacity6);

        results = bayesopt(objconstr_BayesOpt,[cap1, cap2,cap3, cap4,cap5, cap6],'Verbose',0,...
        'AcquisitionFunctionName','expected-improvement-plus','InitialX',C0_BayesOpt,'MaxObjectiveEvaluations',SuggIterations);

        BayesOptCumulativeCPUTimes(:,n,m) = cumsum(results.IterationTimeTrace);
        BayesOptObjData(:,n,m) = -results.ObjectiveMinimumTrace';
        BayesOptSolData(:,n,m) = [results.XAtMinObjective.capacity1; results.XAtMinObjective.capacity2; ...
                                results.XAtMinObjective.capacity3; results.XAtMinObjective.capacity4; ...
                                results.XAtMinObjective.capacity5; results.XAtMinObjective.capacity6; ...
                                ];
        close all;    
        
        
        disp('Progress:')

        (NumberSamples*(m-1)+n)/(NumberSamples*NumberScenarios)
    end
end

save('ModelTwoSimulationData.mat')

