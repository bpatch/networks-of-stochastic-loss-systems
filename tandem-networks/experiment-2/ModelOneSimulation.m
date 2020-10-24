clear all; clc; close all;
rng(2020);
NumberSamples = 10;
NumberScenarios = 20;
deltaTilde = 0.1;
epsilonTilde = 1/10^6;
BatchTime = 10;
DiscardedBatches = 2; 
BlockingSamples = 5;
MaxIterations = 10;
SuggIterations = 100;
k=2;
c1max = 60;
c2max = 60;

Alpha2 = 1;
Alpha3 = 1;
b2 = 150;
b3 = 5;

lb = zeros(1,2);
ub = 100*ones(1,2);
intcon =1:2;

MAMoptData = NaN(NumberScenarios,1);
MaMsolData = NaN(NumberScenarios,2);

FFoptObjData = NaN(NumberScenarios, NumberSamples);
FFoptCPUData = NaN(NumberScenarios, NumberSamples);
FFoptSolData = NaN(NumberScenarios, 2, NumberSamples);

SAoptObjData = NaN(NumberScenarios, NumberSamples);
SAoptCPUData = NaN(NumberScenarios, NumberSamples);
SAoptSolData = NaN(NumberScenarios, 2, NumberSamples);

SurrOptObjData = NaN(SuggIterations, NumberSamples, NumberScenarios);
SurrOptCumulativeCPUTimes = NaN(SuggIterations, NumberSamples, NumberScenarios);
SurrOptSolData = NaN(NumberScenarios, 2, NumberSamples);
SurrOptCPUData = NaN(NumberScenarios, NumberSamples);

BayesOptObjData = NaN(SuggIterations, NumberSamples, NumberScenarios);
BayesOptCumulativeCPUTimes = NaN(SuggIterations, NumberSamples, NumberScenarios);
BayesOptSolData = NaN(NumberScenarios, 2, NumberSamples);
BayesOptCPUData = NaN(NumberScenarios, NumberSamples);


InitialSolData= NaN(NumberScenarios, 2, NumberSamples);
InitialObjData= NaN(NumberScenarios, NumberSamples);

lambdaValues = NaN(NumberScenarios,1);
capacityCostValues = NaN(NumberScenarios,2);
weightValues = NaN(NumberScenarios,1);
muValues = NaN(NumberScenarios,2);


cap1 = optimizableVariable('capacity1', [0,100],'Type','integer');
cap2 = optimizableVariable('capacity2', [0,100],'Type','integer');

for m = 1:NumberScenarios
    c1=c1max;
    c2=c2max;
    RealMax = -1;
    while (abs(c1-c1max) < 2 || abs(c2-c2max) < 2)
        disp('Currently looking for a new model')
        lambda = 20*rand;
        CapacityCost = 0.4*rand(2,1);
        Weights = 5*rand;
        Mu = 3*rand(2,1); 
        [c1, c2, RealMax] = ExactSolutionModelOne(c1max, c2max, Weights, Mu, lambda, CapacityCost);
    end
    disp('Model found, evaluating algorithms!')
    lambdaValues(m,1) = lambda;
    capacityCostValues(m,:) = CapacityCost';
    weightValues(m,1) = Weights;
    muValues(m,:) = Mu';

    MAMoptData(m,1) = RealMax;
    MaMsolData(m,:) = [c1,c2];

    for n = 1:NumberSamples
        c = randi(100, 2, 1);
        
        InitialSolData(m,:,n) = c';
        
        [values1, solutions1, iterationTimes1] = FFOptSamplePathModelOne(k, c, Weights, Mu, lambda, CapacityCost, MaxIterations, BlockingSamples, deltaTilde, epsilonTilde, BatchTime, DiscardedBatches);
        [M, I] = max(values1);
        if M > 0
            I = find(values1 >= 0.99*M,1);
        end
        FFoptCPUData(m, n) = sum(iterationTimes1(1:(I-1)));
        FFoptSolData(m,:,n) = solutions1(I,:);

        [values2, solutions2, iterationTimes2] = StochAppSamplePathModelOne(Weights, CapacityCost, MaxIterations, BlockingSamples, c, lambda, Mu, BatchTime, DiscardedBatches, deltaTilde, epsilonTilde, Alpha2, Alpha3, b2, b3);
        [M, I] = max(values2);
        if M > 0
            I = find(values2 >= 0.99*M,1);
        end
        SAoptCPUData(m, n) = sum(iterationTimes2(1:(I)));
        SAoptSolData(m,:,n) = solutions2(I,:);
      
        objconstr = @(x) -ObjectiveFunctionTandemModelOne(Weights, CapacityCost, BlockingSamples, x', lambda, Mu, BatchTime, DiscardedBatches, epsilonTilde, deltaTilde);

        diary(strcat('SurrOutput_Model1_Scenario_',num2str(m),'Sample_',num2str(n),'.txt'))

        opts = optimoptions('surrogateopt', 'MaxFunctionEvaluations', 100, 'Display', 'iter', 'InitialPoints', c, 'PlotFcn', []);
        [x, fval, ~, output2] = surrogateopt(objconstr, lb, ub, intcon, opts)

        diary off        
        
        SurrOptSolData(m,:,n) = x;

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
        
        values3 = SurrOptObjData(:,n,m);
        [M, I] = max(values3);
        if M > 0
            I = find(values3 >= 0.99*M,1);
        end
        SurrOptCPUData(m, n) = SurrOptCumulativeCPUTimes(I,n,m);
        
        delete(strcat('SurrOutput_Model1_Scenario_', num2str(m),'Sample_',num2str(n),'.txt'))
                
        objconstr_BayesOpt = @(x) -ObjectiveFunctionTandemModelOne(Weights, CapacityCost, BlockingSamples, [x.capacity1; x.capacity2], lambda, Mu, BatchTime, DiscardedBatches, epsilonTilde, deltaTilde);
        
        capacity1 = c(1,1);
        capacity2 = c(2,1);

        C0_BayesOpt = table(capacity1, capacity2);

        results = bayesopt(objconstr_BayesOpt,[cap1, cap2],'Verbose',0,...
        'AcquisitionFunctionName','expected-improvement-plus','InitialX',C0_BayesOpt,'MaxObjectiveEvaluations',SuggIterations)
        close all;

        BayesOptCumulativeCPUTimes(:,n,m) = cumsum(results.IterationTimeTrace);
        
        values4 = -results.ObjectiveMinimumTrace;
        BayesOptObjData(:,n,m) = values4;
        BayesOptSolData(m,:,n) = [results.XAtMinObjective.capacity1, results.XAtMinObjective.capacity2];
        
        
        [M, I] = max(values4);
        if M > 0
            I = find(values4 >= 0.99*M,1);
        end
        BayesOptCPUData(m, n) = BayesOptCumulativeCPUTimes(I,n,m);
        
        (NumberSamples*(m-1)+n)/(NumberSamples*NumberScenarios)
    end
end


save('ModelOneSimulationData.mat')

matrix = [(1:20)', lambdaValues, muValues, capacityCostValues, weightValues];
matrix2latex(matrix, 'out1.tex')