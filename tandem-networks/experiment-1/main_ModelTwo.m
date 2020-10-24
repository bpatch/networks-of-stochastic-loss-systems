clear all; clc; close all;
rng(2020);
NumberSamples = 10;
deltaTilde = 0.1; %0.01;
epsilonTilde = 1/10^6;
BatchTime = 10;
DiscardedBatches = 2; 
BlockingSamples = 5;
MaxIterations = 10;
SuggIterations = 25;
k=2;
c1max = 40;
c2max = 40;

Alpha2 = 1;
Alpha3 = 1;
b2 = 150;
b3 = 5;

lb = zeros(1,2);
ub = 100*ones(1,2);
intcon =1:2;

lambda = 16;
CapacityCost =[0.2; 0.3];
Weights = [1;0.9];
Mu = [0.8;0.6]; 
[c1, c2, RealMax] = ExactSolutionModelTwo(c1max, c2max, Weights, Mu, lambda, CapacityCost);

objconstr = @(x) -ObjectiveFunctionTandemModelTwo(Weights, CapacityCost, BlockingSamples, x', lambda, Mu, BatchTime, DiscardedBatches, epsilonTilde, deltaTilde);
objconstr_BayesOpt = @(x) -ObjectiveFunctionTandemModelTwo(Weights, CapacityCost, BlockingSamples, [x.capacity1; x.capacity2], lambda, Mu, BatchTime, DiscardedBatches, epsilonTilde, deltaTilde);

FFoptObjData = NaN(NumberSamples, MaxIterations+1);
FFoptCPUData = NaN(NumberSamples, MaxIterations+1);
FFoptSolData = NaN(MaxIterations+1,2,NumberSamples);

SAoptObjData = NaN(NumberSamples, MaxIterations+1);
SAoptCPUData = NaN(NumberSamples, MaxIterations+1);
SAoptSolData = NaN(MaxIterations+1,2,NumberSamples);

SurrOptObjData = NaN(NumberSamples, SuggIterations);
SurrOptCPUData = NaN(NumberSamples, SuggIterations);
SurrOptSolData = NaN(NumberSamples,2);

BayesOptObjData = NaN(NumberSamples, SuggIterations);
BayesOptCPUData = NaN(NumberSamples, SuggIterations);
BayesOptSolData = NaN(NumberSamples,2);


for n = 1:NumberSamples
    counter = 1;
    C0 = randi(100, 2, 1);
    [values1, solutions1, iterationTimes1] = FFOptSamplePathModelTwo(k, C0, Weights, Mu, lambda, CapacityCost, MaxIterations, BlockingSamples, deltaTilde, epsilonTilde, BatchTime, DiscardedBatches);
    FFoptObjData(n, :) = values1';
    FFoptCPUData(n, :) = cumsum([0;iterationTimes1])';
    FFoptSolData(:,:,n) = solutions1;

    [values2, solutions2, iterationTimes2] = StochAppSamplePathModelTwo(Weights, CapacityCost, MaxIterations, BlockingSamples, C0, lambda, Mu, BatchTime, DiscardedBatches, deltaTilde, epsilonTilde, Alpha2, Alpha3, b2, b3);
    SAoptObjData(n, :) = values2';
    SAoptCPUData(n, :) = cumsum(iterationTimes2)';
    SAoptSolData(:,:,n) = solutions2;
    
    diary(strcat('SuggOutput_Model2_Sample',num2str(n),'.txt'))

    opts = optimoptions('surrogateopt','MaxFunctionEvaluations',SuggIterations,'Display','iter','InitialPoints',C0, 'PlotFcn',[]);
    [x, fval, ~, output2] = surrogateopt(objconstr, lb, ub, intcon, opts) 

    diary off
    
    SurrOptSolData(n,:) = x;

    temp = textread(strcat('SuggOutput_Model2_Sample',num2str(n),'.txt'),'%s','delimiter','\n'); %#ok<*DTXTRD>
    Index = find(contains(temp,'Surrogate'),1);
    temp = temp(8:(Index-1),1);

    Index = find(contains(temp,'F-'));

    for i = 1:length(Index)
       tempIndex = find(contains(temp,'F-'),1);
       temp = [temp(1:(tempIndex-2),1); temp((tempIndex+2):end,1)];
    end

    for i = 1:length(temp)
       foo = strsplit(cell2mat(temp(i)));
       SurrOptCPUData(n,i) = str2double(foo(2));
       SurrOptObjData(n,i) = -str2double(foo(3));
    end

    delete(strcat('SuggOutput_Model2_Sample',num2str(n),'.txt'))
    
    cap1 = optimizableVariable('capacity1', [0,100],'Type','integer');
    cap2 = optimizableVariable('capacity2', [0,100],'Type','integer');

    capacity1 = C0(1,1);
    capacity2 = C0(2,1);
    
    C0_BayesOpt = table(capacity1, capacity2);
    
    results = bayesopt(objconstr_BayesOpt,[cap1, cap2],'Verbose',0,...
    'AcquisitionFunctionName','expected-improvement-plus','InitialX',C0_BayesOpt,'MaxObjectiveEvaluations',SuggIterations)

    BayesOptCPUData(n,:) = cumsum(results.IterationTimeTrace');
    BayesOptObjData(n,:) = -results.ObjectiveMinimumTrace';
    BayesOptSolData(n,:) = [results.XAtMinObjective.capacity1, results.XAtMinObjective.capacity2];
    close all;
    
    
    disp('Progress:')
    
    n/NumberSamples
end
SurrOptCPUData = [zeros(NumberSamples,1), SurrOptCPUData];
BayesOptCPUData = [zeros(NumberSamples,1), BayesOptCPUData];
save('ModelTwo.mat')


load('ModelTwo.mat')


figure()
plot(0:(MaxIterations),max(FFoptObjData',[],2), 'color', [0, 0.4470, 0.7410],'linestyle', 'none')
hold on
plot(0:(MaxIterations),min(FFoptObjData',[],2), 'color',[0, 0.4470, 0.7410],'linestyle', 'none')
x2 = [0:(MaxIterations), fliplr(0:(MaxIterations))];
inBetween = [min(FFoptObjData',[],2)',fliplr(max(FFoptObjData',[],2)')];
fill(x2, inBetween,[0, 0.4470, 0.7410],'linestyle', 'none')
alpha(0.15)
hold on
plot(0:(MaxIterations), FFoptObjData','color',[0.75 0.75 0.75],'linewidth', 0.2)
plot(0:(MaxIterations), mean(FFoptObjData', 2),'color',[0, 0.4470, 0.7410],'linewidth', 0.5)
plot(0:(MaxIterations),RealMax*ones(1,MaxIterations+1),'k--','linewidth', 0.5)
xlabel('Iteration')
ylabel('Obj. Val.')
title('Alg. 1') 
%axis([IntervalK MaxK 15 RealMax+0.5])
box off


figure()
plot(0:(MaxIterations),max(SAoptObjData',[],2), 'color', [0, 0.4470, 0.7410],'linestyle', 'none')
hold on
plot(0:(MaxIterations),min(SAoptObjData',[],2), 'color',[0, 0.4470, 0.7410],'linestyle', 'none')
x2 = [0:(MaxIterations), fliplr(0:(MaxIterations))];
inBetween = [min(SAoptObjData',[],2)',fliplr(max(SAoptObjData',[],2)')];
fill(x2, inBetween,[0, 0.4470, 0.7410],'linestyle', 'none')
alpha(0.15)
hold on
plot(0:(MaxIterations), SAoptObjData','color',[0.75 0.75 0.75],'linewidth', 0.2)
plot(0:(MaxIterations), mean(SAoptObjData', 2),'color',[0, 0.4470, 0.7410],'linewidth', 0.5)
plot(0:(MaxIterations),RealMax*ones(1,MaxIterations+1),'k--','linewidth', 0.5)
xlabel('Iteration')
ylabel('Obj. Val.')
title('Stoch. App.') 
%axis([IntervalK MaxK 15 RealMax+0.5])
box off

figure()
plot(0:(SuggIterations-1),max(SurrOptObjData',[],2), 'color', [0, 0.4470, 0.7410],'linestyle', 'none')
hold on
plot(0:(SuggIterations-1),min(SurrOptObjData',[],2), 'color',[0, 0.4470, 0.7410],'linestyle', 'none')
x2 = [0:(SuggIterations-1), fliplr(0:(SuggIterations-1))];
inBetween = [min(SurrOptObjData',[],2)',fliplr(max(SurrOptObjData',[],2)')];
fill(x2, inBetween,[0, 0.4470, 0.7410],'linestyle', 'none')
alpha(0.15)
hold on
plot(0:(SuggIterations-1), SurrOptObjData','color',[0.75 0.75 0.75],'linewidth', 0.2)
plot(0:(SuggIterations-1), mean(SurrOptObjData', 2),'color',[0, 0.4470, 0.7410],'linewidth', 0.5)
plot(0:(SuggIterations-1),RealMax*ones(1,SuggIterations-1+1),'k--','linewidth', 0.5)
xlabel('Iteration')
ylabel('Obj. Val.')
title('Surr. Opt.') 
%axis([IntervalK MaxK 15 RealMax+0.5])
box off

figure()
plot(0:(SuggIterations-1),max(BayesOptObjData',[],2), 'color', [0, 0.4470, 0.7410],'linestyle', 'none')
hold on
plot(0:(SuggIterations-1),min(BayesOptObjData',[],2), 'color',[0, 0.4470, 0.7410],'linestyle', 'none')
x2 = [0:(SuggIterations-1), fliplr(0:(SuggIterations-1))];
inBetween = [min(BayesOptObjData',[],2)',fliplr(max(BayesOptObjData',[],2)')];
fill(x2, inBetween,[0, 0.4470, 0.7410],'linestyle', 'none')
alpha(0.15)
hold on
plot(0:(SuggIterations-1), BayesOptObjData','color',[0.75 0.75 0.75],'linewidth', 0.2)
plot(0:(SuggIterations-1), mean(BayesOptObjData', 2),'color',[0, 0.4470, 0.7410],'linewidth', 0.5)
plot(0:(SuggIterations-1),RealMax*ones(1,SuggIterations-1+1),'k--','linewidth', 0.5)
xlabel('Iteration')
ylabel('Obj. Val.')
title('Bayes Opt.') 
%axis([IntervalK MaxK 15 RealMax+0.5])
box off

figure()
%subplot(1,3,2)
plot(0:(MaxIterations),max(cumsum(FFoptCPUData),[],1), 'color', [0, 0.4470, 0.7410],'linestyle', 'none')
hold on
plot(0:(MaxIterations),min(cumsum(FFoptCPUData),[],1), 'color',[0, 0.4470, 0.7410],'linestyle', 'none')
x2 = [0:(MaxIterations), fliplr(0:(MaxIterations))];
inBetween = [min(cumsum(FFoptCPUData),[],1),fliplr(max(cumsum(FFoptCPUData),[],1))];
fill(x2, inBetween,[0, 0.4470, 0.7410],'linestyle', 'none')
alpha(0.15)
hold on
plot(0:(MaxIterations), cumsum(FFoptCPUData),'color',[0.75 0.75 0.75],'linewidth', 0.25)
plot(0:(MaxIterations), mean(cumsum(FFoptCPUData),1),'color',[0, 0.4470, 0.7410],'linewidth', 2)
xlabel('Iteration')
ylabel('Cumulative CPU time')
title('Alg. 1')
%axis([IntervalK MaxK 15 RealMax+0.5])
box off

figure()
%subplot(1,3,3)
plot(0:(MaxIterations),max(cumsum(SAoptCPUData),[],1), 'color', [0, 0.4470, 0.7410],'linestyle', 'none')
hold on
plot(0:(MaxIterations),min(cumsum(SAoptCPUData),[],1), 'color',[0, 0.4470, 0.7410],'linestyle', 'none')
x2 = [0:(MaxIterations), fliplr(0:(MaxIterations))];
inBetween = [min(cumsum(SAoptCPUData),[],1),fliplr(max(cumsum(SAoptCPUData),[],1))];
fill(x2, inBetween,[0, 0.4470, 0.7410],'linestyle', 'none')
alpha(0.15)
hold on
plot(0:(MaxIterations), cumsum(SAoptCPUData),'color',[0.75 0.75 0.75],'linewidth', 0.25)
plot(0:(MaxIterations), mean(cumsum(SAoptCPUData),1),'color',[0, 0.4470, 0.7410],'linewidth', 2)
xlabel('Iteration')
ylabel('CPU (secs)')
title('Stoch. App.')
%axis([IntervalK MaxK 15 RealMax+0.5])
box off


figure()
%subplot(1,3,3)
plot(0:SuggIterations,max(SurrOptCPUData,[],1), 'color', [0, 0.4470, 0.7410],'linestyle', 'none')
hold on
plot(0:SuggIterations,min(SurrOptCPUData,[],1), 'color',[0, 0.4470, 0.7410],'linestyle', 'none')
x2 = [0:SuggIterations, fliplr(0:SuggIterations)];
inBetween = [min(SurrOptCPUData,[],1),fliplr(max(SurrOptCPUData,[],1))];
fill(x2, inBetween,[0, 0.4470, 0.7410],'linestyle', 'none')
alpha(0.15)
hold on
plot(0:SuggIterations, SurrOptCPUData,'color',[0.75 0.75 0.75],'linewidth', 0.25)
plot(0:SuggIterations, mean(SurrOptCPUData,1),'color',[0, 0.4470, 0.7410],'linewidth', 2)
xlabel('Iteration')
ylabel('CPU (secs)')
title('Surr. Opt.')
%axis([IntervalK MaxK 15 RealMax+0.5])
box off

figure()
%subplot(1,3,3)
plot(0:SuggIterations,max(BayesOptCPUData,[],1), 'color', [0, 0.4470, 0.7410],'linestyle', 'none')
hold on
plot(0:SuggIterations,min(BayesOptCPUData,[],1), 'color',[0, 0.4470, 0.7410],'linestyle', 'none')
x2 = [0:SuggIterations, fliplr(0:SuggIterations)];
inBetween = [min(BayesOptCPUData,[],1),fliplr(max(BayesOptCPUData,[],1))];
fill(x2, inBetween,[0, 0.4470, 0.7410],'linestyle', 'none')
alpha(0.15)
hold on
plot(0:SuggIterations, BayesOptCPUData,'color',[0.75 0.75 0.75],'linewidth', 0.25)
plot(0:SuggIterations, mean(BayesOptCPUData,1),'color',[0, 0.4470, 0.7410],'linewidth', 2)
xlabel('Iteration')
ylabel('CPU (secs)')
title('Bayes Opt.')
%axis([IntervalK MaxK 15 RealMax+0.5])
box off

% % save('ModelTwo.mat')