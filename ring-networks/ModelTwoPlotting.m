clear all; close all; clc;
% load('ModelTwoSimulationData.mat')
% 
% BatchTime = 10; 
% DiscardedBatches = 2;
% deltaTilde = 0.01;
% MaxTime = 10^3;
% 
% SAoptObjData2 = NaN(NumberScenarios, NumberSamples);
% SAoptCPUData2 = NaN(NumberScenarios, NumberSamples);
% SAoptSolData2 = NaN(NumberScenarios, NumberSamples, NumberOfStations);
% 
% FFoptObjData2 = NaN(NumberScenarios, NumberSamples);
% FFoptCPUData2 = NaN(NumberScenarios, NumberSamples);
% FFoptSolData2 = NaN(NumberScenarios, NumberSamples, NumberOfStations);
% 
% SurrOptObjData2 = NaN(NumberScenarios, NumberSamples);
% SurrOptCPUData2 = NaN(NumberScenarios, NumberSamples);
% SurrOptSolData2 = NaN(NumberScenarios, NumberSamples, NumberOfStations);
% 
% BayesOptObjData2 = NaN(NumberScenarios, NumberSamples);
% BayesOptCPUData2 = NaN(NumberScenarios, NumberSamples);
% BayesOptSolData2 = NaN(NumberScenarios, NumberSamples, NumberOfStations);
% 
% 
% InitialObjData = NaN(NumberScenarios, NumberSamples);
% 
% for m = 1:NumberScenarios
%     for n = 1:NumberSamples
%         values = SAoptObjData(:,n,m);
%         cpuTimes = cumsum(SAoptCPUData(:,n,m));
%         [M, I] = max(values);
%         if M > 0
%             I = find(values >= 0.99*M,1);
%         end
%         SAoptCPUData2(m,n) = cpuTimes(I);
%         SAoptSolData2(m,n,:) = SAoptSolData(I,:,n,m);
%         rng(2020)
%         SAoptObjData2(m, n) = ObjectiveFunctionCrissCrossModelTwo([weightValues1(m,:); weightValues2(m,:)], capacityCostValues(m,:), BlockingSamples, SAoptSolData(I,:,n,m), lambdaValues(m,:), muValues(m,:), alphaValues(m,1), NumberOfStations, BatchTime, DiscardedBatches, epsilonTilde, deltaTilde,MaxTime);
%     end
% end
% 
% for m = 1:NumberScenarios
%     for n = 1:NumberSamples
%         values = FFoptObjData(:,n,m);
%         cpuTimes = cumsum(FFoptCPUData(:,n,m));
%         [M, I] = max(values);
%         if M > 0
%             I = find(values >= 0.99*M,1);
%         end
%         FFoptCPUData2(m,n) = cpuTimes(I);
%         FFoptSolData2(m,n,:) = FFoptSolData(I,:,n,m);
%         rng(2020)
%         FFoptObjData2(m, n) = ObjectiveFunctionCrissCrossModelTwo([weightValues1(m,:); weightValues2(m,:)], capacityCostValues(m,:), BlockingSamples, FFoptSolData(I,:,n,m), lambdaValues(m,:), muValues(m,:), alphaValues(m,1), NumberOfStations, BatchTime, DiscardedBatches, epsilonTilde, deltaTilde,MaxTime);
%     end
% end
% 
% for m = 1:NumberScenarios
%     for n = 1:NumberSamples
%         values = SurrOptObjData(:,n,m);
%         cpuTimes = SurrOptCumulativeCPUTimes(:,n,m);
%         [M, I] = max(values);
%         if M > 0
%             I = find(values >= 0.99*M,1);
%         end
%         SurrOptCPUData2(m,n) = cpuTimes(I);
%         SurrOptSolData2(m,n,:) = SurrOptSolData(:,n,m);
%         rng(2020)
%         SurrOptObjData2(m, n) = ObjectiveFunctionCrissCrossModelTwo([weightValues1(m,:); weightValues2(m,:)], capacityCostValues(m,:), BlockingSamples, SurrOptSolData(:,n,m)', lambdaValues(m,:), muValues(m,:), alphaValues(m,1), NumberOfStations, BatchTime, DiscardedBatches, epsilonTilde, deltaTilde,MaxTime);
%     end
% end
% 
% for m = 1:NumberScenarios
%     for n = 1:NumberSamples
%         values = BayesOptObjData(:,n,m);
%         cpuTimes = BayesOptCumulativeCPUTimes(:,n,m);
%         [M, I] = max(values);
%         if M > 0
%             I = find(values >= 0.99*M,1);
%         end
%         BayesOptCPUData2(m,n) = cpuTimes(I);
%         BayesOptSolData2(m,n,:) = BayesOptSolData(:,n,m);
%         rng(2020)
%         BayesOptObjData2(m, n) = ObjectiveFunctionCrissCrossModelTwo([weightValues1(m,:); weightValues2(m,:)], capacityCostValues(m,:), BlockingSamples, BayesOptSolData(:,n,m)', lambdaValues(m,:), muValues(m,:), alphaValues(m,1), NumberOfStations, BatchTime, DiscardedBatches, epsilonTilde, deltaTilde,MaxTime);
%     end
% end
% 
% for m = 1:NumberScenarios
%     for n = 1:NumberSamples
%         rng(2020)
%         InitialObjData(m, n) = ObjectiveFunctionCrissCrossModelTwo([weightValues1(m,:); weightValues2(m,:)], capacityCostValues(m,:), BlockingSamples, FFoptSolData(1,:,n,m), lambdaValues(m,:), muValues(m,:), alphaValues(m,1), NumberOfStations, BatchTime, DiscardedBatches, epsilonTilde, deltaTilde,MaxTime);
%     end
% end
% 
% save('ModelTwoPlottingData.mat')

load('ModelTwoPlottingData.mat')

figure()
hold on
plot(1:(NumberScenarios),max(FFoptObjData2,[],2), 'color', [0, 0.5, 0],'linestyle', 'none')
plot(1:(NumberScenarios),min(FFoptObjData2,[],2), 'color',[0, 0.5, 0],'linestyle', 'none')
x2 = [1:(NumberScenarios), fliplr(1:(NumberScenarios))];
inBetween = [min(FFoptObjData2,[],2)',fliplr(max(FFoptObjData2,[],2)')];
fill(x2, inBetween,[0, 0.5, 0],'linestyle', 'none')
alpha(0.15)
plot(1:(NumberScenarios), mean(FFoptObjData2,2),'color',[0, 0.5, 0],'linewidth', 0.2)

plot(1:(NumberScenarios),min(SAoptObjData2,[],2), 'color',[0.83, 0.13, 0.18],'linestyle', 'none')
x2 = [1:(NumberScenarios), fliplr(1:(NumberScenarios))];
inBetween = [min(SAoptObjData2,[],2)',fliplr(max(SAoptObjData2,[],2)')];
fill(x2, inBetween,[0.83, 0.13, 0.18],'linestyle', 'none')
alpha(0.15)
plot(1:(NumberScenarios), mean(SAoptObjData2,2),'color',[0.8, 0.33, 0],'linewidth', 0.2,'linestyle',':')

plot(1:(NumberScenarios),min(SurrOptObjData2,[],2), 'color',[0.8, 0.33, 0],'linestyle', 'none')
x2 = [1:(NumberScenarios), fliplr(1:(NumberScenarios))];
inBetween = [min(SurrOptObjData2,[],2)',fliplr(max(SurrOptObjData2,[],2)')];
fill(x2, inBetween,[0.8, 0.33, 0],'linestyle', 'none')
alpha(0.15)
plot(1:(NumberScenarios), mean(SurrOptObjData2,2),'color',[0.8, 0.33, 0],'linewidth', 0.2,'linestyle','--')
xlabel('Scenario')
ylabel('Obj. Val.')

plot(1:(NumberScenarios),min(InitialObjData,[],2), 'color',[0, 0.28, 0.73],'linestyle', 'none')
x2 = [1:(NumberScenarios), fliplr(1:(NumberScenarios))];
inBetween = [min(InitialObjData,[],2)',fliplr(max(InitialObjData,[],2)')];
fill(x2, inBetween,[0, 0.28, 0.73],'linestyle', 'none')
alpha(0.15)
plot(1:(NumberScenarios), mean(InitialObjData,2),'color',[0, 0.28, 0.73],'linewidth', 0.2,'linestyle','--')

plot(1:(NumberScenarios),min(BayesOptObjData2,[],2), 'color',[0, 0.55, 0.55],'linestyle', 'none')
x2 = [1:(NumberScenarios), fliplr(1:(NumberScenarios))];
inBetween = [min(BayesOptObjData2,[],2)',fliplr(max(BayesOptObjData2,[],2)')];
fill(x2, inBetween,[0, 0.55, 0.55],'linestyle', 'none')
alpha(0.15)
plot(1:(NumberScenarios), mean(BayesOptObjData2,2),'color',[0, 0.55, 0.55],'linewidth', 0.2,'linestyle','-.')
xlabel('Scenario')
ylabel('Obj. Val.')


ax = gca;
ax.YAxis.MinorTick = 'off';
ax.XGrid = 'on';
ax.YGrid = 'on';
box off; 

figure()
plot(1:(NumberScenarios),max(FFoptCPUData2,[],2), 'color', [0, 0.5, 0],'linestyle', 'none')
hold on
plot(1:(NumberScenarios),min(FFoptCPUData2,[],2), 'color',[0, 0.5, 0],'linestyle', 'none')
x2 = [1:(NumberScenarios), fliplr(1:(NumberScenarios))];
inBetween = max([min(FFoptCPUData2,[],2)',fliplr(max(FFoptCPUData2,[],2)')],0.001);
fill(x2, inBetween,[0, 0.5, 0],'linestyle', 'none')
alpha(0.15)
plot(1:(NumberScenarios), mean(FFoptCPUData2,2),'color',[0, 0.5, 0],'linewidth', 0.2)

plot(1:(NumberScenarios),min(max(SAoptCPUData2,0.01),[],2), 'color',[0.83, 0.13, 0.18],'linestyle', 'none')
x2 = [1:(NumberScenarios), fliplr(1:(NumberScenarios))];
inBetween = [min(max(SAoptCPUData2,0.01),[],2)',fliplr(max(max(SAoptCPUData2,0.01),[],2)')];
fill(x2, inBetween,[0.83, 0.13, 0.18],'linestyle', 'none')
alpha(0.15)
plot(1:(NumberScenarios), mean(max(SAoptCPUData2,0.01),2),'color',[0.8, 0.33, 0],'linewidth', 0.2,'linestyle',':')

plot(1:(NumberScenarios),min(SurrOptCPUData2,[],2), 'color',[0.8, 0.33, 0],'linestyle', 'none')
x2 = [1:(NumberScenarios), fliplr(1:(NumberScenarios))];
inBetween = [min(SurrOptCPUData2,[],2)',fliplr(max(SurrOptCPUData2,[],2)')];
fill(x2, inBetween,[0.8, 0.33, 0],'linestyle', 'none')
alpha(0.15)
plot(1:(NumberScenarios), mean(SurrOptCPUData2,2),'color',[0.8, 0.33, 0],'linewidth', 0.2,'linestyle','--')
xlabel('Scenario')
ylabel('CPU (secs)')

plot(1:(NumberScenarios),min(BayesOptCPUData2,[],2), 'color',[0, 0.55, 0.55],'linestyle', 'none')
x2 = [1:(NumberScenarios), fliplr(1:(NumberScenarios))];
inBetween = [min(BayesOptCPUData2,[],2)',fliplr(max(BayesOptCPUData2,[],2)')];
fill(x2, inBetween,[0, 0.55, 0.55],'linestyle', 'none')
alpha(0.15)
plot(1:(NumberScenarios), mean(BayesOptCPUData2,2),'color',[0, 0.55, 0.55],'linewidth', 0.2,'linestyle','--')

set(gca, 'YScale','log')
yticks([0.001 0.01 0.1 1 10 100])
ax = gca;
ax.YAxis.MinorTick = 'off';
ax.XGrid = 'on';
ax.YGrid = 'on';
ax.YMinorGrid = 'off';
box off;
