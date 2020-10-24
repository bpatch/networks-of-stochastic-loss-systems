clear all; close all; clc;
% load('ModelOneSimulationData.mat')
% 
% SurrOptObjData = NaN(NumberScenarios, NumberSamples);
% 
% BayesOptObjData = NaN(NumberScenarios, NumberSamples);
% 
% for m = 1:NumberScenarios
%     for n = 1:NumberSamples
%         SAoptObjData(m, n) = ExactEvaluationModelOne(SAoptSolData(m,1,n),SAoptSolData(m,2,n),weightValues(m,1), muValues(m,:)', lambdaValues(m,1), capacityCostValues(m,:)');
%     end
% end
% 
% for m = 1:NumberScenarios
%     for n = 1:NumberSamples
%         SurrOptObjData(m, n) = ExactEvaluationModelOne(SurrOptSolData(m,1,n),SurrOptSolData(m,2,n),weightValues(m,1), muValues(m,:)', lambdaValues(m,1), capacityCostValues(m,:)');
%     end
% end
% 
% for m = 1:NumberScenarios
%     for n = 1:NumberSamples
%         InitialObjData(m, n) = ExactEvaluationModelOne(InitialSolData(m,1,n),InitialSolData(m,2,n),weightValues(m,1), muValues(m,:)', lambdaValues(m,1), capacityCostValues(m,:)');
%     end
% end
% 
% for m = 1:NumberScenarios
%     for n = 1:NumberSamples
%         if FFoptSolData(m,1,n) > 100 || FFoptSolData(m,2,n)>100
%             FFoptObjData(m, n) = ObjectiveFunctionTandemModelOne(weightValues(m,1),capacityCostValues(m,:)',10,FFoptSolData(m,:,n)',lambdaValues(m,1),muValues(m,:)',10*BatchTime, DiscardedBatches, epsilonTilde, 0.01);
%         else
%             FFoptObjData(m, n) = ExactEvaluationModelOne(FFoptSolData(m,1,n),FFoptSolData(m,2,n),weightValues(m,1), muValues(m,:)', lambdaValues(m,1), capacityCostValues(m,:)');
%         end
%     end
% end
% 
% for m = 1:NumberScenarios
%     for n = 1:NumberSamples
%         BayesOptObjData(m, n) = ExactEvaluationModelOne(BayesOptSolData(m,1,n),BayesOptSolData(m,2,n),weightValues(m,1), muValues(m,:)', lambdaValues(m,1), capacityCostValues(m,:)');
%     end
% end
% 
% save('ModelOnePlottingData.mat')

load('ModelOnePlottingData.mat')

figure()
plot(1:(NumberScenarios), MAMoptData,'color','k','linewidth', 0.5,'linestyle','-.')
hold on
plot(1:(NumberScenarios),max(FFoptObjData,[],2), 'color', [0, 0.5, 0],'linestyle', 'none')
plot(1:(NumberScenarios),min(FFoptObjData,[],2), 'color',[0, 0.5, 0],'linestyle', 'none')
x2 = [1:(NumberScenarios), fliplr(1:(NumberScenarios))];
inBetween = [min(FFoptObjData,[],2)',fliplr(max(FFoptObjData,[],2)')];
fill(x2, inBetween,[0, 0.5, 0],'linestyle', 'none')
alpha(0.15)
plot(1:(NumberScenarios), mean(FFoptObjData,2),'color',[0, 0.5, 0],'linewidth', 0.2)

plot(1:(NumberScenarios),min(SAoptObjData,[],2), 'color',[0.83, 0.13, 0.18],'linestyle', 'none')
x2 = [1:(NumberScenarios), fliplr(1:(NumberScenarios))];
inBetween = [min(SAoptObjData,[],2)',fliplr(max(SAoptObjData,[],2)')];
fill(x2, inBetween,[0.83, 0.13, 0.18],'linestyle', 'none')
alpha(0.15)
plot(1:(NumberScenarios), mean(SAoptObjData,2),'color',[0.8, 0.33, 0],'linewidth', 0.2,'linestyle',':')

plot(1:(NumberScenarios),min(SurrOptObjData,[],2), 'color',[0.8, 0.33, 0],'linestyle', 'none')
x2 = [1:(NumberScenarios), fliplr(1:(NumberScenarios))];
inBetween = [min(SurrOptObjData,[],2)',fliplr(max(SurrOptObjData,[],2)')];
fill(x2, inBetween,[0.8, 0.33, 0],'linestyle', 'none')
alpha(0.15)
plot(1:(NumberScenarios), mean(SurrOptObjData,2),'color',[0.8, 0.33, 0],'linewidth', 0.2,'linestyle','--')


plot(1:(NumberScenarios),min(BayesOptObjData,[],2), 'color',[0, 0.55, 0.55],'linestyle', 'none')
x2 = [1:(NumberScenarios), fliplr(1:(NumberScenarios))];
inBetween = [min(BayesOptObjData,[],2)',fliplr(max(BayesOptObjData,[],2)')];
fill(x2, inBetween,[0, 0.55, 0.55],'linestyle', 'none')
alpha(0.15)
plot(1:(NumberScenarios), mean(BayesOptObjData,2),'color',[0, 0.55, 0.55],'linewidth', 0.2,'linestyle','-.')

xlabel('Scenario')
ylabel('Obj. Val.')

plot(1:(NumberScenarios),min(InitialObjData,[],2), 'color',[0, 0.28, 0.73],'linestyle', 'none')
x2 = [1:(NumberScenarios), fliplr(1:(NumberScenarios))];
inBetween = [min(InitialObjData,[],2)',fliplr(max(InitialObjData,[],2)')];
fill(x2, inBetween,[0, 0.28, 0.73],'linestyle', 'none')
alpha(0.15)
plot(1:(NumberScenarios), mean(InitialObjData,2),'color',[0, 0.28, 0.73],'linewidth', 0.2,'linestyle','--')

ax = gca;
ax.YAxis.MinorTick = 'off';
ax.XGrid = 'on';
ax.YGrid = 'on';
box off; 

figure()
plot(1:(NumberScenarios),max(FFoptCPUData,[],2), 'color', [0, 0.5, 0],'linestyle', 'none')
hold on
plot(1:(NumberScenarios),min(FFoptCPUData,[],2), 'color',[0, 0.5, 0],'linestyle', 'none')
x2 = [1:(NumberScenarios), fliplr(1:(NumberScenarios))];
inBetween = max([min(FFoptCPUData,[],2)',fliplr(max(FFoptCPUData,[],2)')],0.001);
fill(x2, inBetween,[0, 0.5, 0],'linestyle', 'none')
alpha(0.15)
plot(1:(NumberScenarios), mean(FFoptCPUData,2),'color',[0, 0.5, 0],'linewidth', 0.2)

plot(1:(NumberScenarios),min(max(SAoptCPUData,0.01),[],2), 'color',[0.83, 0.13, 0.18],'linestyle', 'none')
x2 = [1:(NumberScenarios), fliplr(1:(NumberScenarios))];
inBetween = [min(max(SAoptCPUData,0.01),[],2)',fliplr(max(max(SAoptCPUData,0.01),[],2)')];
fill(x2, inBetween,[0.83, 0.13, 0.18],'linestyle', 'none')
alpha(0.15)
plot(1:(NumberScenarios), mean(max(SAoptCPUData,0.01),2),'color',[0.8, 0.33, 0],'linewidth', 0.2,'linestyle',':')

plot(1:(NumberScenarios),min(SurrOptCPUData,[],2), 'color',[0.8, 0.33, 0],'linestyle', 'none')
x2 = [1:(NumberScenarios), fliplr(1:(NumberScenarios))];
inBetween = [min(SurrOptCPUData,[],2)',fliplr(max(SurrOptCPUData,[],2)')];
fill(x2, inBetween,[0.8, 0.33, 0],'linestyle', 'none')
alpha(0.15)
plot(1:(NumberScenarios), mean(SurrOptCPUData,2),'color',[0.8, 0.33, 0],'linewidth', 0.2,'linestyle','--')


plot(1:(NumberScenarios),min(BayesOptCPUData,[],2), 'color',[0, 0.55, 0.55],'linestyle', 'none')
x2 = [1:(NumberScenarios), fliplr(1:(NumberScenarios))];
inBetween = [min(SurrOptCPUData,[],2)',fliplr(max(BayesOptCPUData,[],2)')];
fill(x2, inBetween,[0, 0.55, 0.55],'linestyle', 'none')
alpha(0.15)
plot(1:(NumberScenarios), mean(BayesOptCPUData,2),'color',[0, 0.55, 0.55],'linewidth', 0.2,'linestyle','-.')


xlabel('Scenario')
ylabel('CPU (secs)')
set(gca, 'YScale','log')
yticks([0.001 0.01 0.1 1 10 100])
ax = gca;
ax.YAxis.MinorTick = 'off';
ax.XGrid = 'on';
ax.YGrid = 'on';
ax.YMinorGrid = 'off';
box off;

SurrOptSolDataStation1 = [SurrOptSolData(:,1,1),SurrOptSolData(:,1,2), SurrOptSolData(:,1,3)];
SurrOptSolDataStation2 = [SurrOptSolData(:,2,1),SurrOptSolData(:,2,2), SurrOptSolData(:,2,3)];

BayesOptSolDataStation1 = [BayesOptSolData(:,1,1),BayesOptSolData(:,1,2), BayesOptSolData(:,1,3)];
BayesOptSolDataStation2 = [BayesOptSolData(:,2,1),BayesOptSolData(:,2,2), BayesOptSolData(:,2,3)];


SAoptSolDataStation1 = [SAoptSolData(:,1,1),SAoptSolData(:,1,2), SAoptSolData(:,1,3)];
SAoptSolDataStation2 = [SAoptSolData(:,2,1),SAoptSolData(:,2,2), SAoptSolData(:,2,3)];

FFoptSolDataStation1 = [FFoptSolData(:,1,1),FFoptSolData(:,1,2), FFoptSolData(:,1,3)];
FFoptSolDataStation2 = [FFoptSolData(:,2,1),FFoptSolData(:,2,2), FFoptSolData(:,2,3)];

InitialSolDataStation1 = [InitialSolData(:,1,1),InitialSolData(:,1,2), InitialSolData(:,1,3)];
InitialSolDataStation2 = [InitialSolData(:,2,1),InitialSolData(:,2,2), InitialSolData(:,2,3)];

figure()
plot(1:(NumberScenarios),max(FFoptSolDataStation1,[],2), 'color', [0, 0.5, 0],'linestyle', 'none')
hold on

plot(1:(NumberScenarios),min(InitialSolDataStation1,[],2), 'color',[0, 0.28, 0.73],'linestyle', 'none')
x2 = [1:(NumberScenarios), fliplr(1:(NumberScenarios))];
inBetween = [min(InitialSolDataStation1,[],2)',fliplr(max(InitialSolDataStation1,[],2)')];
fill(x2, inBetween,[0, 0.28, 0.73],'linestyle', 'none')
alpha(0.15)
plot(1:(NumberScenarios), mean(InitialSolDataStation1,2),'color',[0, 0.28, 0.73],'linewidth', 0.2,'linestyle','--')


plot(1:(NumberScenarios),min(FFoptSolDataStation1,[],2), 'color',[0, 0.5, 0],'linestyle', 'none')
x2 = [1:(NumberScenarios), fliplr(1:(NumberScenarios))];
inBetween = [min(FFoptSolDataStation1,[],2)',fliplr(max(FFoptSolDataStation1,[],2)')];
fill(x2, inBetween,[0, 0.5, 0],'linestyle', 'none')
alpha(0.15)
plot(1:(NumberScenarios), mean(FFoptSolDataStation1,2),'color',[0, 0.5, 0],'linewidth', 0.2)

plot(1:(NumberScenarios),min(SAoptSolDataStation1,[],2), 'color',[0.83, 0.13, 0.18],'linestyle', 'none')
x2 = [1:(NumberScenarios), fliplr(1:(NumberScenarios))];
inBetween = [min(SAoptSolDataStation1,[],2)',fliplr(max(SAoptSolDataStation1,[],2)')];
fill(x2, inBetween,[0.83, 0.13, 0.18],'linestyle', 'none')
alpha(0.15)
plot(1:(NumberScenarios), mean(SAoptSolDataStation1,2),'color',[0.8, 0.33, 0],'linewidth', 0.2,'linestyle',':')

plot(1:(NumberScenarios),min(SurrOptSolDataStation1,[],2), 'color',[0.8, 0.33, 0],'linestyle', 'none')
x2 = [1:(NumberScenarios), fliplr(1:(NumberScenarios))];
inBetween = [min(SurrOptSolDataStation1,[],2)',fliplr(max(SurrOptSolDataStation1,[],2)')];
fill(x2, inBetween,[0.8, 0.33, 0],'linestyle', 'none')
alpha(0.15)
plot(1:(NumberScenarios), mean(SurrOptSolDataStation1,2),'color',[0.8, 0.33, 0],'linewidth', 0.2,'linestyle','--')

plot(1:(NumberScenarios),min(BayesOptSolDataStation1,[],2), 'color',[0, 0.55, 0.55],'linestyle', 'none')
x2 = [1:(NumberScenarios), fliplr(1:(NumberScenarios))];
inBetween = [min(BayesOptSolDataStation1,[],2)',fliplr(max(BayesOptSolDataStation1,[],2)')];
fill(x2, inBetween,[0, 0.55, 0.55],'linestyle', 'none')
alpha(0.15)
plot(1:(NumberScenarios), mean(BayesOptSolDataStation1,2),'color',[0, 0.55, 0.55],'linewidth', 0.2,'linestyle','-.')




xlabel('Scenario')
ylabel('c_1')

plot(1:(NumberScenarios), MaMsolData(:,1),'color','k','linewidth', 0.5,'linestyle','-.')
ax = gca;
ax.YAxis.MinorTick = 'off';
ax.XGrid = 'on';
ax.YGrid = 'on';
box off; 


figure()
plot(1:(NumberScenarios),max(FFoptSolDataStation2,[],2), 'color', [0, 0.5, 0],'linestyle', 'none')
hold on

plot(1:(NumberScenarios),min(InitialSolDataStation2,[],2), 'color',[0, 0.28, 0.73],'linestyle', 'none')
x2 = [1:(NumberScenarios), fliplr(1:(NumberScenarios))];
inBetween = [min(InitialSolDataStation2,[],2)',fliplr(max(InitialSolDataStation2,[],2)')];
fill(x2, inBetween,[0, 0.28, 0.73],'linestyle', 'none')
alpha(0.15)
plot(1:(NumberScenarios), mean(InitialSolDataStation2,2),'color',[0, 0.28, 0.73],'linewidth', 0.2,'linestyle','--')


plot(1:(NumberScenarios),min(FFoptSolDataStation2,[],2), 'color',[0, 0.5, 0],'linestyle', 'none')
x2 = [1:(NumberScenarios), fliplr(1:(NumberScenarios))];
inBetween = [min(FFoptSolDataStation2,[],2)',fliplr(max(FFoptSolDataStation2,[],2)')];
fill(x2, inBetween,[0, 0.5, 0],'linestyle', 'none')
alpha(0.15)
plot(1:(NumberScenarios), mean(FFoptSolDataStation2,2),'color',[0, 0.5, 0],'linewidth', 0.2)

plot(1:(NumberScenarios),min(SAoptSolDataStation1,[],2), 'color',[0.83, 0.13, 0.18],'linestyle', 'none')
x2 = [1:(NumberScenarios), fliplr(1:(NumberScenarios))];
inBetween = [min(SAoptSolDataStation2,[],2)',fliplr(max(SAoptSolDataStation2,[],2)')];
fill(x2, inBetween,[0.83, 0.13, 0.18],'linestyle', 'none')
alpha(0.15)
plot(1:(NumberScenarios), mean(SAoptSolDataStation2,2),'color',[0.8, 0.33, 0],'linewidth', 0.2,'linestyle',':')

plot(1:(NumberScenarios),min(SurrOptSolDataStation2,[],2), 'color',[0.8, 0.33, 0],'linestyle', 'none')
x2 = [1:(NumberScenarios), fliplr(1:(NumberScenarios))];
inBetween = [min(SurrOptSolDataStation2,[],2)',fliplr(max(SurrOptSolDataStation2,[],2)')];
fill(x2, inBetween,[0.8, 0.33, 0],'linestyle', 'none')
alpha(0.15)
plot(1:(NumberScenarios), mean(SurrOptSolDataStation2,2),'color',[0.8, 0.33, 0],'linewidth', 0.2,'linestyle','--')


plot(1:(NumberScenarios),min(BayesOptSolDataStation2,[],2), 'color',[0, 0.55, 0.55],'linestyle', 'none')
x2 = [1:(NumberScenarios), fliplr(1:(NumberScenarios))];
inBetween = [min(BayesOptSolDataStation2,[],2)',fliplr(max(BayesOptSolDataStation2,[],2)')];
fill(x2, inBetween,[0, 0.55, 0.55],'linestyle', 'none')
alpha(0.15)
plot(1:(NumberScenarios), mean(BayesOptSolDataStation2,2),'color',[0, 0.55, 0.55],'linewidth', 0.2,'linestyle','--')
xlabel('Scenario')
ylabel('c_2')


plot(1:(NumberScenarios), MaMsolData(:,2),'color','k','linewidth', 0.5,'linestyle','-.')
ax = gca;
ax.YAxis.MinorTick = 'off';
ax.XGrid = 'on';
ax.YGrid = 'on';
box off; 


