clear all; close all; clc;
load('ModelOneSimulationData.mat')

SurrOptCPUData = [zeros(NumberSamples,1), SurrOptCPUData];
SurrOptObjData = [zeros(NumberSamples,1), SurrOptObjData];

BayesOptCPUData = [zeros(NumberSamples,1), BayesOptCPUData];
BayesOptObjData = [zeros(NumberSamples,1), BayesOptObjData];

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
xlabel('Iteration')
ylabel('Obj. Val.')
title('Stoch. App.') 
%axis([IntervalK MaxK 15 RealMax+0.5])
box off

figure()
plot(0:(SuggIterations),max(SurrOptObjData',[],2), 'color', [0, 0.4470, 0.7410],'linestyle', 'none')
hold on
plot(0:(SuggIterations),min(SurrOptObjData',[],2), 'color',[0, 0.4470, 0.7410],'linestyle', 'none')
x2 = [0:(SuggIterations), fliplr(0:(SuggIterations))];
inBetween = [min(SurrOptObjData',[],2)',fliplr(max(SurrOptObjData',[],2)')];
fill(x2, inBetween,[0, 0.4470, 0.7410],'linestyle', 'none')
alpha(0.15)
hold on
plot(0:(SuggIterations), SurrOptObjData','color',[0.75 0.75 0.75],'linewidth', 0.2)
plot(0:(SuggIterations), mean(SurrOptObjData', 2),'color',[0, 0.4470, 0.7410],'linewidth', 0.5)
xlabel('Iteration')
ylabel('Obj. Val.')
title('Surr. Opt.') 
%axis([IntervalK MaxK 15 RealMax+0.5])
box off

figure()
plot(0:(SuggIterations),max(BayesOptObjData',[],2), 'color', [0, 0.4470, 0.7410],'linestyle', 'none')
hold on
plot(0:(SuggIterations),min(BayesOptObjData',[],2), 'color',[0, 0.4470, 0.7410],'linestyle', 'none')
x2 = [0:(SuggIterations), fliplr(0:(SuggIterations))];
inBetween = [min(BayesOptObjData',[],2)',fliplr(max(BayesOptObjData',[],2)')];
fill(x2, inBetween,[0, 0.4470, 0.7410],'linestyle', 'none')
alpha(0.15)
hold on
plot(0:(SuggIterations), BayesOptObjData','color',[0.75 0.75 0.75],'linewidth', 0.2)
plot(0:(SuggIterations), mean(BayesOptObjData', 2),'color',[0, 0.4470, 0.7410],'linewidth', 0.5)
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


% save('ModelOne.mat')
% bestSurrSol = [30,30,28,31,33,29,31,31,32,31]; % This was manually copied and pasted (could be automated)
% 
% [~, I ] = max(FFoptObjData(:,end));
% 
% bestFFSol = FFoptSolData(end,:,I)';
% 
% [~, I ] = max(SAoptObjData(:,end));
% 
% bestSASol = SAoptSolData(end,:,I)';
% 
% figure()
% plot(1:10, bestSurrSol,'color',[0.8, 0.33, 0],'linewidth', 0.2,'linestyle','--')
% hold on
% plot(1:10, bestFFSol,'color',[0, 0.5, 0],'linewidth', 0.2)
% plot(1:10, bestSASol,'color',[0.83, 0.13, 0.18],'linewidth', 0.2,'linestyle',':')



