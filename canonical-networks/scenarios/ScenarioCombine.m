clear all; clc;
NumberScenarios = 100;
data = NaN(100,8);

for ScenarioNumber = 1:NumberScenarios
     data(ScenarioNumber, :) = load(strcat('Scenario',num2str(ScenarioNumber),'.csv'));
end
csvwrite('AllScenarios.csv',data)