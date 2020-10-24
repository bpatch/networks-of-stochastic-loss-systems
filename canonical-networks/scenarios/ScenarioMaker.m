
ChosenScenarios = randsample(729,100);
counter = 1;
counter2 = 1;
for Web2 = [2, 3.75, 5.5]
    for Web1 = [2, 3.75, 5.5] 
        for Arr2 = [0.75 2 3.25] 
            for Arr1 = [0.75 2 3.25] 
                for App1 = [1.5 3 4.5]
                    for App2 = [1.5 3 4.5]
                        if ismember(counter, ChosenScenarios)
                            csvwrite(strcat('Scenario',num2str(counter2),'.csv'), [App1, App2, Arr1, Arr2, 1.5, 1.5, Web1, Web2]);
                            counter2 = counter2 + 1;
                        end
                        counter = counter + 1;
                    end
                end
            end
        end
    end
end
