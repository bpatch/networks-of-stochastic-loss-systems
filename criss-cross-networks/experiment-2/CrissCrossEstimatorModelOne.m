function [Blocking1, Blocking2] = CrissCrossEstimatorModelOne(C, Lambda, Mu, Alpha, NumberOfStations, BatchTime, DiscardedBatches, deltaTilde, epsilonTilde, MaxTime) 
    
    % Inputs
%     clear all; clc;
%     NumberOfStations = 5;
%     C = randi(100,1,NumberOfStations);
%     Lambda = [10,20,30];
%     Mu = rand(2,NumberOfStations);
%     alpha = 1;
%     BatchTime = 10;
%     DiscardedBatches = 5;
%     deltaTilde = 0.1;n
%     epsilonTilde = 0.0001;
    %
    TotalArrivals = zeros(2,NumberOfStations);
    TotalLosses = zeros(2,NumberOfStations);
    CurrentlyProcessing = zeros(2,NumberOfStations);
    Cremainders = C - floor(C);
    CurrentEnvironment = 1;
    NextEventTime = -log(rand)/(Lambda(1,1)+Alpha);
    if rand*(Lambda(1,1)+Alpha) < Lambda(1,1)
        NextEventType = 1;
    else
        NextEventType = 3;
    end
    ServiceRates = CurrentlyProcessing.*Mu;
    TotalServiceRate = sum(sum(ServiceRates));
    while NextEventTime <= (BatchTime*DiscardedBatches)
        if NextEventType == 1 % Arrival
            if rand <= 0.5 % Arrival to left-most station
                TotalArrivals(1,1) = TotalArrivals(1,1) + 1;
                if (CurrentlyProcessing(1,1)+CurrentlyProcessing(2,1)) < C(1,1) % Arrival accepted
                    CurrentlyProcessing(1,1) = CurrentlyProcessing(1,1) + 1;
                elseif (CurrentlyProcessing(1,1)+CurrentlyProcessing(2,1)) <= C(1,1) && rand < Cremainders(1,1) % Arrival accepted
                    CurrentlyProcessing(1,1) = CurrentlyProcessing(1,1) + 1;
                else % Arrival blocked
                    TotalLosses(1,1) = TotalLosses(1,1) + 1;
                end
            else % Arrival to right-most station
                TotalArrivals(2,NumberOfStations) = TotalArrivals(2,NumberOfStations) + 1;
                if (CurrentlyProcessing(1,NumberOfStations)+CurrentlyProcessing(2,NumberOfStations)) < C(1,NumberOfStations) % Arrival accepted
                    CurrentlyProcessing(2,NumberOfStations) = CurrentlyProcessing(2,NumberOfStations) + 1;
                elseif (CurrentlyProcessing(1,NumberOfStations)+CurrentlyProcessing(2,NumberOfStations)) <= C(1,NumberOfStations) && rand < Cremainders(1,NumberOfStations) % Arrival accepted
                    CurrentlyProcessing(2,NumberOfStations) = CurrentlyProcessing(2,NumberOfStations) + 1;
                else % Arrival blocked
                    TotalLosses(2,NumberOfStations) = TotalLosses(2,NumberOfStations) + 1;
                end
            end
        elseif NextEventType == 2 % Service
            LeftRightRate = sum(ServiceRates(1,:));
            if rand*TotalServiceRate < LeftRightRate % Left-right route service
                StationOfService = find(rand*LeftRightRate <= cumsum(ServiceRates(1,:)),1);
                CurrentlyProcessing(1,StationOfService) = CurrentlyProcessing(1,StationOfService) - 1; % Job leaves station
                if StationOfService < NumberOfStations % Job tries to move to next station
                    TotalArrivals(1,StationOfService+1) = TotalArrivals(1,StationOfService+1) + 1;
                    if (CurrentlyProcessing(1,StationOfService+1)+CurrentlyProcessing(2,StationOfService+1)) < C(1,StationOfService+1) % Arrival accepted
                        CurrentlyProcessing(1,StationOfService+1) = CurrentlyProcessing(1,StationOfService+1) + 1; 
                    elseif (CurrentlyProcessing(1,StationOfService+1)+CurrentlyProcessing(2,StationOfService+1)) <= C(1,StationOfService+1) && rand < Cremainders(1,StationOfService+1) % Arrival accepted
                        CurrentlyProcessing(1,StationOfService+1) = CurrentlyProcessing(1,StationOfService+1) + 1; 
                    else % Job blocked
                        TotalLosses(1,StationOfService+1) = TotalLosses(1,StationOfService+1) + 1;
                    end
                end
            else % Right-left route service
                RightLeftRate = sum(ServiceRates(2,:));
                StationOfService = find(rand*RightLeftRate <= cumsum(ServiceRates(2,:)),1);
                CurrentlyProcessing(2,StationOfService) = CurrentlyProcessing(2,StationOfService) - 1; % Job leaves station
                if StationOfService > 1 % Job tries to move to next station
                    TotalArrivals(2,StationOfService-1) = TotalArrivals(2,StationOfService-1) + 1;
                    if (CurrentlyProcessing(1,StationOfService-1)+CurrentlyProcessing(2,StationOfService-1)) < C(1,StationOfService-1) % Arrival accepted
                        CurrentlyProcessing(2,StationOfService-1) = CurrentlyProcessing(2,StationOfService-1) + 1; 
                    elseif (CurrentlyProcessing(1,StationOfService-1)+CurrentlyProcessing(2,StationOfService-1)) <= C(1,StationOfService-1) && rand < Cremainders(1,StationOfService-1) % Arrival accepted
                        CurrentlyProcessing(2,StationOfService-1) = CurrentlyProcessing(2,StationOfService-1) + 1; 
                    else % Job blocked
                        TotalLosses(2,StationOfService-1) = TotalLosses(2,StationOfService-1) + 1;
                    end
                end
            end
        else % Change of arrival rate
            if CurrentEnvironment == 3
                CurrentEnvironment = 1;
            else
                CurrentEnvironment = CurrentEnvironment + 1;
            end
        end
        ServiceRates = CurrentlyProcessing.*Mu;
        TotalServiceRate = sum(sum(ServiceRates)); 
        r = rand*(Lambda(1,CurrentEnvironment)+TotalServiceRate+Alpha);
        if r <= Lambda(1,CurrentEnvironment)
            NextEventType = 1;
        elseif r <= (Lambda(1,CurrentEnvironment)+TotalServiceRate)
            NextEventType = 2;
        else 
            NextEventType = 3;
        end
        NextEventTime = NextEventTime - log(rand)/(Lambda(1,CurrentEnvironment)+TotalServiceRate+Alpha);
    end
    TotalArrivals = zeros(2,NumberOfStations);
    TotalLosses = zeros(2,NumberOfStations);
    CIsize = 1;
    while CIsize > deltaTilde && NextEventTime < MaxTime
        ThisBatchMaxTime = NextEventTime + BatchTime;
        while NextEventTime <= ThisBatchMaxTime
            if NextEventType == 1 % Arrival
                if rand <= 0.5 % Arrival to left-most station
                    TotalArrivals(1,1) = TotalArrivals(1,1) + 1;
                    if (CurrentlyProcessing(1,1)+CurrentlyProcessing(2,1)) < C(1,1) % Arrival accepted
                        CurrentlyProcessing(1,1) = CurrentlyProcessing(1,1) + 1;
                    elseif (CurrentlyProcessing(1,1)+CurrentlyProcessing(2,1)) <= C(1,1) && rand < Cremainders(1,1) % Arrival accepted
                        CurrentlyProcessing(1,1) = CurrentlyProcessing(1,1) + 1;
                    else % Arrival blocked
                        TotalLosses(1,1) = TotalLosses(1,1) + 1;
                    end
                else % Arrival to right-most station
                    TotalArrivals(2,NumberOfStations) = TotalArrivals(2,NumberOfStations) + 1;
                    if (CurrentlyProcessing(1,NumberOfStations)+CurrentlyProcessing(2,NumberOfStations)) < C(1,NumberOfStations) % Arrival accepted
                        CurrentlyProcessing(2,NumberOfStations) = CurrentlyProcessing(2,NumberOfStations) + 1;
                    elseif (CurrentlyProcessing(1,NumberOfStations)+CurrentlyProcessing(2,NumberOfStations)) <= C(1,NumberOfStations) && rand < Cremainders(1,NumberOfStations) % Arrival accepted
                        CurrentlyProcessing(2,NumberOfStations) = CurrentlyProcessing(2,NumberOfStations) + 1;
                    else % Arrival blocked
                        TotalLosses(2,NumberOfStations) = TotalLosses(2,NumberOfStations) + 1;
                    end
                end
            elseif NextEventType == 2 % Service
                LeftRightRate = sum(ServiceRates(1,:));
                if rand*TotalServiceRate < LeftRightRate % Left-right route service
                    StationOfService = find(rand*LeftRightRate <= cumsum(ServiceRates(1,:)),1);
                    CurrentlyProcessing(1,StationOfService) = CurrentlyProcessing(1,StationOfService) - 1; % Job leaves station
                    if StationOfService < NumberOfStations % Job tries to move to next station
                        TotalArrivals(1,StationOfService+1) = TotalArrivals(1,StationOfService+1) + 1;
                        if (CurrentlyProcessing(1,StationOfService+1)+CurrentlyProcessing(2,StationOfService+1)) < C(1,StationOfService+1) % Arrival accepted
                            CurrentlyProcessing(1,StationOfService+1) = CurrentlyProcessing(1,StationOfService+1) + 1; 
                        elseif (CurrentlyProcessing(1,StationOfService+1)+CurrentlyProcessing(2,StationOfService+1)) <= C(1,StationOfService+1) && rand < Cremainders(1,StationOfService+1) % Arrival accepted
                            CurrentlyProcessing(1,StationOfService+1) = CurrentlyProcessing(1,StationOfService+1) + 1; 
                        else % Job blocked
                            TotalLosses(1,StationOfService+1) = TotalLosses(1,StationOfService+1) + 1;
                        end
                    end
                else % Right-left route service
                    RightLeftRate = sum(ServiceRates(2,:));
                    StationOfService = find(rand*RightLeftRate <= cumsum(ServiceRates(2,:)),1);
                    CurrentlyProcessing(2,StationOfService) = CurrentlyProcessing(2,StationOfService) - 1; % Job leaves station
                    if StationOfService > 1 % Job tries to move to next station
                        TotalArrivals(2,StationOfService-1) = TotalArrivals(2,StationOfService-1) + 1;
                        if (CurrentlyProcessing(1,StationOfService-1)+CurrentlyProcessing(2,StationOfService-1)) < C(1,StationOfService-1) % Arrival accepted
                            CurrentlyProcessing(2,StationOfService-1) = CurrentlyProcessing(2,StationOfService-1) + 1; 
                        elseif (CurrentlyProcessing(1,StationOfService-1)+CurrentlyProcessing(2,StationOfService-1)) <= C(1,StationOfService-1) && rand < Cremainders(1,StationOfService-1) % Arrival accepted
                            CurrentlyProcessing(2,StationOfService-1) = CurrentlyProcessing(2,StationOfService-1) + 1; 
                        else % Job blocked
                            TotalLosses(2,StationOfService-1) = TotalLosses(2,StationOfService-1) + 1;
                        end
                    end
                end
            else % Change of arrival rate
                if CurrentEnvironment == 3
                    CurrentEnvironment = 1;
                else
                    CurrentEnvironment = CurrentEnvironment + 1;
                end
            end
            ServiceRates = CurrentlyProcessing.*Mu;
            TotalServiceRate = sum(sum(ServiceRates)); 
            r = rand*(Lambda(1,CurrentEnvironment)+TotalServiceRate+Alpha);
            if r <= Lambda(1,CurrentEnvironment)
                NextEventType = 1;
            elseif r <= (Lambda(1,CurrentEnvironment)+TotalServiceRate)
                NextEventType = 2;
            else 
                NextEventType = 3;
            end
            NextEventTime = NextEventTime - log(rand)/(Lambda(1,CurrentEnvironment)+TotalServiceRate+Alpha);
        end
        Blocking1 = TotalLosses./TotalArrivals;
        Blocking2 = sum(TotalLosses)./sum(TotalArrivals);
        CI1 = 1.96.*sqrt(1.96.^2-4.*TotalArrivals.*(Blocking1-1).*Blocking1)./(1.96.^2+TotalArrivals);
        CI2 = [CI1; 1.96.*sqrt(1.96.^2-4.*sum(TotalArrivals).*(Blocking2-1).*Blocking2)./(1.96.^2+sum(TotalArrivals))];
        CIsize = max(max(CI2));
    end
    for i = 1:NumberOfStations
        if Blocking2(1,i) == 0
            Blocking2(1,i) = epsilonTilde;
        end
        for j = 1:2
            if Blocking1(j,i) == 0
                Blocking1(j,i) = epsilonTilde;
            end
        end
    end
    Blocking1 = max(Blocking1, epsilonTilde);
    Blocking2 = max(Blocking2, epsilonTilde);
end
        
