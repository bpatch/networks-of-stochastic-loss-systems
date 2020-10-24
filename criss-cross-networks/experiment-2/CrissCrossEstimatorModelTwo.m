function [Blocking1, Blocking2] = CrissCrossEstimatorModelTwo(c, Lambda, Mu, Alpha, NumberOfStations, BatchTime, DiscardedBatches, deltaTilde, epsilonTilde, MaxTime) 
    
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
    CurrentlyProcessing = zeros(1,NumberOfStations);
    floorC = floor(c);
    Cremainders = c - floorC;
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
                EntryThresholds = rand(1, NumberOfStations);
                DestinationStation = find(((CurrentlyProcessing < floorC) | ((CurrentlyProcessing < c) & (EntryThresholds <= Cremainders))),1);
                if ~isempty(DestinationStation)
                    TotalArrivals(1,1:DestinationStation) = TotalArrivals(1,1:DestinationStation) + 1;
                    TotalLosses(1,1:(DestinationStation-1)) = TotalLosses(1,1:(DestinationStation-1)) + 1;
                    CurrentlyProcessing(1,DestinationStation) = CurrentlyProcessing(1,DestinationStation) + 1;
                else
                    TotalArrivals(1,:) = TotalArrivals(1,:) + 1;
                    TotalLosses(1,:) = TotalLosses(1,:)  + 1;
                end
            else % Arrival to right-most station
                EntryThresholds = rand(1, NumberOfStations);
                DestinationStation = find(((CurrentlyProcessing < floorC) | ((CurrentlyProcessing < c) & (EntryThresholds <= Cremainders))));
                if ~isempty(DestinationStation)
                    DestinationStation = DestinationStation(end); 
                    TotalArrivals(2, DestinationStation:end) = TotalArrivals(2, DestinationStation:end) + 1;
                    TotalLosses(2,(DestinationStation+1):end) = TotalLosses(2, (DestinationStation+1):end) + 1;
                    CurrentlyProcessing(1,DestinationStation) = CurrentlyProcessing(1,DestinationStation) + 1;
                else
                    TotalArrivals(2, :) = TotalArrivals(2, :) + 1;
                    TotalLosses(2, :) = TotalLosses(2, :)  + 1;
                end
            end
        elseif NextEventType == 2 % Service
                StationOfService = find(rand*TotalServiceRate <= cumsum(ServiceRates),1);
                CurrentlyProcessing(1,StationOfService) = CurrentlyProcessing(1,StationOfService) - 1; % Job leaves station
        else % Change of arrival rate
            if CurrentEnvironment == 3
                CurrentEnvironment = 1;
            else
                CurrentEnvironment = CurrentEnvironment + 1;
            end
        end
        ServiceRates = CurrentlyProcessing.*Mu;
        TotalServiceRate = sum(ServiceRates); 
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
                EntryThresholds = rand(1, NumberOfStations);
                DestinationStation = find((CurrentlyProcessing < floorC | (CurrentlyProcessing < c & EntryThresholds <= Cremainders)), 1);
                if ~isempty(DestinationStation)
                    TotalArrivals(1,1:DestinationStation) = TotalArrivals(1,1:DestinationStation) + 1;
                    TotalLosses(1,1:(DestinationStation-1)) = TotalLosses(1,1:(DestinationStation-1)) + 1;
                    CurrentlyProcessing(1,DestinationStation) = CurrentlyProcessing(1,DestinationStation) + 1;
                else
                    TotalArrivals(1,:) = TotalArrivals(1,:) + 1;
                    TotalLosses(1,:) = TotalLosses(1,:)  + 1;
                end
            else % Arrival to right-most station
                EntryThresholds = rand(1, NumberOfStations);
                DestinationStation = find((CurrentlyProcessing < floorC | (CurrentlyProcessing < c & EntryThresholds<= Cremainders)));
                if ~isempty(DestinationStation)
                    DestinationStation = DestinationStation(end); 
                    TotalArrivals(2, DestinationStation:end) = TotalArrivals(2, DestinationStation:end) + 1;
                    TotalLosses(2,(DestinationStation+1):end) = TotalLosses(2, (DestinationStation+1):end) + 1;
                    CurrentlyProcessing(1,DestinationStation) = CurrentlyProcessing(1,DestinationStation) + 1;
                else
                    TotalArrivals(2, :) = TotalArrivals(2, :) + 1;
                    TotalLosses(2, :) = TotalLosses(2, :)  + 1;
                end
            end
        elseif NextEventType == 2 % Service
                StationOfService = find(rand*TotalServiceRate <= cumsum(ServiceRates),1);
                CurrentlyProcessing(1,StationOfService) = CurrentlyProcessing(1,StationOfService) - 1; % Job leaves station
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
%         CI2 = [CI1; 1.96.*sqrt(1.96.^2-4.*sum(TotalArrivals).*(Blocking2-1).*Blocking2)./(1.96.^2+sum(TotalArrivals))];
        CIsize = max(max(CI1));
    end
    for i = 1:NumberOfStations
%         for j = 1:2
            if c(1,i) < 1 && TotalArrivals(1,i) == 0
                Blocking1(1,i) = 1;
            elseif TotalArrivals(1,i) == 0 && Blocking1(2,i)>0
                Blocking1(1,i) = Blocking1(2,i);
            end
            
            if c(1,i) < 1 && TotalArrivals(2,i) == 0
                Blocking1(2,i) = 1;
            elseif TotalArrivals(2,i) == 0 && Blocking1(1,i)>0
                Blocking1(2,i) =  Blocking1(1,i);
            end
            
%         end
    end
    Blocking1 = max(Blocking1, epsilonTilde*(floorC>0));
end
        
