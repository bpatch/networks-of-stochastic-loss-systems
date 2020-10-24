function BlockingEstimate = TandemEstimatorModelOne(c, lambda, Mu, BatchTime, DiscardedBatches, deltaTilde, epsilonTilde)

    NumberStations = length(Mu);
    Cremainders = c-floor(c);
    cFloor = floor(c); %% 
    EventList = [-log(rand)/lambda,1,1;
                 inf,4,4];
    CurrentlyProcessing = zeros(NumberStations,1);

    NextEventTime = EventList(1,1);
    while NextEventTime <= (BatchTime*DiscardedBatches)
        EventType = EventList(1,2);
        if EventType == 1  % Arrival to first station
            %TotalArrivals(1,1) = TotalArrivals(1,1) + 1;
            if CurrentlyProcessing(1,1) < cFloor(1,1)
                CurrentlyProcessing(1,1) = CurrentlyProcessing(1,1) + 1;
                NewServiceTime = NextEventTime - log(rand)/Mu(1,1);
                NewEventPosition = find(EventList(:,1) >= NewServiceTime,1);
                EventList = [EventList(2:(NewEventPosition-1),:);
                            NewServiceTime,2,1;
                            EventList(NewEventPosition:end,:)];
                NewArrivalTime = NextEventTime-log(rand)/lambda;
                NewEventPosition = find(EventList(:,1) >= NewArrivalTime,1);
                EventList = [EventList(1:(NewEventPosition-1),:);
                            NewArrivalTime,1,1;
                            EventList(NewEventPosition:end,:)];   
            elseif CurrentlyProcessing(1,1) <= cFloor(1,1) && rand < Cremainders(1,1)
                CurrentlyProcessing(1,1) = CurrentlyProcessing(1,1) + 1;
                NewServiceTime = NextEventTime -log(rand)/Mu(1,1);
                NewEventPosition = find(EventList(:,1) >= NewServiceTime,1);
                EventList = [EventList(2:(NewEventPosition-1),:);
                            NewServiceTime,2,1;
                            EventList(NewEventPosition:end,:)];
                NewArrivalTime = NextEventTime-log(rand)/lambda;
                NewEventPosition = find(EventList(:,1) >= NewArrivalTime,1);
                EventList = [EventList(1:(NewEventPosition-1),:);
                            NewArrivalTime,1,1;
                            EventList(NewEventPosition:end,:)];   
            else
                %TotalLosses(1,1) = TotalLosses(1,1) + 1;
                NewArrivalTime = NextEventTime-log(rand)/lambda;
                NewEventPosition = find(EventList(:,1) >= NewArrivalTime,1);
                EventList = [EventList(2:(NewEventPosition-1),:);
                            NewArrivalTime,1,1;
                            EventList(NewEventPosition:end,:)];      
            end
        elseif EventType == 2 % Service
            StationOfService = EventList(1,3);
            for j = 1:(NumberStations-1)
                if StationOfService == j
                    %TotalArrivals(j+1,1) = TotalArrivals(j+1,1) + 1;
                    if CurrentlyProcessing(j+1,1) < cFloor(j+1,1)
                        CurrentlyProcessing(j,1) = CurrentlyProcessing(j,1) - 1;
                        CurrentlyProcessing(j+1,1) = CurrentlyProcessing(j+1,1) + 1;
                        NewServiceTime = NextEventTime -log(rand)/Mu(j+1,1);
                        NewEventPosition = find(EventList(:,1) >= NewServiceTime,1);
                        EventList = [EventList(2:(NewEventPosition-1),:);
                                    NewServiceTime,2,j+1;
                                    EventList(NewEventPosition:end,:)];
                    elseif CurrentlyProcessing(j+1,1) <= cFloor(j+1,1) && rand < Cremainders(j+1,1)
                        CurrentlyProcessing(j,1) = CurrentlyProcessing(j,1) - 1;
                        CurrentlyProcessing(j+1,1) = CurrentlyProcessing(j+1,1) + 1;
                        NewServiceTime = NextEventTime -log(rand)/Mu(j+1,1);
                        NewEventPosition = find(EventList(:,1) >= NewServiceTime,1);
                        EventList = [EventList(2:(NewEventPosition-1),:);
                                    NewServiceTime,2,j+1;
                                    EventList(NewEventPosition:end,:)];
                    else
                        CurrentlyProcessing(j,1) = CurrentlyProcessing(j,1) - 1;
                        EventList(1,:) = [];
                        %TotalLosses(j+1,1) = TotalLosses(j+1,1) + 1;
                    end
                end
            end
            if StationOfService == NumberStations
                CurrentlyProcessing(NumberStations,1) = CurrentlyProcessing(NumberStations,1) - 1;
                EventList(1,:) = [];
            end
        end
        NextEventTime = EventList(1,1);
    end
    TotalArrivals = zeros(NumberStations,1);
    TotalLosses = zeros(NumberStations,1);
    CIsize = 1;
    while CIsize > deltaTilde
        ThisBatchMaxTime = NextEventTime + BatchTime;
        while NextEventTime <= ThisBatchMaxTime
            EventType = EventList(1,2);
            if EventType == 1  % Arrival to first station
                TotalArrivals(1,1) = TotalArrivals(1,1) + 1;
                if CurrentlyProcessing(1,1) < cFloor(1,1)
                    CurrentlyProcessing(1,1) = CurrentlyProcessing(1,1) + 1;
                    NewServiceTime = NextEventTime -log(rand)/Mu(1,1);
                    NewEventPosition = find(EventList(:,1) >= NewServiceTime,1);
                    EventList = [EventList(2:(NewEventPosition-1),:);
                                NewServiceTime,2,1;
                                EventList(NewEventPosition:end,:)];
                    NewArrivalTime = NextEventTime-log(rand)/lambda;
                    NewEventPosition = find(EventList(:,1) >= NewArrivalTime,1);
                    EventList = [EventList(1:(NewEventPosition-1),:);
                                NewArrivalTime,1,1;
                                EventList(NewEventPosition:end,:)];   
                elseif CurrentlyProcessing(1,1) <= cFloor(1,1) && rand < Cremainders(1,1)
                    CurrentlyProcessing(1,1) = CurrentlyProcessing(1,1) + 1;
                    NewServiceTime = NextEventTime -log(rand)/Mu(1,1);
                    NewEventPosition = find(EventList(:,1) >= NewServiceTime,1);
                    EventList = [EventList(2:(NewEventPosition-1),:);
                                NewServiceTime,2,1;
                                EventList(NewEventPosition:end,:)];                        
                    NewArrivalTime = NextEventTime-log(rand)/lambda;
                    NewEventPosition = find(EventList(:,1) >= NewArrivalTime,1);
                    EventList = [EventList(1:(NewEventPosition-1),:);
                                NewArrivalTime,1,1;
                                EventList(NewEventPosition:end,:)];   
                else
                    TotalLosses(1,1) = TotalLosses(1,1) + 1;
                    NewArrivalTime = NextEventTime-log(rand)/lambda;
                    NewEventPosition = find(EventList(:,1) >= NewArrivalTime,1);
                    EventList = [EventList(2:(NewEventPosition-1),:);
                                NewArrivalTime,1,1;
                                EventList(NewEventPosition:end,:)];      
                end
            elseif EventType == 2 % Service
                StationOfService = EventList(1,3);
                for j = 1:(NumberStations-1)
                    if StationOfService == j
                        TotalArrivals(j+1,1) = TotalArrivals(j+1,1) + 1;
                        if CurrentlyProcessing(j+1,1) < cFloor(j+1,1)
                            CurrentlyProcessing(j,1) = CurrentlyProcessing(j,1) - 1;
                            CurrentlyProcessing(j+1,1) = CurrentlyProcessing(j+1,1) + 1;
                            NewServiceTime = NextEventTime -log(rand)/Mu(j+1,1);
                            NewEventPosition = find(EventList(:,1) >= NewServiceTime,1);
                            EventList = [EventList(2:(NewEventPosition-1),:);
                                        NewServiceTime,2,j+1;
                                        EventList(NewEventPosition:end,:)];
                        elseif CurrentlyProcessing(j+1,1) <= cFloor(j+1,1) && rand < Cremainders(j+1,1)
                            CurrentlyProcessing(j,1) = CurrentlyProcessing(j,1) - 1;
                            CurrentlyProcessing(j+1,1) = CurrentlyProcessing(j+1,1) + 1;
                            NewServiceTime = NextEventTime -log(rand)/Mu(j+1,1);
                            NewEventPosition = find(EventList(:,1) >= NewServiceTime,1);
                            EventList = [EventList(2:(NewEventPosition-1),:);
                                        NewServiceTime,2,j+1;
                                        EventList(NewEventPosition:end,:)];
                        else
                            CurrentlyProcessing(j,1) = CurrentlyProcessing(j,1) - 1;
                            EventList(1,:) = [];
                            TotalLosses(j+1,1) = TotalLosses(j+1,1) + 1;
                        end
                    end
                end
                if StationOfService == NumberStations
                    CurrentlyProcessing(NumberStations,1) = CurrentlyProcessing(NumberStations,1) - 1;
                    EventList(1,:) = [];
                end
            end
            NextEventTime = EventList(1,1);
        end
        Blocking1 = TotalLosses./TotalArrivals;
        CIsize = max(max(1.96.*sqrt(1.96.^2-4.*TotalArrivals.*(Blocking1-1).*Blocking1)./((1.96.^2+TotalArrivals)))); % See p.175 of Kroese&Chan
    end
    for i = 1:2
        if c(i,1) <= 0.5 && TotalArrivals(i,1) == 0
            Blocking1(i,1) = 1;
        end
    end
    BlockingEstimate = max(epsilonTilde*(cFloor>0), Blocking1);
end