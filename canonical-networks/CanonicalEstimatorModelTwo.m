function Blocking = CanonicalEstimatorModelTwo(m, c, deltaTilde, epsilonTilde, BurnInTime, BatchTime, MaxTime, Mu, Lambda) 

    %%%%%  Inputs  %%%%% Comment these out
%     clear all; clc;
%     ScenarioNumber = 1;
%     Capacity = [10,10,10,10,10,10];
%     CIgap = 0.05;
%     BurnInTime = 10;
%     BatchTime = 10;

    %%%%% Outputs %%%%%%  Do not comment out
%     Blocking = zeros(12,3); % This is per-route per-station blocking

    %%%%%%%%%%%%%%%%%%%
    %Mu = 0.75*rand(1,6); % This vector contains mean service times.  
%     Lambda = [1/30, 1/15]; % This vector contains mean interarrival times. We arbitrarily assume interarrival times for station 1 are 1/30 and for station 2 we assume 1/15
    CoVs = csvread(strcat('Scenarios/Scenario',num2str(m),'.csv'));
    P = 1./(2.*CoVs.^2);
    %%%% P = [Web1, Web2, App1, App2, Dbs1, Dbs2, Arr1, Arr2];
%     c = 10*rand(1,6);
%     c(1,1) = 0;
    %%%% Each row corresponds to a route (there are 12 routes) and says the stations that are visited
    %%%% by that route. 

    Routes = [1, 7, 7, 7;
              1, 3, 7, 7;
              1, 3, 5, 7;
              1, 4, 7, 7;
              1, 4, 5, 7;
              1, 4, 6, 7;
              2, 7, 7, 7;
              2, 4, 7, 7;
              2, 4, 6, 7;
              2, 3, 7, 7;
              2, 3, 6, 7;
              2, 3, 5, 7];
          
    %%%% note: station 7 indicates the job leaves the network

    %%%%% Define anonymous functions to govern interarrival and service times:
    InterArrival1 = @() Lambda(1,1)*(-log(rand)-(rand<P(1,7))*log(rand)/(P(1,7)))/2;
    InterArrival2 = @() Lambda(1,2)*(-log(rand)-(rand<P(1,8))*log(rand)/(P(1,8)))/2;
    Service1 = @() Mu(1,1)*(-log(rand)-(rand<P(1,1))*log(rand)/(P(1,1)))/2;
    Service2 = @() Mu(1,2)*(-log(rand)-(rand<P(1,2))*log(rand)/(P(1,2)))/2;
    % Service3 = @() ;
    % Service4 = @() Mu(1,4)*(-log(rand)-(rand<P(1,4))*log(rand)/(P(1,4)))/2;
    % Service5 = @() Mu(1,5)*(-log(rand)-(rand<P(1,5))*log(rand)/(P(1,5)))/2;
    % Service6 = @() Mu(1,6)*(-log(rand)-(rand<P(1,6))*log(rand)/(P(1,6)))/2;

    %%%% Three types of events: (1) arrival to first station, (2) arrival to second station, (3) service. 
    %%%% For arrivals we specify time and location
    %%%% For services we specify time, route (1-12), progression on route (1,2,3)

    %%%% Initilise %%%%%

%     ArrivalsToStations = zeros(1,6);
%     BlockagesAtStations = zeros(1,6);
    Arrivals = [0, NaN, NaN;
                0, 0, NaN;
                0, 0, 0;
                0, 0, NaN;
                0, 0, 0;
                0, 0, 0;
                0, NaN, NaN;
                0, 0, NaN;
                0, 0, 0;
                0, 0, NaN;
                0, 0, 0;
                0, 0, 0]; 
    Blockages = [0, NaN, NaN;
                 0, 0, NaN;
                 0, 0, 0;
                 0, 0, NaN;
                 0, 0, 0;
                 0, 0, 0;
                 0, NaN, NaN;
                 0, 0, NaN;
                 0, 0, 0;
                 0, 0, NaN;
                 0, 0, 0;
                 0, 0, 0];          

    c = [c, 0];
    floorC = floor(c);
    CurrentlyProcessing = zeros(1,7);

    Cremainders = c - floorC;
    EventList = ...
        [InterArrival1(), 1, NaN, NaN; 
         InterArrival2(), 2, NaN, NaN;
         inf, NaN, NaN, NaN];
    EventList = sortrows(EventList, 1); 
    CIsize = 1; 
    NextEventTime = EventList(1,1);
    %QueuePlot = CurrentSystemOccupancy;
    while NextEventTime <= BurnInTime
        EventType = EventList(1,2);
        switch EventType
            case 1 %%Arrival at station 1
                ArrivalType = randi(6);
                EntryThresholds = rand(1,4);      
                StationIndexOnRoute = find(CurrentlyProcessing(Routes(ArrivalType, :)) < floorC(Routes(ArrivalType, :)) | (CurrentlyProcessing(Routes(ArrivalType, :)) < c(Routes(ArrivalType, :)) & (EntryThresholds <= Cremainders(Routes(ArrivalType, :)))),1);                
                if StationIndexOnRoute == 1
%                     Arrivals(ArrivalType, 1) = Arrivals(ArrivalType, 1) + 1;
                    CurrentlyProcessing(Routes(ArrivalType,1)) = CurrentlyProcessing(Routes(ArrivalType,1))+1;
                    NewServiceTime = NextEventTime + Service1();
                    NewArrivalTime = NextEventTime + InterArrival1();
                    NewServicePosition = find(EventList(:,1) >= NewServiceTime,1);
                    NewArrivalPosition = find(EventList(:,1) >= NewArrivalTime,1);    
                    if NewServiceTime < NewArrivalTime
                        EventList = [EventList(2:(NewServicePosition-1),:);
                                     NewServiceTime, 3, Routes(ArrivalType,1), 1;
                                     EventList(NewServicePosition:(NewArrivalPosition-1),:);
                                     NewArrivalTime, 1, NaN, NaN;
                                     EventList(NewArrivalPosition:end,:)];
                    else
                        EventList = [EventList(2:(NewArrivalPosition-1),:);
                                     NewArrivalTime, 1, NaN, NaN;
                                     EventList(NewArrivalPosition:(NewServicePosition-1),:);
                                     NewServiceTime, 3, Routes(ArrivalType,1), 1;
                                     EventList(NewServicePosition:end,:)];
                    end                        
                elseif  StationIndexOnRoute == 2
%                     Arrivals(ArrivalType, 1:2) = Arrivals(ArrivalType, 1:2) + 1;
%                     Blockages(ArrivalType,1) = Blockages(ArrivalType,1);
                    CurrentlyProcessing(Routes(ArrivalType,2)) = CurrentlyProcessing(Routes(ArrivalType,2))+1;
                    NewServiceTime = NextEventTime + Service1();
                    NewArrivalTime = NextEventTime + InterArrival1();
                    NewServicePosition = find(EventList(:,1) >= NewServiceTime,1);
                    NewArrivalPosition = find(EventList(:,1) >= NewArrivalTime,1);    
                    if NewServiceTime < NewArrivalTime
                        EventList = [EventList(2:(NewServicePosition-1),:);
                                     NewServiceTime, 3, Routes(ArrivalType,2), 1;
                                     EventList(NewServicePosition:(NewArrivalPosition-1),:);
                                     NewArrivalTime, 1, NaN, NaN;
                                     EventList(NewArrivalPosition:end,:)];
                    else
                        EventList = [EventList(2:(NewArrivalPosition-1),:);
                                     NewArrivalTime, 1, NaN, NaN;
                                     EventList(NewArrivalPosition:(NewServicePosition-1),:);
                                     NewServiceTime, 3, Routes(ArrivalType,2), 1;
                                     EventList(NewServicePosition:end,:)];
                    end
                elseif StationIndexOnRoute == 3
%                     Arrivals(ArrivalType, 1:3) = Arrivals(ArrivalType, 1:3) + 1;
%                     Blockages(ArrivalType,1:2) = Blockages(ArrivalType,1:2);
                    CurrentlyProcessing(Routes(ArrivalType,3)) = CurrentlyProcessing(Routes(ArrivalType,3))+1;
                    NewServiceTime = NextEventTime + Service1();
                    NewArrivalTime = NextEventTime + InterArrival1();
                    NewServicePosition = find(EventList(:,1) >= NewServiceTime,1);
                    NewArrivalPosition = find(EventList(:,1) >= NewArrivalTime,1);    
                    if NewServiceTime < NewArrivalTime
                        EventList = [EventList(2:(NewServicePosition-1),:);
                                     NewServiceTime, 3, Routes(ArrivalType,3), 1;
                                     EventList(NewServicePosition:(NewArrivalPosition-1),:);
                                     NewArrivalTime, 1, NaN, NaN;
                                     EventList(NewArrivalPosition:end,:)];
                    else
                        EventList = [EventList(2:(NewArrivalPosition-1),:);
                                     NewArrivalTime, 1, NaN, NaN;
                                     EventList(NewArrivalPosition:(NewServicePosition-1),:);
                                     NewServiceTime, 3, Routes(ArrivalType,3), 1;
                                     EventList(NewServicePosition:end,:)];
                    end
                else
                    Arrivals(ArrivalType, 1:3) = Arrivals(ArrivalType, 1:3) + 1;
                    Blockages(ArrivalType,1:3) = Blockages(ArrivalType,1:3) + 1;
                    NewArrivalTime = NextEventTime + InterArrival1();
                    NewArrivalPosition = find(EventList(:,1) >= NewArrivalTime,1);    
                    EventList = [EventList(2:(NewArrivalPosition-1),:);
                                 NewArrivalTime, 1, NaN, NaN;
                                 EventList(NewArrivalPosition:end,:)];
                end
            case 2 %%Arrival at station 2
                ArrivalType = 6+randi(6);
                EntryThresholds = rand(1,4);      
                StationIndexOnRoute = find(CurrentlyProcessing(Routes(ArrivalType, :)) < floorC(Routes(ArrivalType, :)) | (CurrentlyProcessing(Routes(ArrivalType, :)) < c(Routes(ArrivalType, :)) & (EntryThresholds <= Cremainders(Routes(ArrivalType, :)))),1);                
                if StationIndexOnRoute == 1
%                     Arrivals(ArrivalType, 1) = Arrivals(ArrivalType, 1) + 1;
                        CurrentlyProcessing(Routes(ArrivalType,1)) = CurrentlyProcessing(Routes(ArrivalType,1))+1;
                        NewServiceTime = NextEventTime + Service2();
                        NewArrivalTime = NextEventTime + InterArrival2();
                        NewServicePosition = find(EventList(:,1) >= NewServiceTime,1);
                        NewArrivalPosition = find(EventList(:,1) >= NewArrivalTime,1);    
                    if NewServiceTime < NewArrivalTime
                        EventList = [EventList(2:(NewServicePosition-1),:);
                                     NewServiceTime, 3, Routes(ArrivalType,1), 1;
                                     EventList(NewServicePosition:(NewArrivalPosition-1),:);
                                     NewArrivalTime, 2, NaN, NaN;
                                     EventList(NewArrivalPosition:end,:)];
                    else
                        EventList = [EventList(2:(NewArrivalPosition-1),:);
                                     NewArrivalTime, 2, NaN, NaN;
                                     EventList(NewArrivalPosition:(NewServicePosition-1),:);
                                     NewServiceTime, 3, Routes(ArrivalType,1), 1;
                                     EventList(NewServicePosition:end,:)];
                    end                        
                elseif  StationIndexOnRoute == 2
%                     Arrivals(ArrivalType, 1:2) = Arrivals(ArrivalType, 1:2) + 1;
%                     Blockages(ArrivalType,1) = Blockages(ArrivalType,1);
                    CurrentlyProcessing(Routes(ArrivalType,2)) = CurrentlyProcessing(Routes(ArrivalType,2))+1;
                    NewServiceTime = NextEventTime + Service2();
                    NewArrivalTime = NextEventTime + InterArrival2();
                    NewServicePosition = find(EventList(:,1) >= NewServiceTime,1);
                    NewArrivalPosition = find(EventList(:,1) >= NewArrivalTime,1);    
                    if NewServiceTime < NewArrivalTime
                        EventList = [EventList(2:(NewServicePosition-1),:);
                                     NewServiceTime, 3, Routes(ArrivalType,2), 1;
                                     EventList(NewServicePosition:(NewArrivalPosition-1),:);
                                     NewArrivalTime, 2, NaN, NaN;
                                     EventList(NewArrivalPosition:end,:)];
                    else
                        EventList = [EventList(2:(NewArrivalPosition-1),:);
                                     NewArrivalTime, 2, NaN, NaN;
                                     EventList(NewArrivalPosition:(NewServicePosition-1),:);
                                     NewServiceTime, 3, Routes(ArrivalType,2), 1;
                                     EventList(NewServicePosition:end,:)];
                    end
                elseif StationIndexOnRoute == 3
%                     Arrivals(ArrivalType, 1:3) = Arrivals(ArrivalType, 1:3) + 1;
%                     Blockages(ArrivalType,1:2) = Blockages(ArrivalType,1:2);
                    CurrentlyProcessing(Routes(ArrivalType,3)) = CurrentlyProcessing(Routes(ArrivalType,3))+1;
                    NewServiceTime = NextEventTime + Service2();
                    NewArrivalTime = NextEventTime + InterArrival2();
                    NewServicePosition = find(EventList(:,1) >= NewServiceTime,1);
                    NewArrivalPosition = find(EventList(:,1) >= NewArrivalTime,1);    
                    if NewServiceTime < NewArrivalTime
                        EventList = [EventList(2:(NewServicePosition-1),:);
                                     NewServiceTime, 3, Routes(ArrivalType,3), 1;
                                     EventList(NewServicePosition:(NewArrivalPosition-1),:);
                                     NewArrivalTime, 2, NaN, NaN;
                                     EventList(NewArrivalPosition:end,:)];
                    else
                        EventList = [EventList(2:(NewArrivalPosition-1),:);
                                     NewArrivalTime, 2, NaN, NaN;
                                     EventList(NewArrivalPosition:(NewServicePosition-1),:);
                                     NewServiceTime, 3, Routes(ArrivalType,3), 1;
                                     EventList(NewServicePosition:end,:)];
                    end
                else
                    Arrivals(ArrivalType, 1:3) = Arrivals(ArrivalType, 1:3) + 1;
                    Blockages(ArrivalType,1:3) = Blockages(ArrivalType,1:3) + 1;
                    NewArrivalTime = NextEventTime + InterArrival2();
                    NewArrivalPosition = find(EventList(:,1) >= NewArrivalTime,1);    
                    EventList = [EventList(2:(NewArrivalPosition-1),:);
                                 NewArrivalTime, 2, NaN, NaN;
                                 EventList(NewArrivalPosition:end,:)];
                end
            case 3 %%Service
                StationOfService = EventList(1,3);
                CurrentlyProcessing(StationOfService) = CurrentlyProcessing(StationOfService)-1;
                EventList(1,:) = [];
        end
        NextEventTime = EventList(1,1);
    %     QueuePlot = [QueuePlot; CurrentSystemOccupancy];
    end   

    % TimeSteps = length(QueuePlot);
    % for i = 1:6
    %     subplot(6,1,i)
    %     stairs(1:TimeSteps, QueuePlot(:,i))
    % end

    while CIsize > deltaTilde  && NextEventTime < MaxTime
        ThisBatchMaxTime = NextEventTime + BatchTime;
        while NextEventTime <= ThisBatchMaxTime 
            EventType = EventList(1,2);
            switch EventType
                case 1 %%Arrival at station 1
                    ArrivalType = randi(6);
                    EntryThresholds = rand(1,4);      
                    StationIndexOnRoute = find(CurrentlyProcessing(Routes(ArrivalType, :)) < floorC(Routes(ArrivalType, :)) | (CurrentlyProcessing(Routes(ArrivalType, :)) < c(Routes(ArrivalType, :)) & (EntryThresholds <= Cremainders(Routes(ArrivalType, :)))),1);                
                    if StationIndexOnRoute == 1
                        Arrivals(ArrivalType, 1) = Arrivals(ArrivalType, 1) + 1;
                            CurrentlyProcessing(Routes(ArrivalType,1)) = CurrentlyProcessing(Routes(ArrivalType,1))+1;
                            NewServiceTime = NextEventTime + Service1();
                            NewArrivalTime = NextEventTime + InterArrival1();
                            NewServicePosition = find(EventList(:,1) >= NewServiceTime,1);
                            NewArrivalPosition = find(EventList(:,1) >= NewArrivalTime,1);    
                        if NewServiceTime < NewArrivalTime
                            EventList = [EventList(2:(NewServicePosition-1),:);
                                         NewServiceTime, 3, Routes(ArrivalType,1), 1;
                                         EventList(NewServicePosition:(NewArrivalPosition-1),:);
                                         NewArrivalTime, 1, NaN, NaN;
                                         EventList(NewArrivalPosition:end,:)];
                        else
                            EventList = [EventList(2:(NewArrivalPosition-1),:);
                                         NewArrivalTime, 1, NaN, NaN;
                                         EventList(NewArrivalPosition:(NewServicePosition-1),:);
                                         NewServiceTime, 3, Routes(ArrivalType,1), 1;
                                         EventList(NewServicePosition:end,:)];
                        end                        
                    elseif  StationIndexOnRoute == 2
                        Arrivals(ArrivalType, 1:2) = Arrivals(ArrivalType, 1:2) + 1;
                        Blockages(ArrivalType,1) = Blockages(ArrivalType,1)+1;
                        CurrentlyProcessing(Routes(ArrivalType,2)) = CurrentlyProcessing(Routes(ArrivalType,2))+1;
                        NewServiceTime = NextEventTime + Service1();
                        NewArrivalTime = NextEventTime + InterArrival1();
                        NewServicePosition = find(EventList(:,1) >= NewServiceTime,1);
                        NewArrivalPosition = find(EventList(:,1) >= NewArrivalTime,1);    
                        if NewServiceTime < NewArrivalTime
                            EventList = [EventList(2:(NewServicePosition-1),:);
                                         NewServiceTime, 3, Routes(ArrivalType,2), 1;
                                         EventList(NewServicePosition:(NewArrivalPosition-1),:);
                                         NewArrivalTime, 1, NaN, NaN;
                                         EventList(NewArrivalPosition:end,:)];
                        else
                            EventList = [EventList(2:(NewArrivalPosition-1),:);
                                         NewArrivalTime, 1, NaN, NaN;
                                         EventList(NewArrivalPosition:(NewServicePosition-1),:);
                                         NewServiceTime, 3, Routes(ArrivalType,2), 1;
                                         EventList(NewServicePosition:end,:)];
                        end
                    elseif StationIndexOnRoute == 3
                        Arrivals(ArrivalType, 1:3) = Arrivals(ArrivalType, 1:3) + 1;
                        Blockages(ArrivalType,1:2) = Blockages(ArrivalType,1:2) + 1;
                        CurrentlyProcessing(Routes(ArrivalType,3)) = CurrentlyProcessing(Routes(ArrivalType,3))+1;
                        NewServiceTime = NextEventTime + Service1();
                        NewArrivalTime = NextEventTime + InterArrival1();
                        NewServicePosition = find(EventList(:,1) >= NewServiceTime,1);
                        NewArrivalPosition = find(EventList(:,1) >= NewArrivalTime,1);    
                        if NewServiceTime < NewArrivalTime
                            EventList = [EventList(2:(NewServicePosition-1),:);
                                         NewServiceTime, 3, Routes(ArrivalType,3), 1;
                                         EventList(NewServicePosition:(NewArrivalPosition-1),:);
                                         NewArrivalTime, 1, NaN, NaN;
                                         EventList(NewArrivalPosition:end,:)];
                        else
                            EventList = [EventList(2:(NewArrivalPosition-1),:);
                                         NewArrivalTime, 1, NaN, NaN;
                                         EventList(NewArrivalPosition:(NewServicePosition-1),:);
                                         NewServiceTime, 3, Routes(ArrivalType,3), 1;
                                         EventList(NewServicePosition:end,:)];
                        end
                    else
                        Arrivals(ArrivalType, 1:3) = Arrivals(ArrivalType, 1:3) + 1;
                        Blockages(ArrivalType,1:3) = Blockages(ArrivalType,1:3) + 1;
                        NewArrivalTime = NextEventTime + InterArrival1();
                        NewArrivalPosition = find(EventList(:,1) >= NewArrivalTime,1);    
                        EventList = [EventList(2:(NewArrivalPosition-1),:);
                                     NewArrivalTime, 1, NaN, NaN;
                                     EventList(NewArrivalPosition:end,:)];
                    end
                case 2 %%Arrival at station 2
                    ArrivalType = 6+randi(6);
                    EntryThresholds = rand(1,4);      
                    StationIndexOnRoute = find(CurrentlyProcessing(Routes(ArrivalType, :)) < floorC(Routes(ArrivalType, :)) | (CurrentlyProcessing(Routes(ArrivalType, :)) < c(Routes(ArrivalType, :)) & (EntryThresholds <= Cremainders(Routes(ArrivalType, :)))),1);            
                    if StationIndexOnRoute == 1
                        Arrivals(ArrivalType, 1) = Arrivals(ArrivalType, 1) + 1;
                            CurrentlyProcessing(Routes(ArrivalType,1)) = CurrentlyProcessing(Routes(ArrivalType,1))+1;
                            NewServiceTime = NextEventTime + Service2();
                            NewArrivalTime = NextEventTime + InterArrival2();
                            NewServicePosition = find(EventList(:,1) >= NewServiceTime,1);
                            NewArrivalPosition = find(EventList(:,1) >= NewArrivalTime,1);    
                        if NewServiceTime < NewArrivalTime
                            EventList = [EventList(2:(NewServicePosition-1),:);
                                         NewServiceTime, 3, Routes(ArrivalType,1), 1;
                                         EventList(NewServicePosition:(NewArrivalPosition-1),:);
                                         NewArrivalTime, 2, NaN, NaN;
                                         EventList(NewArrivalPosition:end,:)];
                        else
                            EventList = [EventList(2:(NewArrivalPosition-1),:);
                                         NewArrivalTime, 2, NaN, NaN;
                                         EventList(NewArrivalPosition:(NewServicePosition-1),:);
                                         NewServiceTime, 3, Routes(ArrivalType,1), 1;
                                         EventList(NewServicePosition:end,:)];
                        end                        
                    elseif  StationIndexOnRoute == 2
                        Arrivals(ArrivalType, 1:2) = Arrivals(ArrivalType, 1:2) + 1;
                        Blockages(ArrivalType,1) = Blockages(ArrivalType,1) + 1;
                        CurrentlyProcessing(Routes(ArrivalType,2)) = CurrentlyProcessing(Routes(ArrivalType,2))+1;
                        NewServiceTime = NextEventTime + Service2();
                        NewArrivalTime = NextEventTime + InterArrival2();
                        NewServicePosition = find(EventList(:,1) >= NewServiceTime,1);
                        NewArrivalPosition = find(EventList(:,1) >= NewArrivalTime,1);    
                        if NewServiceTime < NewArrivalTime
                            EventList = [EventList(2:(NewServicePosition-1),:);
                                         NewServiceTime, 3, Routes(ArrivalType,2), 1;
                                         EventList(NewServicePosition:(NewArrivalPosition-1),:);
                                         NewArrivalTime, 2, NaN, NaN;
                                         EventList(NewArrivalPosition:end,:)];
                        else
                            EventList = [EventList(2:(NewArrivalPosition-1),:);
                                         NewArrivalTime, 2, NaN, NaN;
                                         EventList(NewArrivalPosition:(NewServicePosition-1),:);
                                         NewServiceTime, 3, Routes(ArrivalType,2), 1;
                                         EventList(NewServicePosition:end,:)];
                        end
                    elseif StationIndexOnRoute == 3
                        Arrivals(ArrivalType, 1:3) = Arrivals(ArrivalType, 1:3) + 1;
                        Blockages(ArrivalType,1:2) = Blockages(ArrivalType,1:2) + 1;
                        CurrentlyProcessing(Routes(ArrivalType,3)) = CurrentlyProcessing(Routes(ArrivalType,3))+1;
                        NewServiceTime = NextEventTime + Service2();
                        NewArrivalTime = NextEventTime + InterArrival2();
                        NewServicePosition = find(EventList(:,1) >= NewServiceTime,1);
                        NewArrivalPosition = find(EventList(:,1) >= NewArrivalTime,1);    
                        if NewServiceTime < NewArrivalTime
                            EventList = [EventList(2:(NewServicePosition-1),:);
                                         NewServiceTime, 3, Routes(ArrivalType,3), 1;
                                         EventList(NewServicePosition:(NewArrivalPosition-1),:);
                                         NewArrivalTime, 2, NaN, NaN;
                                         EventList(NewArrivalPosition:end,:)];
                        else
                            EventList = [EventList(2:(NewArrivalPosition-1),:);
                                         NewArrivalTime, 2, NaN, NaN;
                                         EventList(NewArrivalPosition:(NewServicePosition-1),:);
                                         NewServiceTime, 3, Routes(ArrivalType,3), 1;
                                         EventList(NewServicePosition:end,:)];
                        end
                    else
                        Arrivals(ArrivalType, 1:3) = Arrivals(ArrivalType, 1:3) + 1;
                        Blockages(ArrivalType,1:3) = Blockages(ArrivalType,1:3) + 1;
                        NewArrivalTime = NextEventTime + InterArrival2();
                        NewArrivalPosition = find(EventList(:,1) >= NewArrivalTime,1);    
                        EventList = [EventList(2:(NewArrivalPosition-1),:);
                                     NewArrivalTime, 2, NaN, NaN;
                                     EventList(NewArrivalPosition:end,:)];
                    end
                case 3 %%Service
                    StationOfService = EventList(1,3);
                    CurrentlyProcessing(StationOfService) = CurrentlyProcessing(StationOfService)-1;
                    EventList(1,:) = [];
            end
            NextEventTime = EventList(1,1);
        end
        Blocking = Blockages./Arrivals;
        CI = 1.96.*sqrt(1.96.^2-4.*Arrivals.*(Blocking-1).*Blocking)./(1.96.^2+Arrivals);
        CIsize = max(max(CI));
    %     QueuePlot = [QueuePlot; CurrentSystemOccupancy];
    end 
    for i = 1:12
        for j = 1:3
            if c(Routes(i,j)) < 1 && Routes(i,j) < 7 &&  Arrivals(i,j) == 0
                Blocking(i,j) = 1;
            elseif Arrivals(i,j) == 0 && min(Blocking( Routes == Routes(i,j)))>0
                tempBlock = Blocking(Routes == Routes(i,j)).*(Arrivals(Routes == Routes(i,j))>0);
                Blocking(i,j) = mean(tempBlock(tempBlock>0));
            end
        end
    end    
    Blocking = max( Blocking, epsilonTilde*(floorC(Routes(1:end, 1:3))>0));
end
