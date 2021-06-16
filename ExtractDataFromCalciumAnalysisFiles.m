%%Code usd to extract data from "[filename] intensity.mat"
%Load "[file name] intensity data.mat" into appropriate directory below
%%% Only run either retroKD or No DG inj code

clear all
includeZero = true; % variable used to include non-firing cells or not into graphs


%% For retrokKD
saline = {...
'C:\Data\Saline M1 intensity data.mat',...
'C:\Data\Saline M2 intensity data.mat',...
'C:\Data\Saline M3 intensity data.mat',...
'C:\Data\Saline M4 intensity data.mat',...
'C:\Data\Saline M5 intensity data.mat',...
'C:\Data\Saline M6 intensity data.mat',...
};

virus = {...
'C:\Data\Virus M1 intensity data.mat',...
'C:\Data\Virus M2 intensity data.mat',...
'C:\Data\Virus M3 intensity data.mat',...
'C:\Data\Virus M4 intensity data.mat',...
'C:\Data\Virus M5 intensity data.mat',...
};

% %% No DG injection
% saline = {...
%     'C:\Data\RetroOnly M1 intensity data.mat',...
%     'C:\Data\RetroOnly M2 intensity data.mat',...
%     'C:\Data\RetroOnly M3 intensity data.mat',...
%     'C:\Data\RetroOnly M4 intensity data.mat',...
%     'C:\Data\RetroOnly M5 intensity data.mat',...
%     'C:\Data\RetroOnly M6 intensity data.mat',...
%     };
% virus = {};
%%


edges = 0:.5:20;
edges = [edges 100]; %just to make sure we don't miss any extremely active cells
xVar = edges(2:end-1);
if includeZero
    xVar = [0 xVar];
end

salineTotCells = [];
virusTotCells  = [];

salineMeanFiringPerAnimal =[];
salineMeanPercentActivePerAnimal = [];
virusMeanFiringPerAnimal = [];
virusMeanPercentActivePerAnimal = [];
salineMeanActivePerAnimal =[];
virusMeanActivePerAnimal=[];
salineMeanTotalEventsPerAnimal =[];
virusMeanTotalEventsPerAnimal =[];


CDFFreqAll = [];
normCDFFreqAll = [];


for i=1:2
    testVar = i;

    if (testVar == 1)
        data2Use = saline;
    else 
        data2Use = virus;
    end


    for ii= 1:length(data2Use)
        directory = string(data2Use(ii));
        load(directory);

    %Use Alex's binary firing for raster
        IntensityRaster = binaryFiring;



    binSize = ceil(3.91*60); %size of 1 minute bin
    continuousFreqAverage =[];
    for iii = 1:size(IntensityRaster,2) - binSize +1
        continuousFreqAverage = [continuousFreqAverage sum(IntensityRaster(:,iii:iii+binSize-1),2)];
    end
    % This give the average firing rate per minute for each animal, this
    % verctor is lost each iteration
        EventsPerMinutePerCell = mean(continuousFreqAverage,2);
        EventsPerMinutePerAnimal = sum(EventsPerMinutePerCell);
        MeanFiringRatePerAnimal = mean(EventsPerMinutePerCell);

        distFreq = histcounts(EventsPerMinutePerCell,edges);
        if includeZero
            distFreq = [size(intensityData,1) distFreq];
        end
        
            CDFFreq = cumsum(distFreq);
                CDFFreqAll = [ CDFFreqAll ; CDFFreq ];
                normCDFFreq = CDFFreq./CDFFreq(end);
                    normCDFFreqAll = [ normCDFFreqAll ; normCDFFreq ];
        
        
     % Percent active per minute 
     %(frequency of firing doesn't matter, only if it is active in any minute
        isActive = continuousFreqAverage >0;
        CellsActivePerMinuteVect = sum(isActive,1);
            MeanCellsActivePerMinute = mean(CellsActivePerMinuteVect);
        PercentCellsActivePerMinuteVect = CellsActivePerMinuteVect./ size(continuousFreqAverage,1);
            MeanPercentCellsActivePerMinute = mean(PercentCellsActivePerMinuteVect);
                                                    
    %Make all numbers accessible after running        
    if testVar ==1
        salineMeanFiringPerAnimal = [salineMeanFiringPerAnimal MeanFiringRatePerAnimal];
        salineMeanActivePerAnimal = [salineMeanActivePerAnimal MeanCellsActivePerMinute];
        salineMeanPercentActivePerAnimal = [salineMeanPercentActivePerAnimal MeanPercentCellsActivePerMinute];
        salineMeanTotalEventsPerAnimal = [salineMeanTotalEventsPerAnimal EventsPerMinutePerAnimal];
        salineTotCells = [salineTotCells size(intensityData,1)];
    elseif testVar ==2
        virusMeanFiringPerAnimal = [virusMeanFiringPerAnimal   MeanFiringRatePerAnimal];
        virusMeanActivePerAnimal = [virusMeanActivePerAnimal MeanCellsActivePerMinute];
        virusMeanPercentActivePerAnimal = [virusMeanPercentActivePerAnimal   MeanPercentCellsActivePerMinute];
        virusMeanTotalEventsPerAnimal = [virusMeanTotalEventsPerAnimal EventsPerMinutePerAnimal];
        virusTotCells = [virusTotCells size(intensityData,1)];
    end



    end

end


% plots: NormalizedFreqDist
figure;
hold on
data = normCDFFreqAll;
for i=1:2
    testVar = i;
    if testVar == 1
        for ii = 1:length(saline)
            plot(xVar, data(ii,1:end-1),'b')
        end
    elseif testVar ==2
        for ii = 1:length(virus)
            plot(xVar, data(ii+length(saline),1:end-1),'r')
        end
    end
end
plot(xVar, mean(data(1:length(saline),1:end-1)),'b','LineWidth',2)
plot(xVar, mean(data(length(saline):end,1:end-1)),'r','LineWidth',2)
title('Cumulative Distribution');
xlabel('events per minute per cell');
ylabel('propotion of cells');

hold off

% plots: NormalizedFreqDist
figure;
hold on
data = CDFFreqAll;
for i=1:2
    testVar = i;
    if testVar == 1
        for ii = 1:length(saline)
            plot(xVar, data(ii,1:end-1),'b')
        end
    elseif testVar ==2
        for ii = 1:length(virus)
            plot(xVar, data(ii+length(saline),1:end-1),'r')
        end
    end
end
plot(xVar, mean(data(1:length(saline),1:end-1)),'b','LineWidth',2)
plot(xVar, mean(data(length(saline):end,1:end-1)),'r','LineWidth',2)
title('Cumulative Distribution');
xlabel('events per minute per cell');
ylabel('propotion of cells');

hold off



ForGraphPad_TotCells = [salineTotCells virusTotCells];

%data=normCDFFreqAll;
%ForGraphPad_NormCDF = transpose([xVar ; 100.*mean(data(1:length(saline),1:end-1)); 100.*mean(data(length(saline):end,1:end-1))]);

%data=CDFFreqAll;
%ForGraphPad_CDF = transpose([xVar ; mean(data(1:length(saline),1:end-1)); mean(data(length(saline):end,1:end-1))]);

ForGraphPad_MeanCellRatePerAnimal = [salineMeanFiringPerAnimal virusMeanFiringPerAnimal];
ForGraphPad_MeanPercentCellActvePerAnimal = [salineMeanPercentActivePerAnimal virusMeanPercentActivePerAnimal];
ForGraphPad_MeanCellActivePerAnimal = [salineMeanActivePerAnimal virusMeanActivePerAnimal];
ForGraphPad_MeanEventsPerMinutePerAnimal = [salineMeanTotalEventsPerAnimal virusMeanTotalEventsPerAnimal];
ForGraphPad_MeanPercentCellActivePerAnimalPerMinute = ForGraphPad_MeanCellActivePerAnimal ./ ForGraphPad_TotCells * 100;