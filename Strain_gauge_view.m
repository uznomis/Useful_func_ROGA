% this program loads data produced by ROGA Terminal developed by James
% Jeffers, display the data according to the user request
%% Initializing
% paramters to change before continuing
cardSN = [1,2];    % which cards to import
chNames = {'J3','J2','J1','I3','I2','I1','H3','H2','H1',...
    'G3','G2','G1','F3','F2','F1','Encoder','Unused';...
    'E3','E2','E1','D3','D2','D1','C3','C2','C1',...
    'B3','B2','B1','A3','A2','A1','Accel','Encoder'};
freq = 5e5;    % freqeuncy in Hz
cardToGetEncoder = 1;    % the number of card whose encoder data is used for velocity and counter calculation; put 0 if you don't want velocity and counter
encoderChannels = [16;17];    % encoder channels for the two cards; should match the ordering in cardSN
encoderAvailable = 1;    % 1 for yes, 0 for no; if no encoder is available, there is no sync between cards
lenPerSector = 1.5e-6;    % encoder sector length in meters
baseLevelFactors = [1,1];
baseLevelOffsets = [2.034,2.034];
defaultBaseVoltages = ones(2,15);
% defaultBaseVoltages = [;];    % custom V0

%% Importing
filename = cell(1,length(cardSN));
filepath = cell(1,length(cardSN));
for i = 1:length(cardSN)
    [filename{i}, filepath{i}] = uigetfile...
        ('C:\Users\User\Desktop\ROGA Terminal\data\*.*',...
        ['Please select ACQ100', num2str(cardSN(i))]);
    if filename{i} == 0
        return
    end
end
test = cell(1,length(cardSN));
raw_test = test;
for i = 1:length(cardSN)
    test{i} = dlmread([filepath{i}, filename{i}]);
    test{i} = reshape(test{i}, [length(test{i})/17 17]);
    raw_test{i} = test{i};
    test{i} = test{i} * baseLevelFactors(i) + baseLevelOffsets(i);
    test{i}(:,16:17) = test{i}(:,16:17) -  baseLevelOffsets(i);
end

%% Velocity and counter calculation
% calculate s and v from encoder; calling another function
if cardToGetEncoder == 0
    s = zeros(1,length(test{1}));
    v = s;
    encoderInd = 1;
else
    encoderInd = find(cardSN == cardToGetEncoder, 1);
    [s,v] = getCounterFromEncoder(test{encoderInd}(:,...
        encoderChannels(encoderInd)),freq,lenPerSector);
end

%% Sync cards
% parameters regarding synchronization; no need to change unless sync is
% failed; if failed, observe the encoder pattern: if the encoder is sparse
% then increase encoderSpacingThrsh and decrease otherwise, and if there is
% a warning saying matrix exceeds dimension then decrease
% encoderXcorrWindow.
encoderShift = 0.5;    % encoder is usually 0/4.6 volts, so pick in between
encoderSpacingThrsh = 5e2;
encoderXcorrWindow = 5e5;

% execution begins here
if ~isempty(filename)
    % calc timeshift of cards
    encodercurrent = 1;
    encoderbefore = 0;
    encoderjump = 1;
    spacing = 0;
    data1 = test{1};
    encoderjumpavg = zeros(1,length(test) - 1);
    if encoderAvailable == 1
        encoderCh1 = encoderChannels(1);
        for j = 2:length(data1)
            if (data1(j-1,encoderCh1) - encoderShift)*...
                    (data1(j,encoderCh1) - encoderShift) < 0
                encoderbefore = encodercurrent;
                encodercurrent = j;
                spacing = encodercurrent - encoderbefore;
                if spacing < encoderSpacingThrsh && encoderbefore > 1
                    encoderjump = encodercurrent;
                    break;
                end
            end
        end
        s1 = data1(encoderjump:encoderjump + encoderXcorrWindow,encoderCh1);
        
        for i=2:length(test)
            s2 = test{i}(encoderjump: encoderjump + encoderXcorrWindow,encoderChannels(i));
            [r,lag] = xcorr(s1,s2);
            [~,I] = max(abs(r));
            encoderjumpavg(i-1) = round(-lag(I));
        end
    end
    % create time col vector
    timecell = cell(1,length(test));
    timecell{1} = 1/freq * (1:length(data1))';
    
    for i=2:length(test)
        timecell{i} = ...
            timecell{1}-encoderjumpavg(i-1)/freq; % apply time shift
    end
    
end

%% Plotting
% can display data median-filtered, smoothened and demeaned, as well as
% scaled. After the above sections are run, just simply run the section to
% plot different versions of data plots catered to your needs.

% parameters to change before executing the following code
cardToShow = [1,2];    % cards to display
% chSN = [16,17;16,17];
% chSN = [6,5,4,3,2,1];
% chSN = [1:17;1:17];
 chSN = [1:15;1:15];    % for velocity field picking
% chSN = [3,6,9;3,6,9];
% chSN = [15,14,13;15,14,13]; % first gage (shear)
% chSN = [13,10,7,4,1,16;13,10,7,4,1,16]; % third gage (shear)
% chSN = [16;16];    % channels to display on each card
% chSN = [14,11,8,5,2;14,11,8,5,2];    % normal load channels to display on each card
% chSN = [8,5,2;8,5,2];
% chSN = [2,5,8,11,14;2,5,8,11,14];    % channels to display on each card
% chSN = [1,3,4,6,7,9,10,12,13,15;1,3,4,6,7,9,10,12,13,15];    % channels to display on each card
smoothSpan = [5,5];    % smoothening window; make it 1 to diable smoothening; note there is no smoothening for AE data
medFilterOn = [1,1];    % choosing 1 makes median filter on; median filter gets rid of spikes
cardOffset = 5;    % offset between cards when plotting
chOffset = 0.075;    % offset between channels within each card when plotting
cardAmp = [10,10];
encoderVelDis = 0;    % velocity & distance from encoder; 0 for not plotting; 1 for only velocity; 2 for v & d
sSlope = 2e2;    % slope of counter for scaling when plotting
vSlope = 2e2;    % slope of velocity for scaling when plotting
encoderSlope = [0.01,0.1];
accelCard = 2;    % only support one card
accelCh = 16;    % only support one channel
accelAmp = 3;
plotWithPeaksAligned = 1;
scalePeaks = 0;
baseVoltageSpan = 100;

% execution begins here
figureTitle = '';
for i = 1:length(filename)
    figureTitle = [figureTitle,' ',filename{i}];
end
figure('Name',figureTitle);
hold on;

lineHandles = cell(length(cardToShow),length(chSN));

for j = 1:length(cardToShow)
    tempInd = find(cardSN == cardToShow(j),1);
    for i = 1:length(chSN(j,:))
        % first do median filter
        if medFilterOn(cardToShow(j))
            tempCh = medfilt2(test{tempInd}(:,chSN(j,i)),[3,1]);
        else
            tempCh = test{tempInd}(:,chSN(j,i));
        end
        
        tempCh(1) = tempCh(2); % for getting rid of the anomaly first point
        
        % apply amplification to individual channel
        if chSN(j,i) == encoderChannels(tempInd)
            tempCh = encoderSlope(j)*tempCh;
        elseif cardToShow(j) == accelCard && chSN(j,i) == accelCh
            tempCh = tempCh * accelAmp;
        else
            tempCh = tempCh * cardAmp(j);
            tempCh = smooth(tempCh, smoothSpan(j));
        end
        
        % plot
        if exist('manualOffsetIndices','var') && ...
                isequal(size(manualOffsetIndices),size(lineHandles)) && ...
                plotWithPeaksAligned
            if length(manualOffsetIndices{j,i}) ~= 2
                lineHandles{j,i} = plot(timecell{tempInd},...
                    j*cardOffset + i*chOffset + tempCh - mean(tempCh));
            else
                xOffset = timecell{tempInd}(manualOffsetIndices{j,i}(2));
                yOffset = tempCh(manualOffsetIndices{j,i}(1));
                scaleValue = tempCh(manualOffsetIndices{j,i}(2)) - ...
                    tempCh(manualOffsetIndices{j,i}(1));
                lineHandles{j,i} = plot(timecell{tempInd} - xOffset,...
                    (tempCh - yOffset)/abs(scaleValue^scalePeaks));
                % register arrival times and base voltages
                arrivalTimes{j,i} = xOffset;
                baseVoltages{j,i} = mean(test{tempInd}(...
                    manualOffsetIndices{j,i}(1)-baseVoltageSpan/2:...
                    manualOffsetIndices{j,i}(1)+baseVoltageSpan/2,...
                    chSN(j,i)));
                baseVoltageIndices{j,i} = manualOffsetIndices{j,i}(1);
            end
        else
            lineHandles{j,i} = plot(timecell{tempInd},...
                j*cardOffset + i*chOffset + tempCh - mean(tempCh));
        end
    end
end

strings = cell(length(cardToShow)*length(chSN) + encoderVelDis,1);

for j = 1:length(cardToShow)
    tempInd = find(cardSN == cardToShow(j),1);
    for i = 1:length(chSN)
        strings{(j - 1)*length(chSN) + i} ...
            = sprintf([num2str(cardToShow(j)),' Ch ',num2str(chSN(j,i)),...
            ' ',chNames{tempInd,chSN(j,i)}]);
    end
end

if encoderVelDis == 1
    plot(timecell{encoderInd},v*vSlope);
    strings{end} = 'velocity';
elseif encoderVelDis == 2
    plot(timecell{encoderInd},v*vSlope);
    plot(timecell{encoderInd},s*sSlope);
    strings{end} = 'counter';
    strings{end-1} = 'velocity';
end

legend (strings);

xlabel ('Time in seconds');
if encoderVelDis == 2
    ylabel (['Velocity * ',num2str(1/vSlope,'%.1e'),...
        ' m/s',newline,'Distance * ',num2str(1/sSlope,'%.1e'),' m']);
end

hold off;
return

%% Start new pick set
picks = [];

%% Pick up points
% go to the figure that you want to pick on, pick, and hit Ctrl + Enter.
try
    dcm_obj = datacursormode(gcf);
    c_info = getCursorInfo(dcm_obj);
    pick = zeros(length(c_info), 2);
    for i = 1:length(c_info)
        pick(i,:) = c_info(i).Position;
    end
    pick = flipud(pick);
catch ME
    error('Please repick.');
end
picks = [picks; pick];

% display picks
format long
disp(picks);

%% Export to file
% export the processed data to a text format

% saving gauge data set to excel file
xRange = xlim;
avg0 = mean(cell2mat(reshape(arrivalTimes,1,numel(arrivalTimes))));
if ~isnan(avg0)
    xRange = xRange + avg0*plotWithPeaksAligned;
end
dt = datestr(now,'mmmm_dd_yyyy_HH_MM_SS');
[~,name,~] = fileparts(filename{1});
for i = 1:length(filename)
    odata = [timecell{i} test{i}];
    leftInd = find(odata(:,1) > xRange(1),1,'first');
    rightInd = find(odata(:,1) < xRange(2),1,'last');
    if rightInd - leftInd > 500000
        msgbox(['too much data to save: ',...
            num2str(rightInd - leftInd),' lines']);
        return
    end
    headers = [{'Time'},chNames(i,:)];
    xlswrite([filepath{1},name,' ',dt,'.xlsx'], headers, [filename{i}(1:10),num2str(i)], 'A1');
    xlswrite([filepath{1},name,' ',dt,'.xlsx'], odata(leftInd:rightInd,:), [filename{i}(1:10),num2str(i)], 'A2');
end

% save encoder velocity and distance
odata = [timecell{encoderInd} v' s'];
headers = {'Time','Velocity','Distance'};
leftInd = find(odata(:,1) > xRange(1),1,'first');
rightInd = find(odata(:,1) < xRange(2),1,'last');
xlswrite([filepath{1},name,' ',dt,'.xlsx'], headers, 'Sheet1','A1');
xlswrite([filepath{1},name,' ',dt,'.xlsx'], odata(leftInd:rightInd,:), 'Sheet1','A2');
msgbox('finished.');

%% Initialize picking for curve aligning
manualOffsetIndices = cell(size(lineHandles,1), size(lineHandles,2));
arrivalTimes = cell(size(lineHandles,1), size(lineHandles,2));
baseVoltages = cell(size(lineHandles,1), size(lineHandles,2));
baseVoltageIndices = cell(size(lineHandles,1), size(lineHandles,2));

%% Picking for curve aligning
try
    dcm_obj = datacursormode(gcf);
    c_info = getCursorInfo(dcm_obj);   
    if mod(length(c_info),2) ~= 0
        error('Repick!');
    end
    c_info = fliplr(c_info);
    for i = 1:length(c_info)
        for ii = 1:size(lineHandles,1)
            for jj = 1:size(lineHandles,2)
                if isequal(c_info(i).Target, lineHandles{ii,jj})
                    xValue = c_info(i).Position(1);
                    manualOffsetIndices{ii,jj}(2-mod(i,2)) = ...
                        find(lineHandles{ii,jj}.XData >= xValue, 1);
                end
            end
        end
    end
catch ME
    error('Please repick.');
end

msgbox(['Please re-plot now. FYI, peak was at ',...
    num2str(mean2(cell2mat(arrivalTimes))),'.']);

%% Export picks of peaks
dt = datestr(now,'mmmm_dd_yyyy_HH_MM_SS');
[~,name,~] = fileparts(filename{1});
save([filepath{1},'picks_',name,' ',dt,'.mat'],'manualOffsetIndices');
msgbox('finished.');

%% Calculate arrival times, velocity field
gaugeDistance = 0.0314;    % in meters

averageArrivalTimes = zeros(1,10);
averageBaseVoltageIndices = zeros(1,10);
cnt = 0;
for i = 1:size(arrivalTimes,1)
    for j = 1:size(arrivalTimes,2)/3
        gaugeCnt = 0;
        gaugeSum = 0;
        gaugeIndSum = 0;
        for k = 1:3
            if isempty(arrivalTimes{i,j+k-1})
                continue
            else
                gaugeSum = gaugeSum + arrivalTimes{i,j+k-1};
                gaugeIndSum = gaugeIndSum + baseVoltageIndices{i,j+k-1};
                gaugeCnt = gaugeCnt + 1;
            end
        end
        cnt = cnt + 1;
        if gaugeCnt > 0
            averageArrivalTimes(cnt) = gaugeSum/gaugeCnt;
            averageBaseVoltageIndices(cnt) = round(gaugeIndSum/gaugeCnt);
        else
            error('You did not pick at least one channel per gauge.');
        end
    end
end

deltaTimes = averageArrivalTimes - circshift(averageArrivalTimes,[0 -1]);
velocityField = gaugeDistance./deltaTimes;

format longE
disp('averageArrivalTimes');
display(averageArrivalTimes');
disp('velocityField');
display(velocityField');

msgbox(['Done. Peak is at ~',num2str(mean(averageArrivalTimes)),...
    '. averageArrivalTimes contains times of peaks. velocityField ',...
    'contains tentatively calculated velocities. See Command Window.']);

%% Export as strain vs distance

GF = 155;
exportOrPlot = 'export';
smoothSpan = 10;
useAverageOfPicks = 1;
usePickedBaseVoltages = 1;
velocityField = ones(1,10);   % for plotting vs time
% velocityField = [];    % custom velocities

xRange = xlim;
avg0 = mean(cell2mat(reshape(arrivalTimes,1,numel(arrivalTimes))));
if ~isnan(avg0)
    xRange = xRange + avg0*plotWithPeaksAligned;
end
dt = datestr(now,'mmmm_dd_yyyy_HH_MM_SS');
[~,name,~] = fileparts(filename{1});
if isequal(exportOrPlot,'plot')
    figure('Name',['strain_',name]);
end
for i = 1:length(filename)
    odata = test{i}(:,1:15);
    leftInd = find(timecell{i} > xRange(1),1,'first');
    rightInd = find(timecell{i} < xRange(2),1,'last');
    if rightInd - leftInd > 500000
        msgbox(['too much data to save: ',...
            num2str(rightInd - leftInd),' lines']);
        return
    end
    tempInd = find(cardToShow == cardSN(i),1);
    % calculate distances
    if useAverageOfPicks
        tempV = averageArrivalTimes(5*(i-1)+1:5*i);
        arrivalTime = reshape([tempV;tempV;tempV],1,15);
    else
        arrivalTime = zeros(1,15);
        for j = 1:15
            if ~isempty(arrivalTimes{tempInd,j})
                arrivalTime(j) = arrivalTimes{tempInd,j};
            end
        end
    end
    velocities = reshape([velocityField(5*(i-1)+1:5*i);
        velocityField(5*(i-1)+1:5*i);
        velocityField(5*(i-1)+1:5*i)], 1, 15);
    distanceData = (timecell{i}(leftInd:rightInd,:)*ones(1,15)).*...
        (ones(rightInd-leftInd+1,1)*velocities);
    distanceAtPeaks = velocities.*arrivalTime;
    distanceData = distanceData - ones(rightInd-leftInd+1,1)*distanceAtPeaks;
    % calculate strains
    if useAverageOfPicks
        tempI = averageBaseVoltageIndices(5*(i-1)+1:5*i);
        baseVoltageInd = reshape([tempI;tempI;tempI],1,15);
    end
    if usePickedBaseVoltages
        baseVoltage = zeros(1,15);
        for j = 1:15
            if useAverageOfPicks
                baseVoltage(j) = test{i}(baseVoltageInd(j),j);
            elseif ~isempty(baseVoltages{tempInd,j})
                baseVoltage(j) = baseVoltages{tempInd,j};
            end
        end
    else 
        baseVoltage = defaultBaseVoltages(i,:);
    end
    strainData123 = (odata(leftInd:rightInd,:)./...
        (ones(1-leftInd+rightInd,1)*baseVoltage)-1)./GF;
    strainDataXYZ = reshape(strainData123,5*(rightInd-leftInd+1),3);
    strainDataXYZ(:,1) = 0.5*(strainDataXYZ(:,3) - strainDataXYZ(:,1));
    strainDataXYZ(:,3) = strainDataXYZ(:,2) - 2*strainDataXYZ(:,1);
    strainDataXYZ = reshape(strainDataXYZ, rightInd-leftInd+1, 15);
    % arrange output data format
    toOutput = [];
    headers = {};
    for j = 1:5
        if useAverageOfPicks
            toOutput = [toOutput distanceData(:,(j-1)*3+1)];
            headers = [headers {['d',chNames{i,(j-1)*3+1}]}];
        end
        for k = 1:3
            if ~useAverageOfPicks
                toOutput = [toOutput distanceData(:,(j-1)*3+k)];
                headers = [headers {['d',chNames{i,(j-1)*3+k}]}];
            end
            toOutput = [toOutput strainData123(:,(j-1)*3+k)];
            headers = [headers {[chNames{i,(j-1)*3+k}]}];
        end
        toOutput = [toOutput strainDataXYZ(:,(j-1)*3+1:(j-1)*3+3)];
        headers = [headers {'XY','YY','XX'}];
    end
    if isequal(exportOrPlot,'export')
        % output
        xlswrite([filepath{1},'strain_',name,' ',dt,'.xlsx'],...
            headers, [filename{i}(1:10),num2str(i)], 'A1');
        xlswrite([filepath{1},'strain_',name,' ',dt,'.xlsx'],...
            toOutput, [filename{i}(1:10),num2str(i)], 'A2');
    elseif isequal(exportOrPlot,'plot')
        hold on
        for i = 1:15
            plot(distanceData(:,i), smooth(strainData123(:,i),smoothSpan));
        end
        hold off
    end
end
if isequal(exportOrPlot,'export')
    msgbox('finished.');
end

%% Truncate original data and save to file
preCut = 1;    % in seconds
postCut = 1;    % in seconds

% execution begins here
dt = datestr(now,'mmmm_dd_yyyy_HH_MM_SS');
for i = 1:length(cardSN)
    tempData  = raw_test{i}(freq*preCut:end-freq*postCut,:);
    tempData = reshape(tempData,1,length(tempData)*17);
    [~,name,~] = fileparts(filename{i});
    dlmwrite([filepath{i},'chopped_',name,' ',dt,'.txt'],tempData);
end
msgbox('finished.');