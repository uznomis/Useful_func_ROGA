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
baseLevelFactors = [20.58,11];
baseLevelOffsets = [3.464,2.036];
% defaultBaseVoltages = ones(2,15);
defaultBaseVoltages = [112.6 123 119 124.2 148.8 129 137.3 134.2 146.3 142.6 162.8 183.1 169.3 171.1 192.1;
    147.1 165.7 148.2 177.4 185.4 183.6 149.7 169 149.8 175.4 179.2 177.6 162.7 188.7 169];    % in mV

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
    test{i}(:,1:15) = 1e3*(test{i}(:,1:15)+baseLevelOffsets(i)) / baseLevelFactors(i);
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
chSN = [1:15;1:15];    % for velocity field picking
% chSN = [1:16:3;1:16:3];
% chSN = [3,6,9;3,6,9];
% chSN = [16,13,10,7,4,1;16,13,10,7,4,1]; % first gage (shear)
% chSN = [13,10,7,4,1,16;13,10,7,4,1,16]; % third gage (shear)
% chSN = [16;16];    % channels to display on each card
% chSN = [14,11,8,5,2;14,11,8,5,2];    % normal load channels to display on each card
% chSN = [8,5,2;8,5,2];
% chSN = [2,5,8,11,14;2,5,8,11,14];    % channels to display on each card
% chSN = [1,3,4,6,7,9,10,12,13,15;1,3,4,6,7,9,10,12,13,15];    % channels to display on each card
smoothSpan = [5,5];    % smoothening window; make it 1 to diable smoothening; note there is no smoothening for AE data
medFilterOn = [1,1];    % choosing 1 makes median filter on; median filter gets rid of spikes
cardOffset = 20;    % offset between cards when plotting
chOffset = 2;    % offset between channels within each card when plotting
cardAmp = [10,10];
encoderVelDis = 0;    % velocity & distance from encoder; 0 for not plotting; 1 for only velocity; 2 for v & d
sSlope = 2e2;    % slope of counter for scaling when plotting
vSlope = 2e2;    % slope of velocity for scaling when plotting
encoderSlope = [0.01,0.1];
accelCard = 2;    % only support one card
accelCh = 16;    % only support one channel
accelAmp = 50;
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
if 1
    return
end

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
if exist('arrivalTimes','var')
    avg0 = mean(cell2mat(reshape(arrivalTimes,1,numel(arrivalTimes))));
    if ~isnan(avg0)
        xRange = xRange + avg0*plotWithPeaksAligned;
    end
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
    msgbox(['Please re-plot now. FYI, last peak you picked was at ',...
        num2str(xValue),' s.']);
catch ME
    error('Please repick.');
end

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
            if isempty(arrivalTimes{i,(j-1)*3+k})
                continue
            else
                gaugeSum = gaugeSum + arrivalTimes{i,(j-1)*3+k};
                gaugeIndSum = gaugeIndSum + baseVoltageIndices{i,(j-1)*3+k};
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

%% Export/plot as strain vs distance/time

GF = 155;
exportOrPlot = 'plot';
smoothSpan = 5;    % used only for 'plot' not 'export'
% !!below used only for strain123!!
picksAvailable = 0;
useAverageOfPicks = 1;    % 1 for only 1 channel picked; if picksAvailable = 0, then this has no effect
usePickedBaseVoltages = 0;    % if picksAvailable = 0, then this has no effect
% !!above used only for strain123!!
detrendLines123 = 1;
detrendLinesXYZ = 1;
% outputFormat = '123';
outputFormat = 'XYZ';
cardOffset = 1e-3;
chOffset = 1e-4;
color = {'r','k','b'};
velocityField = ones(1,10);   % for plotting vs time
% velocityField = [408 -199 -403 628 -21 33 78 234 897 408];    % custom velocities
% velocityField = [200 100 50 100 250 100 500 50 30 100];    % custom velocities
XYZPicksReady = 1;

xRange = xlim;
if exist('arrivalTimes','var')
    avg0 = mean(cell2mat(reshape(arrivalTimes,1,numel(arrivalTimes))));
    if ~isnan(avg0)
        xRange = xRange + avg0*plotWithPeaksAligned;
    end
end
dt = datestr(now,'mmmm_dd_yyyy_HH_MM_SS');
if ~exist('lineHandlesXYZ','var')
    lineHandlesXYZ = cell(2,15);
end
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
    arrivalTime = zeros(1,15);
    if picksAvailable
        if useAverageOfPicks
            tempV = averageArrivalTimes(5*(i-1)+1:5*i);
            arrivalTime = reshape([tempV;tempV;tempV],1,15);
        else            
            for j = 1:15
                if ~isempty(arrivalTimes{tempInd,j})
                    arrivalTime(j) = arrivalTimes{tempInd,j};
                end
            end
        end
    elseif XYZPicksReady && exist('XYZPicks','var')
        arrivalTime = XYZPicks(i,:) - customXYZshifts(i,:);
    end
    velocities = reshape([velocityField(5*(i-1)+1:5*i);
        velocityField(5*(i-1)+1:5*i);
        velocityField(5*(i-1)+1:5*i)], 1, 15);
    distanceData = (timecell{i}(leftInd:rightInd,:)*ones(1,15)).*...
        (ones(rightInd-leftInd+1,1)*velocities);
    distanceAtPeaks = velocities.*arrivalTime;
    distanceData = distanceData - ones(rightInd-leftInd+1,1)*distanceAtPeaks;
    % calculate strains
    if picksAvailable && usePickedBaseVoltages
        baseVoltage = zeros(1,15);
        if useAverageOfPicks
            tempI = averageBaseVoltageIndices(5*(i-1)+1:5*i);
            baseVoltageInd = reshape([tempI;tempI;tempI],1,15);
        end        
        for j = 1:15
            if useAverageOfPicks
                baseVoltage(j) = odata(baseVoltageInd(j),j);
            elseif ~isempty(baseVoltages{tempInd,j})
                baseVoltage(j) = baseVoltages{tempInd,j};
            end
        end
    else
         baseVoltage = odata(leftInd,:);
    end
    strainData123 = (odata(leftInd:rightInd,:)./...
        (ones(1-leftInd+rightInd,1)*defaultBaseVoltages(i,:))-1)./GF;
    strainData123 = strainData123 - detrendLines123*...
        ones(1-leftInd+rightInd,1)*strainData123(1,:);
    strainDataXYZ = [];
    for j = 1:5
        strainDataXYZ = [strainDataXYZ; strainData123(:,(j-1)*3+1:j*3)];
    end    
    strainDataXYZ(:,1) = 0.5*(strainDataXYZ(:,3) - strainDataXYZ(:,1));
    strainDataXYZ(:,3) = strainDataXYZ(:,2) - 2*strainDataXYZ(:,1);
    strainDataXYZ_ = strainDataXYZ;
    strainDataXYZ = [];
    for j = 1:5
        strainDataXYZ = [strainDataXYZ...
            strainDataXYZ_((rightInd - leftInd +1)*(j-1)+1:(rightInd - leftInd +1)*j,:)];
    end
    strainDataXYZ = strainDataXYZ - detrendLinesXYZ*...
        ones(1-leftInd+rightInd,1)*strainDataXYZ(1,:);
    % arrange output data format
    if isequal(outputFormat,'123')
        strainData = strainData123;
    else
        strainData = strainDataXYZ;
    end
    strainDatas{i} = strainData;
    toOutput = [];
    headers = {};
    for j = 1:5
        if useAverageOfPicks || ~picksAvailable || isequal(outputFormat,'XYZ')
            toOutput = [toOutput distanceData(:,(j-1)*3+1)];
            chName = chNames{i,(j-1)*3+1};
            headers = [headers {['d',chName(1)]}];
        end
        for k = 1:3
            if ~useAverageOfPicks && picksAvailable && isequal(outputFormat,'123')
                toOutput = [toOutput distanceData(:,(j-1)*3+k)];
                headers = [headers {['d',chNames{i,(j-1)*3+k}]}];
            end
            toOutput = [toOutput strainData(:,(j-1)*3+k)];
            if isequal(outputFormat,'123')
                headers = [headers {[chNames{i,(j-1)*3+k}]}];
            end
        end
        if isequal(outputFormat,'XYZ')
            headers = [headers {[chName(1),'XY'],[chName(1),'YY'],[chName(1),'XX']}];
        end
    end
    if isequal(exportOrPlot,'export')
        % output
        xlswrite([filepath{1},'strain_',outputFormat,'_',name,' ',dt,'.xlsx'],...
            headers, [filename{i}(1:10),num2str(i)], 'A1');
        xlswrite([filepath{1},'strain_',outputFormat,'_',name,' ',dt,'.xlsx'],...
            toOutput, [filename{i}(1:10),num2str(i)], 'A2');
    elseif isequal(exportOrPlot,'plot')
        hold on
        for j = 1:15
            lineHandlesXYZ{i,j} = plot(distanceData(:,j), i*cardOffset+j*chOffset+...
                smooth(strainData(:,j),smoothSpan),color{mod(j-1,3)+1});
            if XYZPicksReady
                if XYZPicks(i,j) == 0
                    lineHandlesXYZ{i,j}.set('Visible','off');
                end
            end
        end
        hold off
    end
end
if isequal(exportOrPlot,'export')
    msgbox('finished.');
else
    if isequal(outputFormat,'123')
        legend([chNames(1,1:15) chNames(2,1:15)]);
    else
        chNamesXYZ = {};
        chTable = {'XX','YY','XY'};
        for i = 1:30
            chName = chNames{i};
            chNamesXYZ = [chNamesXYZ {[chName(1),chTable{str2double(chName(2))}]}];
        end
        chNamesXYZ = reshape(chNamesXYZ, 2, 15);
        legend([chNamesXYZ(1,:) chNamesXYZ(2,:)]);
    end
end

%% Hide lines for XYZ
hideLines = [1 0 0];    % XY, YY, XX respectively
for i = 1:2
    for j = 1:15
        if hideLines(mod(j-1,3)+1) == 1
            lineHandlesXYZ{i,j}.set('Visible','off');
        else
            lineHandlesXYZ{i,j}.set('Visible','on');
        end
    end
end
%% Re-initialize XYZ picks
XYZPicks = zeros(2,15);

%% Picking/plotting for XYZ strains

% !!!!!!!Remember to get rid of tempInd!!!!!!!

customXYZshifts = zeros(2,15);
velocityXYZ = ones(2,15);

if ~exist('XYZPicks','var')
    XYZPicks = zeros(2,15);
end
if ~exist('lineHandlesXYZ','var')
    return
end
% Getting picks on plot
try
    dcm_obj = datacursormode(gcf);
    c_info = getCursorInfo(dcm_obj);
    if length(c_info) ~= 0
        c_info = fliplr(c_info);
        for i = 1:length(c_info)
            for ii = 1:size(lineHandlesXYZ,1)
                for jj = 1:size(lineHandlesXYZ,2)
                    if isequal(c_info(i).Target, lineHandlesXYZ{ii,jj})
                        XYZPicks(ii,jj) = c_info(i).Position(1);
                    end
                end
            end
        end
    end
catch ME
    error('Please repick.');
end
if exist('figPicked','var')
    figure(figPicked);
else
    figPicked = figure('Name',['strain_picked_',name]);
end
clf(figPicked);
chNamesXYZ = {};
chTable = {'XX','YY','XY'};
hold on
for i = 1:2
    for j = 1:15
        if XYZPicks(i,j) ~= 0
            plot((timecell{i}(leftInd:rightInd)-XYZPicks(i,j)-customXYZshifts(i,j))...
                *velocityXYZ(i,j),...
                strainDatas{i}(:,j));            
            chName = chNames{i,j};
            chNamesXYZ = [chNamesXYZ {[chName(1),chTable{str2double(chName(2))}]}];
        end
    end
end
legend(chNamesXYZ);
hold off
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