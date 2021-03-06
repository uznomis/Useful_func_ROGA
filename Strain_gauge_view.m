% this program loads data produced by ROGA Terminal developed by James
% Jeffers, display the data according to the user request
%% Initializing
% paramters to change before continuing
useLowFrequencyData = 1;
lowFreqDataChannelNames = {'Time','Axial Load','Torque','Eddy2','Encode Volt'};
lowFreqDataChannelFactors = [1e-3, 1, 1, 1, 1];
lowFreqEncoderChannel = 5;

cardSN = [1,2];    % which cards to import
chNames = {'J3','J2','J1','I3','I2','I1','H3','H2','H1',...
    'G3','G2','G1','F3','F2','F1','Encoder','Unused';...
    'E3','E2','E1','D3','D2','D1','C3','C2','C1',...
    'B3','B2','B1','A3','A2','A1','Encoder','Accel'};
freq = 1e6;    % freqeuncy in Hz
downsampleRate = 1;    % remember to also change parameters in Sync cards properly
cardToGetEncoder = 0;    % the number of card whose encoder data is used for velocity and counter calculation; put 0 if you don't want velocity and counter
velocityInterval = 1e6;    % interval to use in calculating velocity from encoder
encoderChannels = [16;16];    % encoder channels for the two cards; should match the ordering in cardSN
encoderAvailable = 1;    % 1 for yes, 0 for no; if no encoder is available, there is no sync between cards
lenPerSector = 1.5e-6;    % encoder sector length in meters
% baseLevelFactors = [20.58,11]; % for PMMA sample
% baseLevelOffsets = [3.464,2.036]; % for PMMA sample
% baseLevelFactors = [20.58 20.58 20.58 20.58 20.58 20.58 20.58 20.58 20.58 20.58 20.58 20.58 20.58 20.58 20.58;
%   11 11 11 11 11 11 11 11 11 11 11 11 11 11 11]; % for PMMA sample gage voltage amplification
% baseLevelOffsets = [3.464 3.464 3.464 3.464 3.464 3.464 3.464 3.464 3.464 3.464 3.464 3.464 3.464 3.464 3.464 ;
%    2.036 2.036 2.036 2.036 2.036 2.036 2.036 2.036 2.036 2.036 2.036 2.036 2.036 2.036 2.036]; % for PMMA sample, voltage shift for gages

% baseLevelFactors = [21.06,11]; % for SWG sample
% baseLevelOffsets = [3.623,2.048]; % for SWG sample
% baseLevelFactors = [11.30,15.60]; % for SWG sample, MEAN VALUES FEB. 9 2018
% baseLevelOffsets = [1.68,2.6]; % for SWG sample, MEAN VALUES FEB. 9 2018
baseLevelFactors = [14.94 11.09 11.10 11.19 11.04 11. 10.97 11.02 11.0 10.97 10.99 11.01 10.99 10.97 11.05;
    14.37 10.73 10.7 10.1 10.79 10.77 10.78 19.44 19.45 19.46 19.47 19.39 19.48 19.44 19.43]; % for SWG sample, Measured VALUES FEB. 9 2018
baseLevelOffsets = [1.68 1.68 1.68 1.68 1.68 1.68 1.68 1.68 1.68 1.68 1.68 1.68 1.68 1.68 1.68;
    2.04 2.04 2.04 2.04 2.03 2.04 1.8 3.29 3.24 3.25 3.30 3.26 3.30 3.28 3.29]; % for SWG sample, Calculated VALUES FEB. 9 2018

% Note: order of individual gages data is J3, J2, J1, I3 ,,,,,,,,,,,,,A3, A2, A1
% order of set of gages is J, I,  H, ,,,,A
% [22.2 22.1 22.2 21.9 21.8 21.8 21.2 21.9 21.8 21.8 10.3 21.8 21.7 21.7 21.9
% 11.0 10.8 11.0 11.1 11.1 11.0 11.0 11.1 11.0 11.1 10.4 11.0 11.1 11.1
% 11.0] % amplification of SWG gages
% defaultBaseVoltages = ones(2,15);

% defaultBaseVoltages = [112.6 123 119 124.2 148.8 129 137.3 134.2 146.3 142.6 162.8 183.1 169.3 171.1 192.1;
%    147.1 165.7 148.2 177.4 185.4 183.6 149.7 169 149.8 175.4 179.2 177.6 162.7 188.7 169];    % PMMMA data in mV no load voltage
% defaultBaseVoltages = [111.1 121.8 120.1 123.8 123.5 125.3 139.2 142.6 147.6 150.8 158.8 160.9 162.8 168.9 174.3;
%     154.2 198.1 151.9 143.9 143 143.6 154 146.2 147.8 147.1 152.4 146.5 140.4 140.5 158.7];    % SWG data in mV
% defaultBaseVoltages = [116.4	127.2	125.7	129.0	127.6	129.1	139.7	148.1	153.1	156.0	117.2	166.7	168.6	174.4	182.1;
%    153.5	189.9	151.3	144.0	143.5	143.2	153.9	146.6	147.5	147.2	143.6	146.3	141.1	141.3	157.2]; % SWG data in mV in run 6151
defaultBaseVoltages = [112.4 115.1 115.4 119.9 119.9 122.7 137.6 141.4 146.9 150 158.6 160.8 161.3 164.5 172.3;
    143.5 189.4 152.5 152 140.5 141.7 156 144.8 138.5 139.4 145.6 141.4 145.6 143.5 144.6]; % SWG data in mV of Feb. 7, 2018
% anglesRelativeToFault = [55.8000 47.9000 43.4000 50.8000 41.6000;
%    54.5000 45.3000 47.2000 46.4000 40.6000]; % data for PMMA
distanceFromFault = [2.1 2.4 2.1 2.1 2.3 2.15 2.2 2.7 2.3 2]; % distance from fault in mm for PMMA
% distanceBetweenGauges = 36.05 32.9 33.15 30.45 32.64 30.85 30.85 31.1
%       31.1 31.8 % distance between gouges in mm for PMMA for the following order: J to A, I to J, .... A to B
anglesRelativeToFault = [54.0000    51.0000 44.0000 51.0000 46.0000;
    49.0000    42.0000 54.0000 48.0000 53.0000]; % data for SWG

%% Importing
filename = cell(1,length(cardSN));
filepath = cell(1,length(cardSN));
for i = 1:length(cardSN)
    if ~exist('temppath','var')
        temppath = '*.txt';
    end
    [filename{i}, filepath{i}] = uigetfile...
        (temppath,...
        ['Please select ACQ100', num2str(cardSN(i))]);
    
    if filename{i} == 0
        return
    else
        temppath = [filepath{i},'*.txt'];
    end
end

% load low frequency file
if useLowFrequencyData
    [filenameLowFreq, filepathLowFreq] = uigetfile...
        ('*.*', 'Please select Low Frequency File');
    if filenameLowFreq ~= 0
        lowFreqRawData = dlmread([filepathLowFreq, filenameLowFreq],',',1,0);
    end
    lowFreq = 1e3/lowFreqRawData(2,1); % get frequency from time(ms)
    lowFreqData = lowFreqRawData.*lowFreqDataChannelFactors;
    lowFreqData(end,1) = lowFreqData(end-1,1) + 1/lowFreq;
end

test = cell(1,length(cardSN));
raw_test = test;
for i = 1:length(cardSN)
    test{i} = dlmread([filepath{i}, filename{i}]);
    test{i} = reshape(test{i}, [length(test{i})/17 17]);
    raw_test{i} = test{i};
    test{i} = downsample(test{i},downsampleRate);
    r_dim = [length(test{1}) 1];
    test{i}(:,1:15) = 1e3*(test{i}(:,1:15)+...
        repmat(baseLevelOffsets(i,:),r_dim))...
        ./ repmat(baseLevelFactors(i,:),r_dim);
end
raw_freq = freq;
freq = freq / downsampleRate;

%% Velocity and counter calculation
% calculate s and v from encoder; calling another function
if cardToGetEncoder == 0
    s = zeros(1,length(test{1}));
    v = s;
    encoderInd = 1;
else
    encoderInd = find(cardSN == cardToGetEncoder, 1);
    [s,v] = getCounterFromEncoder(test{encoderInd}(:,...
        encoderChannels(encoderInd)),freq,lenPerSector, velocityInterval);
end

%% Sync cards
% parameters regarding synchronization; no need to change unless sync is
% failed; if failed, observe the encoder pattern: if the encoder is sparse
% then increase encoderShpacingThrsh and decrease otherwise, and if there is
% a warning saying matrix exceeds dimension then decrease
% encoderXcorrWindow.
encoderShift = 0.5;    % encoder is usually 0/4.6 volts, so pick in between
encoderSpacingThrsh = 1e5; % for high frequency run
encoderXcorrWindow = 1e5; % for high frequency run
% encoderSpacingThrsh = 3.5e4; % encoderXcorrWindow = 7.05e5; % for run 7118
% encoderSpacingThrsh = 1e3; % for diluted runs in which the frequency is reduced
% encoderXcorrWindow = 1e5; % for diluted runs in which the frequency is reduced
useFullTrace = 1;

% execution begins here
if ~isempty(filename)
    % calc timeshift of cards
    trace = cell(1,length(test));
    encoderjumpavg = zeros(1,length(test) - 1);
    encoderjumps = zeros(1,length(test));
    if encoderAvailable == 1
        for ii = 1:length(test)
            data1 = test{ii};
            encoderCh1 = encoderChannels(ii);
            if useFullTrace
                trace{ii} = data1(:,encoderCh1);
                continue
            end
            encodercurrent = 1;
            encoderbefore = 0;
            encoderjump = 1;
            spacing = 0;
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
            encoderjumps(ii) = encoderjump;
            trace{ii} = data1(encoderjump:encoderjump + encoderXcorrWindow,encoderCh1);
        end
        
        for i=2:length(test)
            [r,lag] = xcorr(trace{1},trace{i});
            [~,I] = max(abs(r));
            encoderjumpavg(i-1) = round(-lag(I));
        end
    end
    % create time col vector
    timecell = cell(1,length(test));
    timecell{1} = 1/freq * (1:length(test{1}))';
    
    for i=2:length(test)
        timecell{i} = ...
            timecell{1}-(...
            encoderjumpavg(i-1)+...
            (encoderjumps(i)-encoderjumps(1))*(1-useFullTrace)...
            )/freq; % apply time shift
    end
    
end

%% Sync with Low Frequency Data
if useLowFrequencyData
    dilutedEncoderHighFreq = downsample(trace{1}, freq/lowFreq);
    minLength = min(length(dilutedEncoderHighFreq), length(lowFreqData));
    [r,lag] = xcorr(dilutedEncoderHighFreq(1:minLength), lowFreqData(1:minLength,lowFreqEncoderChannel));
    [~,I] = max(abs(r));
    offsetWithLowFreq = round(-lag(I));
    timeLowFreq = lowFreqData(:,1) - offsetWithLowFreq/lowFreq;
end

%% Plotting
% can display data median-filtered, smoothened and demeaned, as well as
% scaled. After the above sections are run, just simply run the section to
% plot different versions of data plots catered to your needs.

% parameters to change before executing the following code
cardToShow = [1,2];    % cards to display%
% chSN = [16,17,1;16,17,1]; % encoder data
% chSN = [6,5,4,3,2,1];
% chSN = [1:15;1:15];    % for velocity field picking
% chSN = [1:16:3;1:16:3];
% chSN = [3,6,9;3,6,9];
chSN = [16,17,13,10,7,4,1;16,17,13,10,7,4,1]; % first gage (shear)
% chSN = [13,10,7,4,1,16;13,10,7,4,1,16]; % third gage (shear)
% chSN = [16;16];    % channels to display on each card
% chSN = [16,14,11,8,5,2;16,14,11,8,5,2];    % normal load channels to display on each card
% chSN = [8,5,2;8,5,2];
% chSN = [2,5,8,11,14;2,5,8,11,14];    % channels to display on each card
% chSN = [1,3,4,6,7,9,10,12,13,15;1,3,4,6,7,9,10,12,13,15];    % channels to display on each card
lowFreqChannels = [2,3,4,5];
smoothSpan = [50,50];    % smoothening window; make it 1 to diable smoothening; note there is no smoothening for AE data
medFilterOn = [1,1];    % choosing 1 makes median filter on; median filter gets rid of spikes
cardOffset = 0.0;    % offset between cards when plotting
chOffset = 0.05;    % offset between channels within each card when plotting
cardAmp = [1,1];
encoderVelDis = 2;    % velocity & distance from encoder; 0 for not plotting; 1 for only velocity; 2 for v & d
sSlope = 2e2;    % slope of counter for scaling when plotting
vSlope = 2e2;    % slope of velocity for scaling when plotting
encoderSlope = [0.01,0.1];
accelCard = 2;    % only support one card
accelCh = 16;    % only support one channel
accelAmp = 1;
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

if useLowFrequencyData
    for i = 1:length(lowFreqChannels)
        plot(timeLowFreq, lowFreqData(:,lowFreqChannels(i)));
    end
end

strings = cell(length(cardToShow)*length(chSN)...
    + length(lowFreqChannels)*useLowFrequencyData + encoderVelDis,1);

for j = 1:length(cardToShow)
    tempInd = find(cardSN == cardToShow(j),1);
    for i = 1:length(chSN)
        strings{(j - 1)*length(chSN) + i} ...
            = sprintf([num2str(cardToShow(j)),' Ch ',num2str(chSN(j,i)),...
            ' ',chNames{tempInd,chSN(j,i)}]);
    end
end

if useLowFrequencyData
    for i = 1:length(lowFreqChannels)
        strings{length(cardToShow)*length(chSN) + i} = sprintf(['Low ',lowFreqDataChannelNames{lowFreqChannels(i)}]);
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
    % odata = [timecell{i} data1]; % line added by ZR feb 8 2018
    
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
accelInd = find(cardSN == accelCard, 1);
odata = [timecell{encoderInd} v' s'];
odataAccel = [timecell{accelInd} test{accelInd}(:,accelCh)];
headers = {'Time','Velocity','Distance','Time Accel','Accel'};
leftInd = find(odata(:,1) > xRange(1),1,'first');
rightInd = find(odata(:,1) < xRange(2),1,'last');
leftIndAccel = find(odataAccel(:,1) > xRange(1),1,'first');
rightIndAccel = find(odataAccel(:,1) < xRange(2),1,'last');
xlswrite([filepath{1},name,' ',dt,'.xlsx'], headers, 'Sheet1','A1');
xlswrite([filepath{1},name,' ',dt,'.xlsx'],odata(leftInd:rightInd,:), 'Sheet1','A2');
xlswrite([filepath{1},name,' ',dt,'.xlsx'],...
    odataAccel(leftIndAccel:rightIndAccel,:), 'Sheet1','D2');
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
useSimpleXYZConversion = 0;    % make 0 to use accurate XYZ calculation
exportOrPlot = 'plot'; % enter 'export' for saving data
smoothSpan = 50;    % used for BOTH 'plot' and 'export'
detrendLines123 = 1; % setting for clear plotting
detrendLinesXYZ = 1; % setting for clear plotting
%detrendLines123 = 0; % setting for true values relative to zero for export
%detrendLinesXYZ = 0; % setting for true values relative to zero for export
% outputFormat = '123';
outputFormat = 'XYZ';
% cardOffset = 1e-3;
% chOffset = 5e-4;
cardOffset = 2e-5;
chOffset = 3e-3;
chOffset = 0;
% cardOffset = 0;
chAmp = 50; % only for display, not for export
color = {'r','k','b'};
XYZPicksReady = 0;
componentToUse = 1; % 1:XY, 2:YY, 3:XX; this is the component to use to do time/distance shift
componentToReplaceAsRatio = 0;    % 0:no effect, 1:XY, 2:YY, 3:XX; the specified component will be replaced with XY/YY ratio

% usually you don't need to change things below
% !!below used only for strain123!!
picksAvailable = 0;    % default 0; make it 1 ONLY when picking actual voltages
useAverageOfPicks = 1;    % 1 for only 1 channel picked; if picksAvailable = 0, then this has no effect
usePickedBaseVoltages = 0;    % if picksAvailable = 0, then this has no effect
velocityField = ones(1,10);   % for plotting vs time
% velocityField = [408 -199 -403 628 -21 33 78 234 897 408];    % custom velocities
% velocityField = [200 100 50 100 250 100 500 50 30 100];    % custom velocities
% !!above used only for strain123!!

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
    % calculate distances
    arrivalTime = zeros(1,15);
    if picksAvailable
        if useAverageOfPicks
            tempV = averageArrivalTimes(5*(i-1)+1:5*i);
            arrivalTime = reshape([tempV;tempV;tempV],1,15);
        else
            for j = 1:15
                if ~isempty(arrivalTimes{i,j})
                    arrivalTime(j) = arrivalTimes{i,j};
                end
            end
        end
    elseif XYZPicksReady && exist('XYZPicks','var')
        tempV = reshape(XYZPicks(i,:),3,5);
        tempV = tempV(componentToUse,:) - customXYZshifts((i-1)*5+1:i*5);
        arrivalTime = reshape([tempV;tempV;tempV],1,15);
        velocityField = velocityXYZ;
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
            elseif ~isempty(baseVoltages{i,j})
                baseVoltage(j) = baseVoltages{i,j};
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
        if useSimpleXYZConversion
            strainDataXYZ = [strainDataXYZ; strainData123(:,(j-1)*3+1:j*3)];
        else
            theta = anglesRelativeToFault(i,j);
            strainDataXYZ = [strainDataXYZ;
                calculateStrainXYZ(strainData123(:,(j-1)*3+1:j*3),theta)];
        end
    end
    if useSimpleXYZConversion
        strainDataXYZ(:,1) = 0.5*(strainDataXYZ(:,3) - strainDataXYZ(:,1));
        strainDataXYZ(:,3) = strainDataXYZ(:,2) - 2*strainDataXYZ(:,1);
    end
    if componentToReplaceAsRatio ~= 0
        strainDataXYZ(:,componentToReplaceAsRatio) = ...
            strainDataXYZ(:,1)./strainDataXYZ(:,2);
    end
    strainDataXYZ_ = strainDataXYZ;
    strainDataXYZ = [];
    distanceFromFaultBar = mean(distanceFromFault);
    for j = 1:5
        strainDataXYZ = [strainDataXYZ...
            strainDataXYZ_((rightInd - leftInd +1)*(j-1)+1:(rightInd - leftInd +1)*j,:)...
            .*(distanceFromFault((i-1)*5+j)/distanceFromFaultBar)^(-1/2)];
    end
    strainDataXYZ = strainDataXYZ - detrendLinesXYZ*...
        ones(1-leftInd+rightInd,1)*strainDataXYZ(1,:);
    % arrange output data format
    if isequal(outputFormat,'123')
        strainData = strainData123;
    else
        strainData = strainDataXYZ;
    end
    for j = 1:15
        strainData(:,j) = smooth(strainData(:,j),smoothSpan);
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
                strainData(:,j)*chAmp,color{mod(j-1,3)+1});
        end
        hold off
    end
end
% calculate mean of 10 gauges
l1 = length(strainDatas{1});
l2 = length(strainDatas{2});
if l1 > l2
    strainDatas{1} = strainDatas{1}(1:l2,:);
elseif l1 < l2
    strainDatas{2} = strainDatas{2}(1:l1,:);
end
strainDatasMat = cell2mat(strainDatas);
meanOfGauges = zeros(length(strainDatasMat),3);
for i = 1:3
    meanOfGauges(:,i) = mean(strainDatasMat(:,i:3:end),2);
end
% labels
if isequal(exportOrPlot,'export')
    xlswrite([filepath{1},'strain_',outputFormat,'_',name,' ',dt,'.xlsx'],...
        {'avg1','avg2','avg3'}, [filename{1}(1:10),'mean'], 'A1');
    xlswrite([filepath{1},'strain_',outputFormat,'_',name,' ',dt,'.xlsx'],...
        meanOfGauges, [filename{1}(1:10),'mean'], 'A2');
    msgbox('finished.');
else
    hold on
    for i = 1:3
        plot(distanceData(:,i), i*cardOffset+3*chOffset+...
            meanOfGauges(:,i)*chAmp,'color',color{mod(i-1,3)+1},'linewidth',2.5);
    end
    hold off
    if isequal(outputFormat,'123')
        legend([chNames(1,1:15) chNames(2,1:15) {'avg1', 'avg2', 'avg3'}]);
    else
        chNamesXYZ = {};
        chTable = {'XX','YY','XY'};
        for i = 1:30
            chName = chNames{i};
            chNamesXYZ = [chNamesXYZ {[chName(1),chTable{str2double(chName(2))}]}];
        end
        chNamesXYZ = reshape(chNamesXYZ, 2, 15);
        legend([chNamesXYZ(1,:) chNamesXYZ(2,:) {'avg1', 'avg2', 'avg3'}]);
    end
end

%% Hide lines for XYZ
hideLines = [0 1 1];    % XY, YY, XX respectively
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
clear('figPicked');

%% Fix zoom
fixedZoom = {xlim,ylim};

%% Auto zoom
xlim auto
ylim auto
clear('fixedZoom');

%% Picking/plotting for XYZ strains

% customXYZshifts = zeros(1,10);
customXYZshifts = [0 0 0 0 0 0 0 0 0 0];    % custom time shift
velocityXYZ = ones(1,10);    % custom velocity
% velocityXYZ   I-J,  H-I, G-H, F-G,   E-F,  D-E,  C-D,  B-C,  A-B, J-A
%velocityXYZ = [5 5 5 5 5 5 3 3 5 5];    % custom velocities
% velocityXYZ = [700 700 2545 2040 1000 907 1928 2591 700 2253];    % custom velocities for secondary event PMMA 7039
% velocityXYZ = [1316 1200 3045 2040 1000 907 1928 2591 1110 2253];
% velocityXYZ = [2350 2762 2175 3626 2057 -2204 -3456 -2221 -1000 -1803]; % velocities for 7282


componentOffset = 2e-4;    % offset between groups of lines in XY, YY, and XX
componentToDisplay = 1;  % 1:XY, 2:YY, 3:XX for display in command line only
alignOption = 'right'; % 'left'; 'right'; 'none' ('none' means XX YY XY are separated by componentOffset)
widthToAlign = 100; % how long of a segment (num of points) on the left (right) is used to aligning
% ampArray = [0.5 2 1 1 1 1 1 1 1 0.5]; For 7118
autoAmp = 1; % set to 1 to automatically calculates amplification instead of using ampArray
% ampArray = [1 1.3 1 1 1 1 .7 .7 1 .7]; % for PMMA 7039
plotOpt = '.'; % for plotting in dots
plotOpt = '-'; % for plotting in lines
smoothSpan = 1;

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
if exist('figPicked','var') && isvalid(figPicked)
    figure(figPicked);
else
    figPicked = figure('Name',['strain_picked_',name]);
end
clf(figPicked);
chNamesXYZ = {};
chTable = {'XX','YY','XY'};
alignOffsets = zeros(1,10);
hold on
for i = 1:2
    leftInd = find(timecell{i} > xRange(1),1,'first');
    rightInd = find(timecell{i} < xRange(2),1,'last');
    for j = 1:15
        if XYZPicks(i,j) ~= 0
            % magic happens in the next sentence, read carefully before
            % changing any of it
            if autoAmp
                ampArray(ceil(((i-1)*15+j)/3)) = 1/abs(mean(strainDatas{i}(end-widthToAlign+1:end,j))...
                    -mean(strainDatas{i}(1:widthToAlign,j)));
            end
            if isequal(alignOption, 'left')
                alignOffsets(ceil(((i-1)*15+j)/3)) = mean(strainDatas{i}(end-widthToAlign+1:end,j));
            elseif isequal(alignOption, 'right')
                alignOffsets(ceil(((i-1)*15+j)/3)) = mean(strainDatas{i}(1:widthToAlign,j));
            elseif isequal(alignOption, 'none')
                alignOffsets(ceil(((i-1)*15+j)/3)) = -mod(j-1,3)*componentOffset;
            end
            plot(0-(timecell{i}(leftInd:rightInd)-XYZPicks(i,j)...
                -customXYZshifts(ceil(((i-1)*15+j)/3)))...
                *velocityXYZ(ceil(((i-1)*15+j)/3)),...
                (smooth(strainDatas{i}(:,j),smoothSpan) ...
                - alignOffsets(ceil(((i-1)*15+j)/3)))...
                *ampArray(ceil(((i-1)*15+j)/3)), plotOpt);
            chName = chNames{i,j};
            chNamesXYZ = [chNamesXYZ {[chName(1),chTable{str2double(chName(2))}]}];
        end
    end
end
legend(chNamesXYZ);
hold off
if exist('fixedZoom','var')
    xlim(fixedZoom{1});
    ylim(fixedZoom{2});
end
format LongE
disp('Picks are:');
tempPicks = [XYZPicks(1,:) XYZPicks(2,:)]';
tempPicks = reshape(tempPicks,3,10);
disp(tempPicks(componentToDisplay,:)');

%% Truncate original data and save to file
preCut = 1;    % in seconds
postCut = 1;    % in seconds

% execution begins here
dt = datestr(now,'mmmm_dd_yyyy_HH_MM_SS');
for i = 1:length(cardSN)
    tempData  = raw_test{i}(raw_freq*preCut+1:end-raw_freq*postCut,:);
    tempData = reshape(tempData,1,length(tempData)*17);
    [~,name,~] = fileparts(filename{i});
    dlmwrite([filepath{i},'chopped_',name,' ',dt,'.txt'],tempData);
end
msgbox('finished.');