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

%% Importing
filename = cell(1,length(cardSN));
filepath = cell(1,length(cardSN));
for i = 1:length(cardSN)
    [filename{i}, filepath{i}] = uigetfile...
        ('C:\Users\User\Documents\Test files results on fast cards\*.*',...
        ['Please select ACQ100', num2str(cardSN(i))]);
    if filename{i} == 0
        return
    end
end
test = cell(1,length(cardSN));
for i = 1:length(cardSN)
    test{i} = dlmread([filepath{i}, filename{i}]);
    test{i} = reshape(test{i}, [length(test{i})/17 17]);
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
encoderShift = 3;    % encoder is usually 0/4.6 volts, so pick in between
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
cardToShow = [1,2];    % cards to display% 
% chSN = [16,16];
% chSN = [6,5,4,3,2,1];
% chSN = [1:17;1:17];
% chSN = [3,2,1,16;3,2,1,16];
% chSN = [15,12,9,6,3,16;15,12,9,6,3,16]; % first gage (shear)
% chSN = [13,10,7,4,1,16;13,10,7,4,1,16]; % third gage (shear)
% chSN = [16;16];    % channels to display on each card
chSN = [14,11,8,5,2,16;14,11,8,5,2,16];    % normal load channels to display on each card
% chSN = [2,5,8,11,14;2,5,8,11,14];    % channels to display on each card
% chSN = [1,3,4,6,7,9,10,12,13,15;1,3,4,6,7,9,10,12,13,15];    % channels to display on each card
smoothSpan = [10,3];    % smoothening window; make it 1 to diable smoothening; note there is no smoothening for AE data
medFilterOn = [1,1];    % choosing 1 makes median filter on; median filter gets rid of spikes
cardOffset = 0.1;    % offset between cards when plotting
chOffset = 0.005;    % offset between channels within each card when plotting
cardAmp = [4,1];
encoderVelDis = 0;    % velocity & distance from encoder; 0 for not plotting; 1 for only velocity; 2 for v & d
sSlope = 3e2;    % slope of counter for scaling when plotting
vSlope = 3e2;    % slope of velocity for scaling when plotting
encoderSlope = [0.01,0.1];
accelCard = 2;    % only support one card
accelCh = 16;    % only support one channel
accelAmp = 5;

% execution begins here
figureTitle = '';
for i = 1:length(filename)
    figureTitle = [figureTitle,' ',filename{i}];
end
figure('Name',figureTitle);
hold on;

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
        plot(timecell{tempInd},j*cardOffset + i*chOffset + tempCh - mean(tempCh));
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
%% Pick up points
% go to the figure that you want to pick on, pick, and hit Ctrl + Enter.

% Review the parameters below before running this section.
startNewPickSet = 1;    % 1 for starting a new set of points; 0 for appending to the current one

% gather picks
if ~exist('picks','var') || startNewPickSet
    picks = [];
end
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
disp(picks);

%% Export to file
% export the processed data to a text format

% saving gauge data set to excel file
xRange = xlim;
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