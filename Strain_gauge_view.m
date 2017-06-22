% this program loads data produced by ROGA Terminal developed by James
% Jeffers, display the data according to the user request
%% Initializing
% paramters to change before continuing
cardSN = [1, 2];    % which cards to import; numbers should match file names
freq = 5e5;    % freqeuncy in Hz
cardToGetEncoder = 1;    % the number of card whose encoder data is used for velocity and counter calculation
encoderChannels = [16;17];
encoderAvailable = 1;    % 1 for yes, 0 for no; if no encoder is available, there is no sync between cards
lenPerSector = 1.5e-6;    % encoder sector length in meters

%% Importing
for i = 1:2
    [filename{i}, filepath{i}] = uigetfile...
        ('C:\Users\User\Documents\Test files results on fast cards\*.*',...
        ['Please select ACQ100', num2str(i)]);
    if filename{i} == 0
        return
    end
end
for i = 1:2
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
    [s,v] = getCounterFromEncoder(test{cardToGetEncoder}(:,...
        encoderChannels(cardToGetEncoder)),freq,lenPerSector);
end

%% Sync cards
% please provide some important info here
encoderAvailable = 1;
encoderChannel = encoderChannels(1);
encoderChannel2 = encoderChannels(2);
% parameters regarding synchronization; no need to change
encoderShift = 3;
encoderSpacingThrsh = 5e2;
encoderXcorrWindow = 1e5;

% execution begins here
if ~isempty(filename)
    % calc timeshift of cards
    encodercurrent = 1;
    encoderbefore = 0;
    encoderjump = 1;
    spacing = 0;
    data1 = test{1};
    encoderjumpavg = zeros(1,length(filename) - 1);
    if encoderAvailable == 1
        for j = 2:length(data1)
            if (data1(j-1,encoderChannel) - encoderShift)*...
                    (data1(j,encoderChannel) - encoderShift) < 0
                encoderbefore = encodercurrent;
                encodercurrent = j;
                spacing = encodercurrent - encoderbefore;
                if spacing < encoderSpacingThrsh && encoderbefore > 1
                    encoderjump = encodercurrent;
                    break;
                end
            end
        end
        s1 = data1(encoderjump:encoderjump + encoderXcorrWindow,encoderChannel);
        
        for i=2:length(filename)
            s2 = test{i}(encoderjump: encoderjump + encoderXcorrWindow,encoderChannel2);
            [r,lag] = xcorr(s1,s2);
            [~,I] = max(abs(r));
            encoderjumpavg(i-1) = round(-lag(I));
        end
    end
    % create time col vector
    timecell = cell(1,length(filename));
    rawTcell{1} = 1/freq * (1:length(data1))';
    timecell{1} = rawTcell{1};
    
    for i=2:length(filename)
        timecell{i} = ...
            timecell{1}-encoderjumpavg(i-1)/freq; % apply time shift
    end
    
end

%% Plotting
% can display data median-filtered, smoothened and demeaned; there is no
% scaling for strain gage data; there is scaling for AE data; no
% smoothening for AE; smoothening for strain gage data; after importing the
% data file, simply change parameters and run this section with 'Run and
% Advance' button;

% parameters to change before executing the following code
cardToShow = [1,2];    % cards to display
% chSN = [14:17;14:17];
% chSN = [16,15,12,9,6,3;16,15,12,9,6,3]; % first gage (shear) 
chSN = [16,13,10,7,4,1;16,13,10,7,4,1]; % third gage (shear)
% chSN = [16;16];    % channels to display on each card
% chSN = [16,14,11,8,5,2;16,14,11,8,5,2];    % normal load channels to display on each card
% chSN = [2,5,8,11,14;2,5,8,11,14];    % channels to display on each card
% chSN = [1,3,4,6,7,9,10,12,13,15;1,3,4,6,7,9,10,12,13,15];    % channels to display on each card
smoothSpan = 1;    % smoothening window; make it 1 to diable smoothening; note there is no smoothening for AE data
cardOffset = 0.2;    % offset between cards when plotting
chOffset = 0.03;    % offset between channels within each card when plotting
sSlope = 1e3;    % slope of counter for scaling when plotting
vSlope = 1e3;    % slope of velocity for scaling when plotting
encoderSlope = 0.01;
ACQ1001Amp = 10;
accelCh = 16;
accelAmp = 1;
encoderVelDis = 2;    % a constant; don't change
medFilterOn = 1;    % choosing 1 makes median filter on; median filter gets rid of spikes

% execution begins here
figure('Name',[filename{1},' ',filename{2}]);
hold on;

for j = 1:length(cardToShow)
    tempInd = find(cardSN == cardToShow(j),1);
    for i = 1:length(chSN(j,:))
        if medFilterOn
            tempCh = medfilt2(test{tempInd}(:,chSN(j,i)),[3,1]);
        else
            tempCh = test{tempInd}(:,chSN(j,i));
            
        end
        tempCh(1) = tempCh(2); % for getting rid of the anomaly first point
        
        if j == 1 && i ~= encoderChannels(1)
            tempCh = tempCh*ACQ1001Amp;
        end
        
        if j == 2 && chSN(j,i) == accelCh
            tempCh = tempCh*accelAmp;
        end
              
        if ismember(chSN(j,i),encoderChannels) || (j == 2 && chSN(j,i) == accelCh)
            plot(timecell{tempInd},j*...
                cardOffset + i*chOffset + encoderSlope*(tempCh - mean(tempCh)));
        
        else
            plot(timecell{tempInd},j*...
                cardOffset + i*chOffset + smooth(tempCh,smoothSpan) - mean(tempCh));
        end
    end
end

% plot(timecell{encoderInd},s*sSlope);
% plot(timecell{encoderInd},v*vSlope);

strings = cell(length(cardToShow)*length(chSN) + encoderVelDis,1);

for j = 1:length(cardToShow)
    for i = 1:length(chSN)
        strings{(j - 1)*length(chSN) + i} ...
            = sprintf(['DTAQC ',num2str(cardToShow(j)),' Ch ',num2str(chSN(j,i))]);
    end
end

% strings{end-2} = [['smooth = ';
strings{end-1} = 'counter';
strings{end} = 'velocity';

legend (strings);

xlabel ('Time in seconds');
ylabel (['Velocity * ',num2str(1/vSlope,'%.1e'),...
    ' m/s',sprintf('\n'),'Distance * ',num2str(1/sSlope,'%.1e'),' m']);

hold off;