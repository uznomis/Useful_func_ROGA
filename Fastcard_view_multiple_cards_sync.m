% Fast card data viewing app; automatic sync between cards if there is
% encoder channel 1)WITH slip, 2)on channel 1;
%% Initializing
% paramters to change before continuing
cardSN = [2 3 4];    % which cards to import; numbers should match file names
freq = 1e6;    % freqeuncy in Hz
chopLastSec = 0;    % 1 for yes, 0 for no; sometimes the last sec of data is not recorded
cardToGetEncoder = 2;    % the number of card whose encoder data is used for velocity and counter calculation
encoderAvailable = 1;    % 1 for yes, 0 for no; if no encoder is available, there is no sync between cards
lenPerSector = 1.5e-6;    % encoder sector length in meters

%% Importing
[filename, filepath] = uigetfile...
    ('C:\Users\User\Documents\Test files results on fast cards\*.*');
if filename == 0
    return
end
[~,fname,ext] = fileparts(filename);
test = cell(1,length(cardSN));
for i = 1:length(cardSN)
    tmp = dlmread([filepath, fname(1:end-1), num2str(cardSN(i)),ext]);
    test{i} = tmp(1:end-freq*chopLastSec,:);
end
%% Velocity and counter calculation
% calculate s and v from encoder; calling another function
if cardToGetEncoder == 0
    s = zeros(1,length(test{1}));
    v = s;
    encoderInd = 1;
else
    encoderInd = find(cardSN == cardToGetEncoder, 1);
    [s,v] = getCounterFromEncoder(test{encoderInd},freq,lenPerSector);
end
%% Sync cards
% parameters regarding synchronization; no need to change
encoderShift = 3;
encoderSpacingThrsh = 5e4;
encoderXcorrWindow = 5e5;

% execution begins here
if filename ~= 0
    % calc timeshift of cards
    encodercurrent = 1;
    encoderbefore = 0;
    encoderjump = 1;
    spacing = 0;
    data1 = test{1};
    encoderjumpavg = zeros(1,length(cardSN) - 1);
    if encoderAvailable == 1
        for j = 2:length(data1)
            if (data1(j-1,1) - encoderShift)*...
                    (data1(j,1) - encoderShift) < 0
                encoderbefore = encodercurrent;
                encodercurrent = j;
                spacing = encodercurrent - encoderbefore;
                if spacing < encoderSpacingThrsh && encoderbefore > 1
                    encoderjump = encodercurrent;
                    break;
                end
            end
        end
        s1 = data1(encoderjump:encoderjump + encoderXcorrWindow,1);
        
        for i=2:length(cardSN)
            s2 = test{i}(encoderjump: encoderjump + encoderXcorrWindow,1);
            [r,lag] = xcorr(s1,s2);
            [~,I] = max(abs(r));
            encoderjumpavg(i-1) = round(-lag(I));
        end
    end
    % create time col vector
    timecell = cell(1,length(cardSN));
    rawTcell{1} = 1/freq * (1:length(data1))';
    timecell{1} = rawTcell{1};
    
    for i=2:length(cardSN)
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
cardToShow = [2 3 4];    % cards to display
chSN = [1 2 3 4];    % channels to display on each card
smoothSpan = 1;    % smoothening window; make it 1 to diable smoothening; note there is no smoothening for AE data
cardOffset = 0.5;    % offset between cards when plotting
chOffset = 0.05;    % offset between channels within each card when plotting
sSlope = 1e2;    % slope of counter for scaling when plotting
vSlope = 1e2;    % slope of velocity for scaling when plotting
encoderVelDis = 2;    % a constant; don't change
accelCard = [0];    % card numbers of AE channels; make [0] if no AE data available
accelSlope = 3;    % slope of velocity for scaling when plotting
medFilterOn = 0;    % choosing 1 makes median filter on; median filter gets rid of spikes

% execution begins here
figure('Name',filename);
hold on;

for j = 1:length(cardToShow)
    tempInd = find(cardSN == cardToShow(j),1);
    for i = 1:length(chSN)
        if medFilterOn
            tempCh = medfilt2(test{tempInd}(:,chSN(i)),[3,1]);
        else
            tempCh = test{tempInd}(:,chSN(i));
        end
        if ismember(cardToShow(j),accelCard)            
            plot(timecell{tempInd},j*...
                cardOffset + i*chOffset + accelSlope*(tempCh - mean(tempCh)));
        else
            plot(timecell{tempInd},j*...
                cardOffset + i*chOffset + smooth(tempCh,smoothSpan) - mean(tempCh));
        end
    end
end

plot(timecell{encoderInd},s*sSlope);
plot(timecell{encoderInd},v*vSlope);

strings = cell(length(cardToShow)*length(chSN) + encoderVelDis,1);

for j = 1:length(cardToShow)
    for i = 1:length(chSN)
        strings{(j - 1)*length(chSN) + i} ...
            = sprintf(['Card ',num2str(cardToShow(j)),' Ch ',num2str(chSN(i))]);
    end
end
strings{end-1} = 'counter';
strings{end} = 'velocity';

legend (strings);

xlabel ('Time in seconds');
ylabel (['Velocity * ',num2str(1/vSlope,'%.1e'),...
    ' m/s',sprintf('\n'),'Distance * ',num2str(1/sSlope,'%.1e'),' m']);

hold off;