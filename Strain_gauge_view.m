% this program loads data produced my ROGA Terminal developed by James
% Jeffers, display the data according to the user request
%% Importing
for i = 1:2
    [filename{i}, filepath{i}] = uigetfile...
        ('C:\Users\User\Documents\Test files results on fast cards\*.*',...
        ['Please select ACQ100', num2str(i)]);
    if filename(i) == 0
        return
    end
end
for i = 1:2
    test{i} = dlmread([filepath{i}, filename{i}]);
    test{i} = reshape(test{i}, [length(test{i})/17 17]);
end
msgbox('Import successful!');

%% Sync cards
% please provide some important info here
encoderAvailable = 1;
encoderChannel = 16;
freq = 1e6;
% parameters regarding synchronization; no need to change
encoderShift = 3;
encoderSpacingThrsh = 5e4;
encoderXcorrWindow = 5e5;

% execution begins here
if filename{1} ~= 0 && filename{2} ~= 0
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
            s2 = test{i}(encoderjump: encoderjump + encoderXcorrWindow,1);
            [r,lag] = xcorr(s1,s2);
            [~,I] = max(abs(r));
            encoderjumpavg(i-1) = round(-lag(I));
        end
    end
    % create time col vector
    timecell = cell(1,length(filname));
    rawTcell{1} = 1/freq * (1:length(data1))';
    timecell{1} = rawTcell{1};
    
    for i=2:length(filename)
        timecell{i} = ...
            timecell{1}-encoderjumpavg(i-1)/freq; % apply time shift
    end
    
end

%% Plotting
% after importing the data file, simply change parameters and run this
% section with 'Run and Advance' button;

% change paramters below accordingly before continuing
plotInd = [1:17];    % channels to plot; if channel 1 is encoder, [2 3 4] will display data without encoder
legendNames = {'encoder'};    % channels names
medFilterOn = 1;    % choosing 1 makes median filter on; median filter gets rid of spikes

% plot
figure('Name',filename);
for i = 1:2
    if medFilterOn
        toPlot = medfilt2(test{i}(:,plotInd),[3,1]);
    else
        toPlot = test{i}(:,plotInd);
    end
    plot(toPlot);
    hold on
end
hold off
legend(legendNames);