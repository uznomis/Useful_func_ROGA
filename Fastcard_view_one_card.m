% load one file from fast card recording and display the four data
% channels;
%% Importing
[filename, filepath] = uigetfile...
    ('C:\Users\User\Documents\Test files results on fast cards\*.*');
if filename == 0
    return
end
test = dlmread([filepath, filename]);
%% Plotting
% after importing the data file, simply change parameters and run this
% section with 'Run and Advance' button;

% change paramters below accordingly before continuing
plotInd = [8];    % channels to plot; if channel 1 is encoder, [2 3 4] will display data without encoder
legendNames = {'encoder'};    % channels names
medFilterOn = 0;    % choosing 1 makes median filter on; median filter gets rid of spikes

% plot
figure('Name',filename);
if medFilterOn
    toPlot = medfilt2(test(:,plotInd),[3,1]);
else
    toPlot = test(:,plotInd);
end
plot(toPlot);
legend(legendNames);