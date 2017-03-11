% This script contains three sections; first section loads accel data
% single file and plots raw data; second section lets you select a window
% and show the intergration of the data within the window; third section is
% a feature of exporting data indices with cursor_info; see descriptions in
% the beginning of each section for details;
%%
plotInd = [1 2 3 4];
medFilterOn = 0;
[filename, filepath] = uigetfile...
    ('C:\Users\User\Documents\Test files results on fast cards\*.*');

if filename ~= 0
    test = dlmread([filepath, filename]);
    figure;
    if medFilterOn
        test = medfilt2(test(:,plotInd),[3,1]);
    else
        test = test(:,plotInd);
    end
    plot(test);
    legend('show');
end
return

%%
% Plots velocity curves of the original signal. To use, two right clicks to
% select waveform window, and right click + left click to select background
% level, and that gives a new plot of velocity. Must have 4 channels.

encoderCh = 1;    % encoder channel number, usually 1

if filename ~= 0
    
    x1 = 1;  % set up a flag to indicate exit or not
    pickedIntvl = 0; % a flag to show whether backgrd is picked already
    
    while (x1 > 0)
        
        [x,y,ibut] = ginput (2);
        
        if (ibut == 1)   % left clicks to zoom horizontally
            xlim([min(x) max(x)]);
        end
        
        if (ibut == 2)  % middle click to exit
            x1=0;
        end
        
        if (ibut(1) == 1 & ibut(2) == 3)    % left + right clicks to auto zoom
            xlim('auto');
        end
        
        if (ibut == 3)   % right clicks to select interval to display
            initial = round(min(x)); final = round(max(x));
            pickedIntvl = 1;
        end
        
        if (ibut(1) == 3 & ibut(2) == 1)    % right + left clicks to display velocity curves
            
            offset = 20;    % offset between channels on plot for easy viewing
            f = 1e6;    % frequency of recording
            test_v = [];    %initialize an array to hold vel curves
            startbackgrd = round(min(x)); endbackgrd = round(max(x));
            
            if pickedIntvl == 0
                warning('Please pick intervals first!');
                continue
            end
            
            for j = 1:length(plotInd) % calculates velocity curves and demean
                if (j == encoderCh)
                    test_v(:,j) = test(initial:final,j).*offset/f;
                else
                    test_v(:,j) = cumtrapz(test(initial:final,j)...
                        -mean(test(startbackgrd:endbackgrd,j)))/f;
                end
            end
            
            figure; % nes figure
            
            for k=1:length(plotInd)   % offset curves
                test_v(:,k)=test_v(:,k)+(k-1)*offset/f;
            end
            
            plot(test_v);
            
            for l=1:length(plotInd)   % draw baselines
                hold on;
                temp_len = length(test_v);
                y_line = (l-1)*offset/f*ones(1,temp_len);
                plot(1:temp_len,y_line,'LineStyle','--','Color','k');
                hold off;
            end
            
            legend('show');
            
            x1=0;   % set flag to 0 and exit
        end
        
        
    end
end
return
%%
% output 4 indices to variable time_raw
Index = [];
for i=1:length(cursor_info);
    Index(i,1) = cursor_info(i).DataIndex;
end;
time_raw = flipud(Index)    % wiht no semicolon just to display the indices
return