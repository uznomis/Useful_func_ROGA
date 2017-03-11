% This function returns a timeseries (see help) of time,distance,axial,shear,
% dilation,friction as a function of distance (streched evenly).

% Modified on 10072016

% This script should be always kept in a same folder with importfile_for_distance_stretch.m!
% This code stretches the distance data column to an evenly spaced format for spatial Fourier
% Transform. Spacing is 0.1 mm and can be changed below. It extracts time,distance,axial,shear,
% dilation,friction and generates a .c6 file
% with those parameters. The .c6 file has a row number that is truncated to a power of 2 that
% is nearest to the original row number. A MAT file can also be generated if the line with save
% command is uncommented.

function data_res = stretchDistanceEvenly()

[name,path] = uigetfile('*.c40','Select a C40 file');
if name == 0
    return;
end

[~,dryName,~] = fileparts(name);
[time,distance,axial,shear,dilation,friction]...
    = importfile_for_distance_stretch([path,name]);

% get rid of NaN's in the beggining and end of the column vectors (file
% specific NaN's in *.c40)
t = time(2:end,1);
d_end = distance(end-1,1);
d = [0;distance(3:end-1,1);d_end];
ax = axial(2:end,1);
sh = shear(2:end,1);
di = dilation(2:end,1);
fr = friction(2:end,1);

data = [t ax sh di fr];

% forming a timeseries as a function of distance and resample evenly
ts = timeseries(data,d);

dd = 0.1e-3; % spacing of distance stretch
new_d = d(1):dd:d(end); % new vector of time with fixed interval

data_res = resample(ts,new_d); % data resampled at constant interval

% export to a .mat file
% distance = data_res.time;
% time = data_res.data(:,1);
% axial = data_res.data(:,2);
% shear = data_res.data(:,3);
% dilation = data_res.data(:,4);
% friction = data_res.data(:,5);
% save([path,'Distance_domain_',dryName,'.mat'],...
% 'distance','time','axial','shear','dilation','friction');   % commented unless MAT output is needed

% export to a txt file *.c6
cut_l = 2^(nextpow2(length(new_d))-1);  % cut to power of 2 for fft purpose

fid = fopen([path,'Distance_domain_',dryName,'.c6'],'w');
if fid == -1
    return;
end

fprintf(fid,'%s,%s,%s,%s,%s,%s\n','D_Distance_m','D_Time_s','D_AxialStress_MPa',...
    'D_ShearStress_MPa','D_Dilation_microns','D_Friction'); % printing titles

fprintf(fid,'%10e,%10e,%10e,%10e,%10e,%10e\n',[data_res.time(1:cut_l,:)';...
    data_res.data(1:cut_l,:)']);    % printing data

fclose(fid);

% below is an auto generated importing function by MATLAB
function [Time_s1,Distance_m1,AxialStress_MPa1,ShearStress_MPa,Dilation_microns,Friction] = importfile_for_distance_stretch(filename, startRow, endRow)
%% Initialize variables.
delimiter = '\t';
if nargin<=2
    startRow = 1;
    endRow = inf;
end

%% Read columns of data as strings:
% For more information, see the TEXTSCAN documentation.
formatSpec = '%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%s%s%*s%s%s%s%s%[^\n\r]';

%% Open the text file.
fileID = fopen(filename,'r');

%% Read columns of data according to format string.
% This call is based on the structure of the file used to generate this
% code. If an error occurs for a different file, try regenerating the code
% from the Import Tool.
dataArray = textscan(fileID, formatSpec, endRow(1)-startRow(1)+1, 'Delimiter', delimiter, 'HeaderLines', startRow(1)-1, 'ReturnOnError', false);
for block=2:length(startRow)
    frewind(fileID);
    dataArrayBlock = textscan(fileID, formatSpec, endRow(block)-startRow(block)+1, 'Delimiter', delimiter, 'HeaderLines', startRow(block)-1, 'ReturnOnError', false);
    for col=1:length(dataArray)
        dataArray{col} = [dataArray{col};dataArrayBlock{col}];
    end
end

%% Close the text file.
fclose(fileID);

%% Convert the contents of columns containing numeric strings to numbers.
% Replace non-numeric strings with NaN.
raw = repmat({''},length(dataArray{1}),length(dataArray)-1);
for col=1:length(dataArray)-1
    raw(1:length(dataArray{col}),col) = dataArray{col};
end
numericData = NaN(size(dataArray{1},1),size(dataArray,2));

for col=[1,2,3,4,5,6]
    % Converts strings in the input cell array to numbers. Replaced non-numeric
    % strings with NaN.
    rawData = dataArray{col};
    for row=1:size(rawData, 1);
        % Create a regular expression to detect and remove non-numeric prefixes and
        % suffixes.
        regexstr = '(?<prefix>.*?)(?<numbers>([-]*(\d+[\,]*)+[\.]{0,1}\d*[eEdD]{0,1}[-+]*\d*[i]{0,1})|([-]*(\d+[\,]*)*[\.]{1,1}\d+[eEdD]{0,1}[-+]*\d*[i]{0,1}))(?<suffix>.*)';
        try
            result = regexp(rawData{row}, regexstr, 'names');
            numbers = result.numbers;
            
            % Detected commas in non-thousand locations.
            invalidThousandsSeparator = false;
            if any(numbers==',');
                thousandsRegExp = '^\d+?(\,\d{3})*\.{0,1}\d*$';
                if isempty(regexp(thousandsRegExp, ',', 'once'));
                    numbers = NaN;
                    invalidThousandsSeparator = true;
                end
            end
            % Convert numeric strings to numbers.
            if ~invalidThousandsSeparator;
                numbers = textscan(strrep(numbers, ',', ''), '%f');
                numericData(row, col) = numbers{1};
                raw{row, col} = numbers{1};
            end
        catch me
        end
    end
end


%% Replace non-numeric cells with NaN
R = cellfun(@(x) ~isnumeric(x) && ~islogical(x),raw); % Find non-numeric cells
raw(R) = {NaN}; % Replace non-numeric cells

%% Allocate imported array to column variable names
Time_s1 = cell2mat(raw(:, 1));
Distance_m1 = cell2mat(raw(:, 2));
AxialStress_MPa1 = cell2mat(raw(:, 3));
ShearStress_MPa = cell2mat(raw(:, 4));
Dilation_microns = cell2mat(raw(:, 5));
Friction = cell2mat(raw(:, 6));