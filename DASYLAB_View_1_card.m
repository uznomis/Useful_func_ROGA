% this script reads data produced by DASYLAB on 'one' MCC card only

%% Import
headerlinesToSkip = 10;

[filename, filepath] = uigetfile...
    ('*.*', 'Please select file from DASYLAB');
if filename == 0
    return
end

raw_data = dlmread([filepath,filename],'',headerlinesToSkip,0);

%% Plot
plot(raw_data);
return

%% Export to PSI
try
    dcm_obj = datacursormode(gcf);
    c_info = getCursorInfo(dcm_obj);
    c_info = fliplr(c_info);
    xInd = ones(1,2);
    for i = 1:2
        xInd(i) = c_info(i).Position(1);
    end
catch ME
    error('Error. Please see ME for details.');
end

[~,name,~] = fileparts(filename);
dt = datestr(now,'mm_dd_yyyy_HH_MM_SS');
toExport = raw_data(xInd(1):xInd(2),:);

save([filepath, 'Crop_', name, '_', dt, '.txt'], 'toExport', '-ascii');
msgbox('Done.');