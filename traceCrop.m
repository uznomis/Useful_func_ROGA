% pick a segment from a full trace; timeCol is the column number of time
% vector, zero if no time vector; if import is a structure then it doesn't
% plot and use the structure as the pick (beware whether import is time
% based or indix based
function [segTrace] = traceCrop(trace,timeCol,import)
% create time vector if necessary and plot if import is not struct
if ~timeCol
    if ~isstruct(import)
        figure;
        plot(trace);
    end
else
    time = trace(:,timeCol);
    if ~isstruct(import)
        figure;
        plot(time,[trace(:,1:timeCol - 1) trace(:,timeCol + 1: end)]);
    end
end

% read cursor info from either workspace or input
if ~isstruct(import)
    pause();    % pause for manual picking
    cursor = evalin('base','cursor_info');
else
    cursor = import;
end

% crop trace
x = ones(1,length(cursor));
for i = 1:length(cursor)
    x(i) = cursor(i).Position(1);
end

xMin = min(x);
xMax = max(x);
if ~timeCol
    segTrace = trace(xMin:xMax,:);
else
    minInd = find(time >= xMin, 1);
    maxInd = find(time >= xMax, 1);
%     segTrace = trace(minInd:maxInd,...
%         [1:timeCol-1 timeCol+1:end]);   % excluding time
    segTrace = trace(minInd:maxInd,:);
end