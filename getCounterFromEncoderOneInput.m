% calculates encoder counts and encoder velocity using csvread which is a
% relatively fast method, please write your own code because this only
% serves as an algorithm reference

function [counter,v] = getCounterFromEncoderOneInput(test)

lenPerSector = 1.5e-6;  % in meters
f = 1e6;    % frequency in Hz

count = 0;
encoderShift = 3;
temp = zeros(1,length(test));
shiftedEncoder = test(:,1) - encoderShift;
for j=1:length(test(:,1))-1
    if (shiftedEncoder(j))*(shiftedEncoder(j+1))<0
        count = count + 1;
    end
    temp(j+1) = count;
end

counter = temp * lenPerSector;

%%
% averaging method of calc velocity
intvl = 1e5;
tempv = zeros(1,length(test));

denom = intvl/(f*lenPerSector);
for k=intvl+1:length(temp)-intvl
    tempv(1,k) = -(temp(k-intvl)-temp(k+intvl))/denom;
end

v = tempv;
