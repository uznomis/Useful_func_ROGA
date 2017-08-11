function output = calculateStrainXYZ(input,theta)
exx = input(:,3);
eyy = input(:,1);
exy = (2*input(:,2)-input(:,1)-input(:,3))./2;
t = 180 - theta;
L = [cosd(t) sind(t);-sind(t) cosd(t)];
[XX, YY, XY] = arrayfun(@(xx,yy,xy) transformXYZ(xx,yy,xy,L), exx, eyy, exy);
output = [XY YY XX];
    
function [XX,YY,XY] = transformXYZ(exx,eyy,exy,L)
e_ = L*[exx exy;exy eyy]*L';
XX = e_(1,1);
YY = e_(1,2);
XY = e_(2,2);