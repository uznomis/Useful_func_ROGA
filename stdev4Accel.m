% calculates standard deviation of t0 in locating ball drop test with 4 3D
% accels, assuming a circular geometry and 4 orthogonal accels at the
% interface of experimental fault (or a fixed depth z), i.e. a 2D problem.

function b = stdev4Accel(xy,t,v,D,z)
R = D/2;
accelloc = [R 0 -R 0;0 R 0 -R];
b = std(t - sqrt((xy(1)-accelloc(1,:)).^2+(xy(2)-accelloc(2,:)).^2+z^2)...
    ./v);