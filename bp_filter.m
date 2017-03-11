function  [o] =  bp_filter(d,dt,f1,f2,f3,f4);
% BP_Filter: apply a BP filter with corner freq. given
%            ny f1,f2,f3,f4. FFT construction (not good
%            for short time series)
% 
%  [o] = bp_filter(d,dt,f1,f2,f3,f4);
%
%  IN     d:    data (columns are traces)
%         dt:   sampling interval in sec
%         f1:   freq. in Hz
%         f2:   freq. in Hz
%         f3:   freq. in Hz
%         f4:   freq. in Hz
%
%   ^
%   |     ___________
%   |    /           \   Amplitude spectrum
%   |   /             \
%   |  /               \
%   |------------------------>
%      f1 f2        f3 f4
%
%  OUT    o:    output  (columns are traces)
%
%
%  Author(s): M.D.Sacchi (sacchi@phys.ualberta.ca)
%  Copyright 1998-2003 SeismicLab
%  Revision: 1.2  Date: Dec/2002
%
%  Signal Analysis and Imaging Group (SAIG)
%  Department of Physics, UofA
%


 [nt,nx] = size(d);
 k = nextpow2(nt);
 nf = 4*(2^k);

 i1 = floor(nf*f1*dt)+1;
 i2 = floor(nf*f2*dt)+1;
 i3 = floor(nf*f3*dt)+1;
 i4 = floor(nf*f4*dt)+1;

 up =  (1:1:(i2-i1))/(i2-i1);
 down = (i4-i3:-1:1)/(i4-i3);
 aux = [zeros(1,i1), up, ones(1,i3-i2), down, zeros(1,nf/2+1-i4) ];
 aux2 = fliplr(aux(1,2:nf/2));

 c = 0; % zero phase (could apply rotations as well)
 F = ([aux,aux2]');
 Phase = (pi/180.)*[0.,-c*ones(1,nf/2-1),0.,c*ones(1,nf/2-1)];
 Transfer = F.*exp(-i*Phase');


 D = fft(d,nf,1);

 for k = 1:nx
  Do(:,k) = Transfer.*D(:,k);
 end

 o = ifft(Do,nf,1);

 o = real(o(1:nt,:));


 
