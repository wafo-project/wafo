function Steep = lwav2frontback(L)
%LWAV2FRONTBACK Gives front/back crest periods/wavelength of wave data
%
%CALL: Steep =lwav2frontback(L)
%
%   Steep = struct array with fields
%      .ffull = full front period/wavelength
%      .fhalf = half front period/wavelength
%      .bfull = full back period/wavelength
%      .bhalf = half back period/wavelength
%
%   L     = 2D wave L.t (L.u) and L.Z

% Tested on Matlab 8.1, 8.6
% History
%   Created 2014-02-26 by Georg Lindgren

Steep=[];
if isfield(L,'t')
    x=L.t;
    Steep.type='time';
elseif isfield(L,'u')
    x=L.u;
    Steep.type='space';
else
    Steep=[];    
    disp('unknown type')
end
x=reshape(x,length(x),1);
y=reshape(L.Z,length(L.Z),1);
data=[x y-mean(y)];
[TC, tc_ind, v_ind]=dat2tc(data,0,'tw');

% First crest-trough must be a trough (automatic with option 'tw')
if TC(1,2)>TC(2,2),
    TC=TC(2:end,:);
    tc_ind=tc_ind(2:end);
    v_ind=v_ind(2:end);
end
% First v-crossing must be between trough and crest
if v_ind(1)<tc_ind(1),
    v_ind=v_ind(2:end);
end
% Last crest-trough must be a trough
if TC(end,2)>TC(end-1,2),
    TC=TC(1:end-1,:);
    tc_ind=tc_ind(1:end-1);
    v_ind=v_ind(1:end-1);
end
% Last v-crossing must be between crest and trough
if v_ind(end)>=tc_ind(end),
    v_ind=v_ind(1:end-1);
end

Tr=TC(1:2:end,:);
t_ind=tc_ind(1:2:end);
Cr=TC(2:2:end,:);
c_ind=tc_ind(2:2:end);

Up_ind=v_ind(1:2:end);
Do_ind=v_ind(2:2:end);

try
Steep.ffull=x(c_ind)-x(t_ind(1:end-1));
Steep.bfull=x(t_ind(2:end))-x(c_ind);
Steep.fhalf=x(c_ind)-x(v_ind(1:2:end));
Steep.bhalf=x(v_ind(2:2:end))-x(c_ind);
catch
    size(c_ind)
    size(t_ind)
    size(v_ind)
    [c_ind(1) c_ind(end)]
    [t_ind(1) t_ind(end)]
    [v_ind(1) v_ind(end)]
end


