function Slope = spec2slopedat3D(Spec,Nsim,Points,options,varargin)
%SPEC2SLOPEDAT3D Simulates values and slopes in 3D Lagrange field 
%
%     with choice between 
%     GL, fft2 over space, stepping over t, if Spec is dirspec
%     MP, fft in time, stepping over space, if Spec is onedim with field D               

% Tested on Matlab 8-1, 8.6
% History:
% By GL October 2015 for simple comparison between fft-methods

if nargin>4,
    opt = simoptset(options,varargin{:});
else
    opt = options;
end
opt3D = genoptset('type','timeseries','PP',Points)
S = Spec;
PP = Points;

if isfield(S,'D'),
    method = 'MP';
elseif strcmpi(S.type(1:3),'dir')
    method = 'GL';
else
    disp('Unknown spec type')
end

[M,Np] = size(PP); 
if (M ~= 2), 
    disp('Points must be size  2 x npoints')
    return
end

Slope = struct('dat',[],'der',[]);
if Np > 1, 
    for np = 1:Np,
        Slope.dat{np} = [];
        Slope.der{np}=[];
    end
end

if strcmpi(method,'MP'),
    for nsim=1:Nsim,
        [W,X,Y] = spec2ldat3DP(S,1,opt);
        L = ldat2lwav3D(W,X,Y,opt3D);
        for np=1:Np,
            Der = gradient(L.Z{np},L.t(2)-L.t(1));
            Slope.dat{np} = [Slope.dat{np}; L.Z{np}];
            Slope.der{np} = [Slope.der{np}; Der];
        end
    end
end

if strcmpi(method,'GL'),
    for nsim=1:Nsim,
        L = spec2lseries(S,PP,opt);
        for np=1:Np,
            Der = gradient(L.Z{np},L.t(2)-L.t(1));
            Slope.dat{np} = [Slope.dat{np}; L.Z{np}];
            Slope.der{np} = [Slope.der{np}; Der];
        end
    end
end



