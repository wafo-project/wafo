function [rww,rwx,rxx]=spec2lcov(spec,t,u,type,alpha,beta)
%SPEC2LCOV Calculates auto- and cross-covariance functions 
%       between W(0,0) and X(t,u) 
%       or selected derivatives in the 2D Lagrange model 
%
%CALL: [rww,rwx,rxx]=spec2lcov(spec,t,u,type,alpha,beta)%
%   For type==1
%       rww      = structure with fields  rww.R, rww.t, rww.u
%                  with covariance values Cov(W(0,0),W(t,u))
%       rwx, rxx = similar structures with  Cov(W(0,0),X(t,u)) and 
%                  Cov(X(0,0),X(t,u))
%   For the other types the three covariance structures contain the
%   covariance and cross-covariances indicated below
%       spec     = orbital spectrum structure with depth  h
%       t        = vector of time values
%       u        = vector of space values
%       type     = 1,2,3,4,5,6,7
%       alpha,beta = parameters in the linked model 
%
%   if type=1: W(0,0),X(0,0) and W(t,u),X(t,u) 
%           gives r^(ww),r^(wx),r^(xx)
%   if type=2: dW/dt(0,0),dX/dt(0,0) and W(t,u),X(t,u) 
%           gives r^(ww)_(t0),r^(wx)_(t0),r^(xx)_(t0)
%   if type=3: dW/du(0,0),dX/du(0,0) and W(t,u),X(t,u) 
%           gives r^(ww)_(u0),r^(wx)_(u0),r^(xx)_(u0)
%   if type=4: dW/dt(0,0),dX/dt(0,0) and dW/dt(t,u),dX/dt(t,u) 
%           gives r^(ww)_(tt),r^(wx)_(tt),r^(xx)_(tt)
%   if type=5: dW/du(0,0),dX/du(0,0) and dW/du(t,u),dX/du(t,u) 
%           gives r^(ww)_(uu),r^(wx)_(uu),r^(xx)_(uu)
%   if type=6: dW/dt(0,0),dX/du(0,0) and dW/dt(t,u),dX/du(t,u) 
%           gives r^(ww)_(tu),r^(wx)_(tu),r^(xx)_(tu)
%
% NOTE: This routine works only for one-sided spectra

% Tested on Matlab 8.1, 8.6
% History
% Written by Sofia Äberg and Georg Lindgren for use in Lagrange papers

S=spec;
if strcmpi(S.type,'k1d'),
    S=spec2spec(S,'freq');
end

if ~isfield(S,'h')|isempty('S.h')
    h=Inf;
else
    h=S.h;
end

if ~isfield(S,'g')|isempty('S.g')
    g=gravity;
else
    g=S.g;
end

m=spec2mom(S,4);
if (nargin < 2)|isempty(t),
    tmax=2*pi*sqrt(m(1)/m(2));
    t=(0:100)/100*3*floor(tmax);
end
t=t';

w=S.w;
k=w2k(w,0,h,g);
w_=w*ones(1,length(t));
k_=k*ones(1,length(t));

% Compute transfer function
H1=(tanh(h*k_(2:end,:))).^(-1);
%H1=(tanh(h*k_(2:end,:))).^(-1);
H2=(-alpha+1i*beta*k_(2:end,:))./(w_(2:end,:).^2);
H=1i*H1-H2;
%plot(real(H(20:end)),imag(H(20:end)),'.')
H=[zeros(1,length(t));H];
theta=angle(H);
theta(1,:)=pi;
rho=abs(H);
S_=S.S*ones(1,length(t));
S_div_tanh=(S.S(2:end)*ones(1,length(t)))./tanh(h*k_(2:end,:));
S_div_tanh=[zeros(1,length(t));S_div_tanh];
S_div_tanh_2=(S.S(2:end)*ones(1,length(t)))./tanh(h*k_(2:end,:)).^2;
S_div_tanh_2=[zeros(1,length(t));S_div_tanh_2];

%M1=max(abs(S_.*rho-S_div_tanh))
%M2=max(abs(S_.*rho.^2-S_div_tanh_2))
Sx=createspec('freq','w');
Sx.w=w;
Sx.S=S_div_tanh(:,1);
mx=spec2mom(Sx,2);

if (nargin < 3)||isempty(u),
    umax=2*pi*sqrt(mx(1)/mx(2));
    u=(-100:100)/100*5*floor(umax);
    u=u';
end

% For each  u  make a vector function of  w  with 
% number of columns = length of  t
% Then loop over  u  and integrate

rww_=[]; rwx_=[];rxx_=[];
rwx__=rwx_;

if type==6,
  for iu=1:length(u),
    ui=u(iu);
    f1=-(w_.*k_).*cos(k_*ui-w*t').*S_;
    f2=-(w_.*k_).*cos(k_*ui-w*t'+theta).*rho.*S_;
    f3=-(w_.*k_).*cos(k_*ui-w*t').*rho.^2.*S_;
    I1=simpson(w,f1);
    I2=simpson(w,f2);
    I3=simpson(w,f3);
    rww_=[rww_;I1];
    rwx_=[rwx_;I2];
    rxx_=[rxx_;I3];
  end

elseif type==5,
  for iu=1:length(u),
    ui=u(iu);
    f1=(k_.^2).*cos(k_*ui-w*t').*S_;
    f2=(k_.^2).*cos(k_*ui-w*t'+theta).*rho.*S_;
    f3=(k_.^2).*cos(k_*ui-w*t').*rho.^2.*S_;
    I1=simpson(w,f1);
    I2=simpson(w,f2);
    I3=simpson(w,f3);
    rww_=[rww_;I1];
    rwx_=[rwx_;I2];
    rxx_=[rxx_;I3];
  end

elseif type==4,
  for iu=1:length(u),
    ui=u(iu);
    f1=(w_.^2).*cos(k_*ui-w*t').*S_;
    f2=(w_.^2).*cos(k_*ui-w*t'+theta).*rho.*S_;
    f3=(w_.^2).*cos(k_*ui-w*t').*rho.^2.*S_;
    I1=simpson(w,f1);
    I2=simpson(w,f2);
    I3=simpson(w,f3);
    rww_=[rww_;I1];
    rwx_=[rwx_;I2];
    rxx_=[rxx_;I3];
  end

elseif type==3,
  for iu=1:length(u),
    ui=u(iu);
    f1=k_.*sin(k_*ui-w*t').*S_;
    f2=k_.*sin(k_*ui-w*t'+theta).*rho.*S_;
    f3=k_.*sin(k_*ui-w*t').*rho.^2.*S_;
    I1=simpson(w,f1);
    I2=simpson(w,f2);
    I3=simpson(w,f3);
    rww_=[rww_;I1];
    rwx_=[rwx_;I2];
    rxx_=[rxx_;I3];
  end

elseif type==2,
  for iu=1:length(u),
    ui=u(iu);
    f1=-w_.*sin(k_*ui-w*t').*S_;  
    f2=-w_.*sin(k_*ui-w*t'+theta).*rho.*S_;
    f3=-w_.*sin(k_*ui-w*t').*rho.^2.*S_;
    I1=simpson(w,f1);
    I2=simpson(w,f2);
    I3=simpson(w,f3);
    rww_=[rww_;I1];
    rwx_=[rwx_;I2];
    rxx_=[rxx_;I3];
  end

elseif type==1,
  for iu=1:length(u),
    ui=u(iu);   
    f1=cos(k_*ui-w*t').*(S.S*ones(1,length(t)));
    f2=cos(k_*ui-w*t'+theta).*rho.*S_;
    f3=cos(k_*ui-w*t').*rho.^2.*S_;
    I1=simpson(w,f1);
    I2=simpson(w,f2);
    I3=simpson(w,f3);
    rww_=[rww_;I1];
    rwx_=[rwx_;I2];
    rxx_=[rxx_;I3];
  end
 
else disp('Unknown type')
end

rww.R=rww_;
rwx.R=rwx_;
rxx.R=rxx_;
rww.t=t;
rwx.t=t;
rxx.t=t;
rww.u=u;
rwx.u=u;
rxx.u=u;

