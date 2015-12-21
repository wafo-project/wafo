function [Qr,QQ,FF,FFr] = mctp2reverse(F)
% MCTP2REVERSE  Calculates the time-reversed MCTP for a SMCTP.
%
% CALL: [Qr] = mctp2reverse(F);
%       [Qr,QQ,FF,FFr] = mctp2reverse(F);
%
% Qr     = Cell array of min-max and max-min transition 
%          matrices for time-reversed MCTP.           {1x2}
% QQ     = Cell array of min-max and max-min transition 
%          matrices for MCTP.                         {1x2}
% FF     = From-To matrix calculated from QQ.         [nxn]
% FFr    = From-To matrix calculated from Qr.         [nxn]
%
% F{1,1} = min-Max matrix                             [nxn]
% F{1,2} = Max-min matrix                             [nxn]
%
% If a matrix F{1,2}=[], then the process will
% be assumed to be time-reversible.
%
% Se also 

% Copyright (c) 1999 by Pär Johannesson
% Toolbox: Rainflow Cycles for Switching Processes V.1.0, 2-Oct-1997

% Check input arguments

ni = nargin;
no = nargout;
error(nargchk(1,1,ni));

% Define 

n = length(F{1,1});  % Number of levels

% Normalize the rowsums of F{1,1} to 1  ==>  QQ{1,1}
% Normalize the rowsums of F{1,2} to 1  ==>  QQ{1,2}

QQ= cell(1,2);

QQ{1,1} = F{1,1};      % min-max matrix
if isempty(F{1,2})     % Time-reversible?
  QQ{1,2} = F{1,1}';   % max-min matrix
else                   % F{i,2} is given
  QQ{1,2} = F{1,2};    % max-min matrix
end
    
QQ{1,1} = mat2tmat(QQ{1,1},1);  % normalize min-max matrix
QQ{1,2} = mat2tmat(QQ{1,2},-1); % normalize max-min matrix

%
% Create Transition matrices for time-reversed MCTP
%

Qr = cell(1,2);

% Calculate stationary distribution of minima and maxima
[ro,roh] = mctp2stat(F);

% Backward min-to-max
I1 = find(ro>0); I2 = find(ro<=0);
ro_inv = ro; ro_inv(I1) = 1./ro(I1); ro_inv(I2) = zeros(1,length(I2));
Qr{1,1} = QQ{1,2}' .* (ro_inv'*roh);

% Backward max-to-min
I1 = find(roh>0); I2 = find(roh<=0);
roh_inv = roh; roh_inv(I1) = 1./roh(I1); roh_inv(I2) = zeros(1,length(I2));
Qr{1,2} = QQ{1,1}' .* (roh_inv'*ro);

% Make the frequency matrix FF for the joint min-Max and Max-min
% distribution (from Q)

FF = QQ{1,1}.*(ro'*ones(1,n)) + QQ{1,2}.*(roh'*ones(1,n));

% Make the frequency matrix FF for the joint min-Max and Max-min
% distribution (from Qr)

FFr = Qr{1,1}.*(ro'*ones(1,n)) + Qr{1,2}.*(roh'*ones(1,n));

