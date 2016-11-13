function [F_rfc,FF,FFr,T] = mctp2arfm(F,c_m)
%MCTP2ARFM  Calculates asymmetric rainflow matrix for a MCTP.
%
% CALL:  [F_rfc]        = mctp2arfm(F);
%        [F_rfc,FF,FFr] = mctp2arfm(F,c_m);
%
% F_rfc  = Asymmetric Rainflow Matrix                 [n,n]
% FF     = Cell array of min-Max and Max-min matrices {1,2}
% FFr    = Cell array of min-Max and Max-min matrices {1,2}
%
% F      = Cell array of min-Max and Max-min matrices {1,2}
% F{1,1} = min-Max matrix                             [n,n]
% F{1,2} = Max-min matrix                             [n,n]
% c_m    = Intensity of local minima, switching proc. [1,1]
%          (Default: 1)
%
% Calculates the expected Asymmetric RainFlow Matrix (ARFM) for a 
% Markov Chain of Turning Points.  Recursive formulation a'la Igor
% If a matrix F{1,2}=[], then the process will be assumed to be 
% time-reversible.
%
% Example: 
%   param = [-1 1 32]; u = levels(param);
%   F = mktestmat(param,[-0.2 0.2],0.15,2);
%   Frfc = mctp2rfm({F []});
%   Farfc = mctp2arfm({F []});
%   cmatplot(u,u,{Frfc Farfc},3);
%   assert(sum(sum(abs((Farfc+Farfc')-(Frfc+Frfc')))), 0, 1e-10) % should be zero
%
%  close all;
%
% See also  arfm2mctp, smctp2arfm, mctp2rfm, cmatplot

% References:
%  
%  P. Johannesson (1999):
%  Rainflow Analysis of Switching Markov Loads.
%  PhD thesis, Mathematical Statistics, Centre for Mathematical Sciences,
%  Lund Institute of Technology.
  
% Tested  on Matlab  5.3
%
% History:
% Revised by PJ  09-Apr-2001
%   updated for WAFO
% Created by PJ (P�r Johannesson) 1998
% Copyright (c) 1997-1998 by P�r Johannesson
% Toolbox: Rainflow Cycles for Switching Processes V.1.1, 22-Jan-1998

% This program used the formulation where the probabilities
% of the events are calculated using "elementary" events for
% the MCTP.
%
% Standing
%    pS = Max*pS1*pS2*pS3;
%    F_rfc(i,j) = pS;
% Hanging
%    pH = Min*pH1*pH2*pH3;
%    F_rfc(j,i) = pH;
%
% The cond. prob. pS1, pS2, pS3, pH1, pH2, pH3 are calculated using
% the elementary cond. prob. C, E, R, D, E3, Ch, Eh, Rh, Dh, E3h. 

T(1,:)=clock;

% Check input arguments

ni = nargin;
no = nargout;
error(nargchk(1,2,ni));

if ni < 2
  c_m=1;
end

% Define 

n = length(F{1,1});  % Number of levels

% Normalize the rowsums of F{1,1} to 1  ==>  Q

Q = F{1,1};     
Q = mat2tmat(Q,1); % Convert to min-max transition matrix

% Normalize the rowsums of F{1,2} to 1  ==>  Qh

if isempty(F{1,2})       % Time-reversible
  Qh = F{1,1}';
else                     % Fh is given
  Qh = F{1,2}; 
end

Qh = mat2tmat(Qh,-1); % Convert to max-min transition matrix

T(2,:)=clock;

% Stationary distribution (=ro) of local minima with transition matrix
% Qt = Q*Qh = "Transition matrix for min-to-min"

Qt = Q*Qh;
ro = mc2stat(Qt(1:(n-1),1:(n-1)));  % Stationary distr., row vector  
ro = [ro 0];  % Minimum can't reach the highest level
I = find(ro<0);
if ~isempty(I)
  warning(['Negative elements in ro. Setting to zero!']);
  ro(I) = zeros(size(I));
end

% Stationary distribution (=roh) of local maxima with transition matrix
% Qt = Qh*Q = "Transition matrix for max-to-max"

Qth = Qh*Q;
roh = mc2stat(Qth(2:n,2:n));  % Stationary distr., row vector  
roh = [0 roh];  % Maximum can't reach the highest level
I = find(roh<0);
if ~isempty(I)
  warning(['Negative elements in roh. Setting to zero!']);
  roh(I) = zeros(size(I));
end

T(3,:)=clock;

% Create Transition matrices for time-reversed MC

% Backward min-to-max
I1 = find(ro>0); 
ro_inv = zeros(1,n); ro_inv(I1) = 1./ro(I1); 
Qr = Qh' .* (ro_inv'*roh);

% Backward max-to-min
I1 = find(roh>0); 
roh_inv = zeros(1,n); roh_inv(I1) = 1./roh(I1); 
Qrh = Q' .* (roh_inv'*ro);

T(4,:)=clock;

% Define matrices for going to or above a level
%R   = fliplr(cumsum(fliplr(Q),2));
%Rr  = fliplr(cumsum(fliplr(Qr),2));

% Define matrices for going to or below a level
%Rh  = cumsum(Qh,2);
%Rrh = cumsum(Qrh,2);

% Compure min-max and max-min matrices for model
FF{1,1} = diag(ro)*Q;     % min-max matrix
FF{1,2} = diag(roh)*Qh;   % max-min matrix
FFr{1,1} = diag(roh)*Qrh; % min-max matrix (from reversed chain)
FFr{1,2} = diag(ro)*Qr;   % min-max matrix (from reversed chain)
% FF and FFr should be identical

T(5,:)=clock;

% Calculate min-max and max-min matrices from input F
%F1 = F{1,1}; 
%F1 = F1/sum(sum(F1));
%F1h = F{1,2}; 
%F1h = F1h/sum(sum(F1h));
% The model FF need not be equal to F
F1=Q; F1h=Qh;

% Calculate rainflow matrix 

F_rfc = zeros(n,n);
EYE = eye(n,n);

%fprintf(1,'Calculating row ');
for i=1:n-1
  %  fprintf(1,'-%1d',i);
  
  for j=i+1:n
    
    xx = FF{1,1}(i,j);
    yy = FF{1,2}(j,i);
    
    x = Q(i,j)*ro(i);
    y = Qh(j,i)*roh(j);
    
    %fprintf(1,'x=%f, y=%f, xx=%f, yy=%f\n',x,y,xx,yy);
    
    %    Min = sum(F1(i,:));
    %    Max = sum(F1h(j,:));
    Min = sum(F1(i,:)+F1h(:,i)')/2;
    Max = sum(F1h(j,:)+F1(:,j)')/2;
    
    Ro = ro(i);    % Probability of "min=i"
    Roh = roh(j);  % Probability of "max=j"
    
    %    fprintf(1,'Min=%f, Max=%f, Ro=%f, Roh=%f\n',Min,Max,Ro,Roh);
    
    Min = Ro; Max = Roh; % Just to be sure
    
    if j-i == 1  % First subdiagonal
      
      C = y/Min; Ch = x/Max;
      D = 0; Dh = 0;
      E = 1-y/Min; Eh = 1-x/Max;
      R = 0; Rh = 0;
      E3 = 1; E3h = 1;
      
    elseif j-i == 2  % Second subdiagonal
      
      % For Part 1 & 2 of cycle
      
      I  = i+1:j-2;
      J  = i+2:j-1;
      
      e  = 1 - sum(Qr(I,i+2:j),2);
      eh = 1 - sum(Qrh(J,i:j-2),2);
      
      C   = y/Min;
      Ch  = x/Max;
      d1  = Qr(i,i+1)*(1-Qrh(i+1,i));
      D   = d1;
      d1h = Qrh(j,j-1)*(1-Qr(j-1,j));
      Dh  = d1h;
      d0  = sum(Qr(i,i+1:j-1));
      E   = 1-d0-y/Min;
      d0h = sum(Qrh(j,i+1:j-1));
      Eh  = 1-d0h-x/Max;
      r1  = Qr(i,i+1)*Qrh(i+1,i);
      R   = r1;
      r1h = Qrh(j,j-1)*Qr(j-1,j);
      Rh  = r1h;
      
      % For Part 3 of cycle
      
      d3  = sum(Q(i,i+1:j-1));
      E3 = 1-d3;
      d3h = sum(Qh(j,i+1:j-1));
      E3h = 1-d3h;
      
    else
      
      Eye = EYE(1:j-i-2,1:j-i-2);
      
      % For Part 1 & 2 of cycle
      
      I  = i+1:j-2;
      J  = i+2:j-1;
      A  = Qr(I,J);
      Ah = Qrh(J,I);
      a  = Qr(i,J);
      ah = Qrh(j,I);
      b  = Qr(I,j);
      bh = Qrh(J,i);
      
      e  = 1 - sum(Qr(I,i+2:j),2);
      eh = 1 - sum(Qrh(J,i:j-2),2);
      
      Inv = inv(Eye-A*Ah);
      C   = y/Min + a*Ah*Inv*b;
      Ch  = x/Max + ah*Inv*A*bh;
      d1  = Qr(i,i+1)*(1-Qrh(i+1,i));
      D   = d1+a*eh+a*Ah*Inv*A*eh;
      d1h = Qrh(j,j-1)*(1-Qr(j-1,j));
      Dh  = d1h+ah*Inv*e;
      d0  = sum(Qr(i,i+1:j-1));
      E   = 1-d0-y/Min+a*Ah*Inv*e;
      d0h = sum(Qrh(j,i+1:j-1));
      Eh  = 1-d0h-x/Max+ah*Inv*A*eh;
      r1  = Qr(i,i+1)*Qrh(i+1,i);
      R   = r1+a*bh+a*Ah*Inv*A*bh;
      r1h = Qrh(j,j-1)*Qr(j-1,j);
      Rh  = r1h+ah*Inv*b;
      
      % For Part 3 of cycle
      
      A3  = Q(I,J);
      A3h = Qh(J,I);
      a3  = Q(i,J);
      a3h = Qh(j,I);
      e3  = 1 - sum(Q(I,i+2:j-1),2);
      e3h = 1 - sum(Qh(J,i+1:j-2),2);
      
      Inv3 = inv(Eye-A3*A3h);
      d3  = sum(Q(i,i+1:j-1));
      E3 = 1-d3 + a3*A3h*Inv3*e3;
      d3h = sum(Qh(j,i+1:j-1));
      E3h = 1-d3h + a3h*Inv3*A3*e3h;
      
    end
    
    if ~(i == 1 && j == n)
      
      % Standing
      if j == n
        pS1 = 0; pS2=0; pS3=0;
      else
        
        % Part 1 of cycle
        ES = E + C*Dh/(1-Rh);
        RS = R + C*Ch/(1-Rh);
        pS1 = ES/(1-RS);
        
        % Part 2 of cycle
        pS2 = Ch/(1-Rh);
        
        % Part 3 of cycle
        pS3 = E3h;
      end
      
      % Hanging
      if i == 1
        pH1=0; pH2=0; pH3=0;
      else
        
        % Part 1 of cycle
        EH = Eh + Ch*D/(1-R);
        RH = Rh + Ch*C/(1-R);
        pH1 = EH/(1-RH);
        
        % Part 2 of cycle
        pH2 = C/(1-R);
        
        % Part 3 of cycle
        pH3 = E3;
      end
      
      % Check 
      
    end
    
    if Min == 0 || Max == 0
      pS = 0;
      pH = 0;
    elseif ~(i == 1 & j == n)
      pS = Max*pS1*pS2*pS3;
      pH = Min*pH1*pH2*pH3;
    else % i == 1 & j == n
      pS = Min*E3; % = Max*E3h;
      pH = 0;
      %     Old code      
      %      pS = 1/2*Min*E3;
      %      pH = 1/2*Max*E3h;
    end
    
    F_rfc(i,j) = pS; % Standing
    F_rfc(j,i) = pH; % Hanging
    
  end
end
%fprintf(1,'\n');

% Check negative elements

[I,J] = find(F_rfc<0);
if ~isempty(I)
  warning(['Negative elements in F_rfc. Setting to zero!']);
  for k = 1: length(I)
    F_rfc(I(k),J(k)) = 0;
  end
end

% Multiply with  the intensity of local minima for the swithching process

F_rfc = c_m*F_rfc;

T(6,:)=clock;


