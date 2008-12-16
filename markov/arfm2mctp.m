function [F,T] = arfm2mctp(Frfc)
%ARFM2MCTP  Calculates the markov matrix given an asymmetric rainflow matrix. 
%
% CALL:  F = arfm2mctp(Frfc);
%
% F      = Markov matrix (from-to-matrix)             [n,n]
% Frfc   = Rainflow Matrix                            [n,n]
%
% Example: 
%   param = [-1 1 32]; u = levels(param);
%   F = mktestmat(param,[-0.2 0.2],0.15,2);
%   F = F/sum(sum(F));
%   Farfc = mctp2arfm({F []});
%   F1 = arfm2mctp(Farfc);
%   cmatplot(u,u,{F+F' F1},3);
%   sum(sum(abs(F1-(F+F')))) % should be zero
%
% See also  rfm2mctp, mctp2arfm, smctp2arfm, cmatplot

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
% Created by PJ (Pär Johannesson) 1998
% Copyright (c) 1997-1998 by Pär Johannesson
% Toolbox: Rainflow Cycles for Switching Processes V.1.1, 22-Jan-1998

% Recursive formulation a'la Igor
%
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

ni = nargin;
no = nargout;
error(nargchk(1,1,ni));

T(1,:)=clock;

N = sum(sum(Frfc));
Frfc = Frfc/N;

n = length(Frfc);  % Number of levels

T(7,:)=clock;

% Transition matrices for MC

Q = zeros(n,n);
Qh = zeros(n,n);

% Transition matrices for time-reversed MC

Qr = zeros(n,n);
Qrh = zeros(n,n);

% Probability of minimum and of maximun

MIN = sum(triu(Frfc)') + sum(tril(Frfc));
MAX = sum(triu(Frfc))  + sum(tril(Frfc)');

% Calculate rainflow matrix 

F = zeros(n,n);
EYE = eye(n,n);

%fprintf(1,'Calculating row ');
for k=1:n-1  % k = subdiagonal
%  fprintf(1,'-%1d',i);

  for i=1:n-k  % i = minimum

    j = i+k; % maximum;

    pS = Frfc(i,j);  % Standing cycle
    pH = Frfc(j,i);  % Hanging cycle

    Min = MIN(i);
    Max = MAX(j);

%   fprintf(1,'Min=%f, Max=%f\n',Min,Max);


    if j-i == 2  % Second subdiagonal

      % For Part 1 & 2 of cycle

      %C   = y/Min;
      c0  = 0;
      c1  = 1/Min;
      %Ch  = x/Max;
      c0h = 0;
      c1h = 1/Max;
      d1  = Qr(i,i+1)*(1-Qrh(i+1,i));
      D   = d1;
      d1h = Qrh(j,j-1)*(1-Qr(j-1,j));
      Dh  = d1h;
      d0  = sum(Qr(i,i+1:j-1));
      %E   = 1-d0-y/Min;
      e0  = 1-d0;
      e1  = -1/Min;
      d0h = sum(Qrh(j,i+1:j-1));
      %Eh  = 1-d0h-x/Max;
      e0h = 1-d0h;
      e1h = -1/Max;
      r1  = Qr(i,i+1)*Qrh(i+1,i);
      R   = r1;
      r1h = Qrh(j,j-1)*Qr(j-1,j);
      Rh  = r1h;

      % For Part 3 of cycle

      d3h = sum(Qh(j,i+1:j-1));
      E3h = 1-d3h;
      d3  = sum(Q(i,i+1:j-1));
      E3 = 1-d3;

      % Define coeficients for equation system
      a0 = -pS+2*pS*Rh-pS*Rh^2+pS*R-2*pS*Rh*R+pS*Rh^2*R;
      a1 = -E3h*Max*c1h*e0*Rh+E3h*Max*c1h*e0;
      a3 = -E3h*Max*c1h*e1*Rh+E3h*Max*c1h*Dh*c1+E3h*Max*c1h*e1+pS*c1h*c1-pS*c1h*c1*Rh;

      b0 = -pH+2*pH*R+pH*Rh-2*pH*Rh*R-pH*R^2+pH*Rh*R^2;
      b2 = -Min*E3*e0h*R*c1+Min*E3*e0h*c1;
      b3 = Min*E3*e1h*c1+Min*E3*D*c1h*c1-pH*c1h*c1*R-Min*E3*e1h*R*c1+pH*c1h*c1;

      C2 = a3*b2;
      C1 = (-a0*b3+a1*b2+a3*b0);
      C0 = a1*b0;
      % Solve: C2*z^2 + C1*z + C0 = 0
      z1 = -C1/2/C2 + sqrt((C1/2/C2)^2-C0/C2);
      z2 = -C1/2/C2 - sqrt((C1/2/C2)^2-C0/C2);

      % Solution 1
      x1 = -(b0+b2*z1)/(b3*z1);
      y1 = z1;
      % Solution 2
      x2 = -(b0+b2*z2)/(b3*z2);
      y2 = z2;

      x = x2;
      y = y2;

%      fprintf(1,'2nd: i=%d, j=%d: x1=%f, y1=%f, x2=%f, y2=%f\n',i,j,x1,y1,x2,y2);

      % Test Standing cycle: assume x=y

      C0 = a0; C1 = a1; C2 = a3;
      z1S = -C1/2/C2 + sqrt((C1/2/C2)^2-C0/C2);
      z2S = -C1/2/C2 - sqrt((C1/2/C2)^2-C0/C2);

      % Test Hanging cycle: assume x=y

      C0 = b0; C1 = b2; C2 = b3;
      z1H = -C1/2/C2 + sqrt((C1/2/C2)^2-C0/C2);
      z2H = -C1/2/C2 - sqrt((C1/2/C2)^2-C0/C2);

%      fprintf(1,'2nd: i=%d, j=%d: z1S=%f,: z2S=%f, z1H=%f, z2H=%f\n',i,j,z1S,z2S,z1H,z2H);

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
      %C   = y/Min + a*Ah*Inv*b;
      c0  = a*Ah*Inv*b;
      c1  = 1/Min;
      %Ch  = x/Max + ah*Inv*A*bh;
      c0h = ah*Inv*A*bh;
      c1h = 1/Max;
      d1  = Qr(i,i+1)*(1-Qrh(i+1,i));
      D   = d1+a*eh+a*Ah*Inv*A*eh;
      d1h = Qrh(j,j-1)*(1-Qr(j-1,j));
      Dh  = d1h+ah*Inv*e;
      d0  = sum(Qr(i,i+1:j-1));
      %E   = 1-d0-y/Min+a*Ah*Inv*e;
      e0  = 1-d0+a*Ah*Inv*e;
      e1  = -1/Min;
      d0h = sum(Qrh(j,i+1:j-1));
      %Eh  = 1-d0h-x/Max+ah*Inv*A*eh;
      e0h = 1-d0h+ah*Inv*A*eh;
      e1h = -1/Max;
      r1  = Qr(i,i+1)*Qrh(i+1,i);
      R   = r1+a*bh+a*Ah*Inv*A*bh;
      r1h = Qrh(j,j-1)*Qr(j-1,j);
      Rh  = r1h+ah*Inv*b;

      % For Part 3 of cycle

      A3  = Q(I,J);
      A3h = Qh(J,I);
      Inv3 = inv(Eye-A3*A3h);

      % For Standing cycle
      d3h = sum(Qh(j,i+1:j-1));
      c3h = Qh(j,I);
      e3h = 1 - sum(Qh(J,i+1:j-2),2);
      E3h = 1-d3h + c3h*Inv3*A3*e3h;

      % For Hanging cycle
      d3  = sum(Q(i,i+1:j-1));
      c3  = Q(i,J);
      e3  = 1 - sum(Q(I,i+2:j-1),2);
      E3  = 1-d3 + c3*A3h*Inv3*e3;

    end

    if j-i == 1  % First subdiagonal

        if i == 1
        x = Max;
        y = Max;
      elseif j == n
        x = Min;
        y = Min;
      else
        if pS == 0
          x = 0;
          y = pH;
        elseif pH == 0
          x = pS;
          y = 0;
        else
          x = Min*pS/(Min-pH);
          y = Max*pH/(Max-pS);
        end
      end

    elseif j-i >= 2
      if i == 1
        x = Max*(1-sum(Qh(j,2:j-1)));
        y = Max*(1-sum(Qrh(j,2:j-1)));
      elseif j == n
        x = Min*(1-sum(Q(i,i+1:n-1)));
        y = Min*(1-sum(Qr(i,i+1:n-1)));
      else
        if pS == 0
          x = 0;
          y = pH;
        elseif pH == 0
          x = pS;
          y = 0;
        else
          % Define coeficients for equation system
          a0 = pS*c0h*c0+pS*Rh^2*R-2*pS*Rh*R-E3h*Max*c0h*e0*Rh+E3h*Max*c0h*e0+2*pS*Rh+pS*R-pS*c0h*c0*Rh-pS-pS*Rh^2+E3h*Max*c0h*Dh*c0;
          a1 = pS*c1h*c0+E3h*Max*c1h*Dh*c0-E3h*Max*c1h*e0*Rh-pS*c1h*c0*Rh+E3h*Max*c1h*e0;
          a2 = pS*c0h*c1+E3h*Max*c0h*e1-pS*c0h*c1*Rh+E3h*Max*c0h*Dh*c1-E3h*Max*c0h*e1*Rh;
          a3 = -E3h*Max*c1h*e1*Rh+E3h*Max*c1h*Dh*c1+E3h*Max*c1h*e1+pS*c1h*c1-pS*c1h*c1*Rh;

          b0 = pH*c0h*c0+pH*Rh*R^2-pH+pH*Rh-2*pH*Rh*R-pH*c0h*c0*R+Min*E3*e0h*c0-Min*E3*e0h*R*c0+Min*E3*D*c0h*c0+2*pH*R-pH*R^2;
          b1 = Min*E3*D*c1h*c0+Min*E3*e1h*c0+pH*c1h*c0-Min*E3*e1h*R*c0-pH*c1h*c0*R;
          b2 = -pH*c0h*c1*R-Min*E3*e0h*R*c1+Min*E3*D*c0h*c1+Min*E3*e0h*c1+pH*c0h*c1;
          b3 = Min*E3*e1h*c1+Min*E3*D*c1h*c1-pH*c1h*c1*R-Min*E3*e1h*R*c1+pH*c1h*c1;

          C2 = a2*b3-a3*b2;
          C1 = a0*b3-a1*b2+a2*b1-a3*b0;
          C0 = a0*b1-a1*b0;
%fprintf(1,'i=%d, j=%d, C0/C2=%f,C1/C2=%f,C2=%f\n',i,j,C0/C2,C1/C2,C2);
          % Solve: C2*z^2 + C1*z + C0 = 0
          z1 = -C1/2/C2 + sqrt((C1/2/C2)^2-C0/C2);
          z2 = -C1/2/C2 - sqrt((C1/2/C2)^2-C0/C2);

          % Solution 1
          x1 = -(b0+b2*z1)/(b1+b3*z1);
          y1 = z1;
          % Solution 2
          x2 = -(b0+b2*z2)/(b1+b3*z2);
          y2 = z2;

          x = x2;
          y = y2;

%          fprintf(1,'End: i=%d, j=%d: x1=%f, y1=%f, x2=%f, y2=%f\n',i,j,x1,y1,x2,y2);
        end
      end
    end

%    fprintf(1,'i=%d, j=%d: x=%f, y=%f\n',i,j,x,y);

    % min-max
    F(i,j) = x;

    % max-min
    F(j,i) = y;

    % Fill the transitions matrices
    Q(i,j)   = x/Min;
    Qh(j,i)  = y/Max;
    Qr(i,j)  = y/Min;
    Qrh(j,i) = x/Max;
 
  end
end
%fprintf(1,'\n');


T(8,:)=clock;




