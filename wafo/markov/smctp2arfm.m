function [F_rfc,F_rfc_z,T] = smctp2arfm(P,F,c_m,SideInfo)
%SMCTP2ARFM  Calculates the asymmetric rainflow matrix for a SMCTP.
%
% CALL:  [Frfc]        = smctp2arfm(P,F);
%        [Frfc,Frfc_z] = smctp2arfm(P,F,c_m,SideInfo);
%
% Frfc    = Rainflow Matrix                            [n,n]
% Frfc_z  = Rainflow Matrix, with side info z    {r,r1}[n,n]
%
% P       = Transition matrix for regime process.      [r,r]
% F       = Cell array of min-Max and Max-min matrices {r,2}
% F{i,1}  = min-Max matrix, process i                  [n,n]
% F{i,2}  = Max-min matrix, process i                  [n,n]
% c_m     = Intensity of local minima, switching proc. [1,1]
%           (Default: 1)
% SideInfo = Which type of side information
%          0: No side info
%          1: Mark min & max, r1=r
%          2: Mark when counted, r1=1
%
% Calculates the asymmetric rainflow matrix for a switching process
% with a Markov chain of turning points within each regime.
% If a matrix F{i,2}=[], then the process will be assumed to be 
% time-reversible.
%
% Example: (Two processes as in Example 4.1 in PhD thesis)
%   P = [0.9 0.1; 0.05 0.95];
%   param = [-1 1 32]; u = levels(param);
%   F1 = mktestmat(param,[-0.4 -0.3],0.15,1);
%   F2 = mktestmat(param,[0.3 0.4],0.15,1);
%   [Frfc,Frfc_z] = smctp2arfm(P,{F1 F1'; F2 F2'},1,1);
%   figure(1),cmatplot(u,u,Frfc_z,3)
%   Frfc1 = Frfc_z{1,1}+Frfc_z{1,2}+Frfc_z{2,1}+Frfc_z{2,2};
%   figure(2),cmatplot(u,u,{Frfc Frfc1},3) % Shall be identical
%
% See also  mctp2arfm, smctp2rfm, dtp2arfm_sid

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
%error(nargchk(2,4,ni));
narginchk(2,4)
if ni < 3
  c_m=[];
end

if ni < 4
  SideInfo = [];
end

if isempty(c_m)
  c_m=1;
end
if isempty(SideInfo)
  SideInfo=0;
end

% Define 

Zstr = '123456789'; Zstr=Zstr';

r = length(P);   % Number of regime states
n = length(F{1,1});  % Number of levels

% Check that the rowsums of P are equal to 1

P = mat2tmat(P);

T(2,:)=clock;

% Normalize the rowsums of F{1,1},...,F{r,1} to 1
%  ==>  Q{1,1},...,Q{r,1}

for i = 1:r
  QQ{i,1} = F{i,1};
  QQ{i,1} = mat2tmat(QQ{i,1},1);
end

T(3,:)=clock;

% Normalize the rowsums of F{1,2},...,F{r,2} to 1
%  ==>  Q{1,2},...,Q{r,2}

for i = 1:r
  
  if isempty(F{i,2})        % Time-reversible
    QQ{i,2} = F{i,1};
    QQ{i,2} = QQ{i,2}';  
  else                   % Fhi is given
    QQ{i,2} = F{i,2}; 
  end
    
  QQ{i,2} = mat2tmat(QQ{i,2},-1);

end

T(4,:)=clock;

% Make the transition matrix Q for the joint min-Max process

Q = zeros(n*r,n*r);
I = 0:r:(n-1)*r;
for z=1:r
  Q0 = kron(QQ{z,1},P);
  Q(I+z,:) = Q0(I+z,:);
end

T(5,:)=clock;

% Make the transition matrix Qh for the joint Max-min process

Qh = zeros(n*r,n*r);
I = 0:r:(n-1)*r;
for z=1:r
  Q0 = kron(QQ{z,2},P);
  Qh(I+z,:) = Q0(I+z,:);
end

T(6,:)=clock;

% Stationary distribution (=ro) of local minima with transition matrix
% Qt = Q*Qh = "Transition matrix for min-to-min"

Qt = Q*Qh;
ro = mc2stat(Qt(1:r*(n-1),1:r*(n-1)));  % Stationary distr., row vector  
ro = [ro zeros(1,r)];  % Minimum can't reach the highest level

% Stationary distribution (=roh) of local maxima with transition matrix
% Qt = Qh*Q = "Transition matrix for max-to-max"

Qth = Qh*Q;
roh = mc2stat(Qth(r+1:r*(n),r+1:r*(n)));  % Stationary distr., row vector  
roh = [zeros(1,r) roh];  % Maximum can't reach the highest level

T(7,:)=clock;

% Make the frequency matrix FF for the joint min-Max and Max-min
% distribution

FF = Q.*(ro'*ones(1,n*r)) + Qh.*(roh'*ones(1,n*r));

% Create Transition matrices for time-reversed MC

% Backward min-to-max
I1 = find(ro>0); I2 = find(ro<=0);
ro_inv = ro; ro_inv(I1) = 1./ro(I1); ro_inv(I2) = zeros(1,length(I2));
Qr = Qh' .* (ro_inv'*roh);

% Backward max-to-min
I1 = find(roh>0); I2 = find(roh<=0);
roh_inv = roh; roh_inv(I1) = 1./roh(I1); roh_inv(I2) = zeros(1,length(I2));
Qrh = Q' .* (roh_inv'*ro);

% Make the frequency matrix FF for the joint min-Max and Max-min
% distribution

FF1 = Qr.*(ro'*ones(1,n*r)) + Qrh.*(roh'*ones(1,n*r));

T(8,:)=clock;

% Initiation of matrices

F_rfc = zeros(n,n);
EYE = eye(n*r,n*r); Eye1  = eye(r,r);
Zero1 = zeros(r,1); Zero2 = zeros(r,r);
One1 = ones(r,1);   One2  = ones(r,r);

if SideInfo == 1  % Rainflow matrix with side information
  [F_rfc_z{1:r,1:r}] = deal(zeros(n,n));
end % Rainflow matrix with side information
    
% Calculate rainflow matrix 


%fprintf(1,'Calculating row ');
for i=1:n-1
%  fprintf(1,'-%1d',i);

  for j=i+1:n

    I = r*(i-1)+1:r*i;
    J = r*(j-1)+1:r*j;
    I1 = r*(i+1-1)+1:r*(i+1);
    J1 = r*(j-1-1)+1:r*(j-1);
    I1J1 = r*(i+1-1)+1:r*(j-1);

    x = FF(I,J);
    y = FF(J,I);

    Min = sum(FF(I,:),2)';
    Max = sum(FF(J,:),2)';

    Ro = ro(I);    % Probability of "min=i"
    Roh = roh(J);  % Probability of "max=j"
%    fprintf(1,'Min=%f, Max=%f, Ro=%f, Roh=%f\n',Min,Max,Ro,Roh);

    Min = Ro; Max = Roh; % Just to be sure
    
    Min = Min'*One1';
    Max = Max'*One1';
    
    % For all subdiagonals

    C = y'./Min; Ch = x'./Max;
    E = 1-sum(y'./Min,2); Eh = 1-sum(x'./Max,2);

    
    if j-i >= 2  % Second and higher subdiagonals

      % For Part 1 & 2 of cycle

      d1  = Qr(I,I1)*(1-sum(Qrh(I1,I),2));
      D   = d1;
      d1h = Qrh(J,J1)*(1-sum(Qr(J1,J),2));
      Dh  = d1h;
      d0  = sum(Qr(I,I1J1),2); 
      E   = E-d0; 
      d0h = sum(Qrh(J,I1J1),2); % ???
      Eh  = Eh-d0h;
      r1  = Qr(I,I1)*Qrh(I1,I);
      R   = r1;
      r1h = Qrh(J,J1)*Qr(J1,J);
      Rh  = r1h;

      % For Part 3 of cycle

      d3  = sum(Q(I,I1J1),2);
      E3  = 1-d3;
      d3h = sum(Qh(J,I1J1),2);
      E3h = 1-d3h;

    else % First subdiagonal
      
      D  = Zero1; Dh  = Zero1;
      R  = Zero2; Rh  = Zero2;
      E3 = One1;  E3h = One1;
      
    end
    
      
    if j-i >= 3  % Third and higher subdiagonals

      Eye = EYE(1:r*(j-i-2),1:r*(j-i-2));

      % For Part 1 & 2 of cycle

      II  = r*(i+1-1)+1:r*(j-2);  % i+1:j-2
      JJ  = r*(i+2-1)+1:r*(j-1);  % i+2:j-1;
      
      A  = Qr(II,JJ);
      Ah = Qrh(JJ,II);
      Inv = inv(Eye-A*Ah);
      
      a  = Qr(I,JJ);
      ah = Qrh(J,II);
      b  = Qr(II,J);
      bh = Qrh(JJ,I);

      e  = 1 - sum(Qr(II,r*(i+2-1)+1:r*(j)),2);   % i+2:j
      eh = 1 - sum(Qrh(JJ,r*(i-1)+1:r*(j-2)),2);  % i:j-2

      C   = C + a*Ah*Inv*b;
      Ch  = Ch + ah*Inv*A*bh;
      D   = D + a*eh+a*Ah*Inv*A*eh;
      Dh  = Dh + ah*Inv*e;
      E   = E + a*Ah*Inv*e;
      Eh  = Eh + ah*Inv*A*eh;
      R   = R + a*bh+a*Ah*Inv*A*bh;
      Rh  = Rh + ah*Inv*b;

      % For Part 3 of cycle

      A3  = Q(II,JJ);
      A3h = Qh(JJ,II);
      Inv3 = inv(Eye-A3*A3h);
      c3  = Q(I,JJ);
      c3h = Qh(J,II);
      e3  = 1 - sum(Q(II,r*(i+2-1)+1:r*(j-1)),2);   % i+2:j-1
      e3h  = 1 - sum(Qh(JJ,r*(i+1-1)+1:r*(j-2)),2); % i+1:j-2

      E3  = E3  + c3*A3h*Inv3*e3;
      E3h = E3h + c3h*Inv3*A3*e3h;

    end

    if ~(i == 1 && j == n)

      % Standing
      if j == n
        pS1 = Zero1; pS2=Zero2; pS3=Zero1;
      else

        % Part 1 of cycle
	AS = (Eye1-Rh)\Dh;
        ES = E + C*AS;
	BS = (Eye1-Rh)\Ch;
        RS = R + C*BS;
        pS1 = (Eye1-RS)\ES;

        % Part 2 of cycle
        pS2 = (Eye1-Rh)\Ch;

        % Part 3 of cycle
        pS3 = E3h;
	
      end

      % Hanging
      if i == 1
        pH1 = Zero1; pH2=Zero2; pH3=Zero1;
      else

        % Part 1 of cycle
	AH = (Eye1-R)\D;
        EH = Eh + Ch*AH;
	BH = (Eye1-R)\C;
        RH = Rh + Ch*BH;
        pH1 = (Eye1-RH)\EH;

        % Part 2 of cycle
        pH2 = (Eye1-R)\C;

        % Part 3 of cycle
        pH3 = E3;
      end

    else % i == 1 & j == n
      
      % Standing
      pS1 = One1; 
      pS2 = Eye1;
      pS3 = E3h*sum(Roh)/sum(Ro+Roh);
      
      % Hanging
      pH1 = One1; 
      pH2 = Eye1;
      pH3 = E3*sum(Ro)/sum(Ro+Roh);
      
    end

    % Standing RFC
    F_rfc(i,j) = Roh*diag(pS3)*pS2*pS1;
    
    % Hanging RFC
    F_rfc(j,i) = Ro*diag(pH3)*pH2*pH1;
    
    if SideInfo == 1  % Rainflow matrix with side information
      if (i==1) && (j==n)
	
        b3  = Q(II,J);
        EE  = Q(I,J) + c3*A3h*Inv3*b3;
        for z=1:r
  	  for w=1:r
	  
            % Standing RFC
	    F_rfc_z{z,w}(i,j) = 1/2*Ro(z)*EE(z,w);
	  
            % Hanging RFC
	    F_rfc_z{z,w}(j,i) = 1/2*Ro(z)*EE(z,w);
	  end
	end
	
      else
	
        for z=1:r
  	  for w=1:r
	  
            % Standing RFC
	    F_rfc_z{z,w}(i,j) = Roh(w)*pS3(w)*pS2(w,z)*pS1(z);
	  
            % Hanging RFC
	    F_rfc_z{z,w}(j,i) = Ro(w)*pH3(w)*pH2(w,z)*pH1(z);
	  
	  end
	end
	
      end
    elseif SideInfo ~= 0
      error(['SideInfo = ' num2str(SideInfo) ' not implemented']);
    end % Rainflow matrix with side information
    
  end 
end
%fprintf(1,'\n');

% Multiply with the intensity of local minima for the swithching process

F_rfc = c_m*F_rfc;

if SideInfo == 1
  for z=1:r
    for w=1:r
      F_rfc_z{z,w} = c_m*F_rfc_z{z,w}; 
    end
  end
end

T(9,:)=clock;

