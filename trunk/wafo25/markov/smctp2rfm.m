function [F_rfc,mu_rfc,T] = smctp2rfm(P,F,c_m)
% SMCTP2RFM  Calculates the rainflow matrix for a SMCTP.
%
% CALL: [Frfc,mu_rfc] = smctp2rfm(P,F,c_m);
%
% Frfc   = Rainflow Matrix                            [n,n]
% mu_rfc = Rainflow Counting intensity                [n,n]
%
% P      = Transition matrix for regime process.      [r,r]
% F      = Cell array of min-Max and Max-min matrices {r,2}
% F{i,1} = min-Max matrix, process i                  [n,n]
% F{i,2} = Max-min matrix, process i                  [n,n]
% c_m    = Intensity of local minima, switching proc. [1,1]
%          (Default: 1)
% If a matrix F{i,2}=[], then the process will
% be assumed to be time-reversible.
%
% Calculates the rainflow matrix for a switching process
%   with a Markov chain of turning points within each regime.
%
% Example: (Example 4.1 in PhD thesis)
%   P = [0.9 0.1; 0.05 0.95];
%   param = [-1 1 32]; u = levels(param);
%   F1 = mktestmat(param,[-0.4 -0.3],0.15,1);
%   F2 = mktestmat(param,[0.3 0.4],0.15,1);
%   Frfc = smctp2rfm(P,{F1 F1'; F2 F2'});
%   cmatplot(u,u,Frfc)
%
% See also  mctp2rfm, smctp2arfm
  
% References
%  
%  P. Johannesson (1999):
%  Rainflow Analysis of Switching Markov Loads.
%  PhD thesis, Mathematical Statistics, Centre for Mathematical Sciences,
%  Lund Institute of Technology.
%  
%  P. Johannesson (1998):
%  Rainflow Cycles for Switching Processes with Markov Structure.
%  Probability in the Engineering and Informational Sciences, 
%  Vol. 12, No. 2, pp. 143-175.
  
% Tested  on Matlab  5.3
%
% History:
% Revised by PJ  23-Nov-1999
%   updated for WAFO
% Created by PJ (Pär Johannesson) 1997
%   Copyright (c) 1997-1998 by Pär Johannesson
%   Toolbox: Rainflow Cycles for Switching Processes V.1.1, 22-Jan-1998


% Check input arguments

ni = nargin;
no = nargout;
error(nargchk(2,3,ni));

if ni < 3
  c_m=[];
end

if isempty(c_m)
  c_m=1;
end

% Define 

Zstr = '123456789'; Zstr=Zstr';

T(1,:)=clock;

r = length(P);   % Number of regime states

n = length(F{1,1});  % Number of levels

% Check that the rowsums of P are equal to 1

sumP = sum(P');
if sum(sumP == 1) ~= r
  warning(['Rowsums of P not equal to 1. Renormalizing!']);
  for i = 1:length(P)
    P(i,:) = P(i,:)/sumP(i);
  end
end

T(2,:)=clock;

% Normalize the rowsums of F{1,1},...,F{r,1} to 1
%  ==>  Q{1,1},...,Q{r,1}

for i = 1:r
  QQ{i,1} = triu(F{i,1},1); % Set zeros below diagonal and diagonal
  % Set negative elements to zero
  [I,J] = find(QQ{i,1}<0);
  if length(I) ~= 0
    warning(['Negative elements in Q' Zstr(i) '. Setting to zero!']);
    for k = 1:length(I)
      QQ{i,1}(I(k),J(k)) = 0;
    end
  end

  sumQQi = sum(QQ{i,1}');
  % Normalize rowsums
  if sum(sumQQi == 1) ~= length(QQ{i,1})
    %disp(['Warning: Rowsums of Q' Zstr(i) ' not equal to 1. Renormalizing!']);
    for j = 1:n-1
      if sumQQi(j)~=0, QQ{i,1}(j,:) = QQ{i,1}(j,:)/sumQQi(j); end
    end
  end
end

T(3,:)=clock;

% Normalize the rowsums of F{1,2},...,F{r,2} to 1
%  ==>  Q{1,2},...,Q{r,2}

% Normalize the rowsums of Fh1,...,Fhr to 1  ==>  Q1,...,Qr

for i = 1:r
  if isempty(F{i,2})        % Time-reversible
    QQ{i,2} = F{i,1}';
  else                   % Fhi is given
    QQ{i,2} = F{i,2}; 
  end
    
  QQ{i,2} = tril(QQ{i,2},-1); % Set zeros above diagonal and diagonal
  % Set negative elements to zero
  [I,J] = find(QQ{i,2}<0);
  if length(I) ~= 0
    warning(['Negative elements in Qh' Zstr(i) '. Setting to zero!']);
    for k = 1:length(I)
      QQ{i,2}(I(k),J(k)) = 0;
    end
  end

  sumQQi = sum(QQ{i,2}');
  if sum(sumQQi == 1) ~= length(QQ{i,2})
    %disp(['Warning: Rowsums of Qh' Zstr(i) ' not equal to 1. Renormalizing!']);
    for j = 2:n
      if sumQQi(j)~=0, QQ{i,2}(j,:) = QQ{i,2}(j,:)/sumQQi(j); end
    end
  end
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

T(7,:)=clock;

% Calculate rainflow matrix 

mu_rfc = zeros(n,n);
EYE = eye(n*r,n*r);
Qcumsum = fliplr(cumsum(fliplr(Q)')');

%fprintf(1,'Calculating row ');
for i=2:n
  %fprintf(1,'-%1d',i);
  for j=i-1:min([i n-1])
    q = Qcumsum(:,r*j+1);

    d = q(1:r*(i-1));                           % 1:i-1
    Ro = ro(1:r*(i-1));                         % 1:i-1

    mu_rfc(i,j) = Ro*d;
  end
  for j=i+1:n-1
    q = Qcumsum(:,r*j+1);

    A = Q(r*(i-1)+1:r*(j-1),r*(i+1-1)+1:r*j);   % i:j-1, i+1:j
    Ah = Qh(r*(i+1-1)+1:r*j,r*(i-1)+1:r*(j-1)); % i+1:j, i:j-1
    C = Q(1:r*(i-1),r*(i+1-1)+1:r*j);           % 1:i-1, i+1:j
    Eye = EYE(1:r*(j-i),1:r*(j-i));             % eye (size of A*Ah)
    d = q(1:r*(i-1));                           % 1:i-1
    e = q(r*(i-1)+1:r*(j-1));                   % i:j-1
    Ro = ro(1:r*(i-1));                         % 1:i-1

    mu_rfc(i,j) = Ro*(d+C*Ah*inv(Eye-A*Ah)*e);
  end
end
%fprintf(1,'\n');


% Calculation of the intensity of local minima for the swithching process

mu_rfc = c_m*mu_rfc;

T(8,:)=clock;

% Convert: mu_rfc  -->  F_rfc

F_rfc = zeros(n,n);
for i= 1:n-1
  for j= i+1:n
    F_rfc(i,j) = mu_rfc(i+1,j-1) - mu_rfc(i,j-1) - mu_rfc(i+1,j) + mu_rfc(i,j);
  end
end

I = F_rfc<0;
if sum(sum(I)) ~= 0
  warning(['Negative elements in calculated rainflow matrix F_rfc. Setting to zero!']);
  F_rfc(I) = 0;
end


if no>1
  mu_rfc = cmat2nt(F_rfc); % Calculate rainflow counting intensity
end


T(9,:)=clock;

