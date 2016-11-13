function [Fest,Est,OPTIONS] = estsmctp(Fobs,whatEst,method,known,whatKnown,init,OPTIONS)
%ESTSMCTP  Estimate SMCTP model from an observed rainflow matrix.
%
%   Estimates parameters in a Switching Markov Chain
%   of Turning Points from an observation of the rainflow matrix.
%
% CALL:  [Fest,Est] = estsmtp(Fobs,whatEst,method,known,whatKnown,init,OPTIONS)
%
% Fest      = Estimated SMCTP model.                     [SA]
% Est       = Estimated parameters.                      [SA]
%
% Fobs      = Observation of rainflow matrix.            [nxn]
% whatEst   = See below.
% method    = 'ML' / 'chi2' / 'HD' / 'KL' (See below.)
% known     = Values of known parameters of the model.   [SA]
% whatKnown = Which parameters are known? (Not used!)    [SA]
% init      = Initial guess. (for optimization)          [SA]
% OPTIONS   = Options to optimization routine. (Optional)
%             (see 'help foptions')
%
% [SA]=[structure array]
%
% method:
%   'ML'   : Approximate Maximum Likelihood (Multinomial)
%   'chi2' : Chi-square distance
%   'HD'   : Hellinger distance
%   'KL'   : Kullback-Leibler distance
%
% whatEst:
%  'P'         : Estimate P-matrix, min-max matrices for the
%     subloads are known [known.F].
%  'MeanStd'   : Estimate mean and std of subloads.
%     The shape of min-max matrices for the subloads are known
%     [known.F]. P-matrix known [known.P].
%  'P,MeanStd' : Also estimate P-matrix otherwise as above.
%  'CalcF'     :
%  'P,CalcF'   :
%  'SimF'      :
%  'P,SimF'    :
%
% Side Information: known.SideInfo = 
%   0:  No side information 
%   11: Mark min&max, y = 'regime process'
%   12: Mark min&max, y = 'scrambled regime process'
%   21: Mark when counted, y = 'regime process'
%   22: Mark when counted, y = 'scrambled regime process'
%   (Optional, Default = 0, No side information)
%
% known.NOsubzero = Number of subdiagonals that are zero
%   (Optional, Default = 0, only the diagonal is zero)
%
% Example:
%   M1.x0=[-0.4 -0.3]; M1.s=0.15; M1.lam=1; 
%   M2.x0=[0.3 0.4]; M2.s=0.15; M2.lam=1;
%   F1 = mktestmat([-1 1 8],M1.x0,M1.s,M1.lam);
%   F2 = mktestmat([-1 1 8],M2.x0,M2.s,M2.lam);
%   P=[1-0.1 0.1; 0.05 1-0.05];             % Transition matrix
%   [xD,z] = smctpsim(P,{F1 []; F2 []},5000); % Simulate
%   Fobs = dtp2rfm(xD,8);
%
%   known.F = {F1 []; F2 []};   % known min-max and max-min matrices
%   init.P = P;                 % initial guess of P-matrix
%   [Fest,Est] = estsmctp(Fobs,'P','ML',known,[],init);
%
%   known.Ffun = 'f_funm';      % Function for calculating a submodel
%   known.trModel2X = 'tr_m2x'; % transform from Model to X-vector
%   known.trX2Model = 'tr_x2m'; % transform from X-vector to model
%   known.param = [-1 1 8];
%   init.P = P;       % initial guess of P-matrix
%   init.M = {M1 M2}; % initial guess of Models for min-max mat
%   [Fest,Est] = estsmctp(Fobs,'P,CalcF','ML',known,[],init);
%
%  % Further examples of using ESTSMCP can be found in WDEMOS/ITMKURS,
%  % especially in the script ITMKURS_LAB2.
%
% See also

% Updated by PJ 07-Apr-2005
%   Adaptation for Matlab 7.
%   Changed 'fmins' to 'fminsearch'.

% Check input aruments

ni = nargin;
no = nargout;
error(nargchk(6,7,ni));

if ni < 7
  OPTIONS = struct();
end

if ~isfield(known,'NOsubzero')
  known.NOsubzero = 0;
end

if ~isfield(known,'SideInfo')
  known.SideInfo = 0;
end

% Options to fmins

%  OPTIONS(1)  = 0;       %  1 = intermediate steps in the solution are displayed
%  OPTIONS(2)  = 1e-2;    % the termination tolerance for x;
%  OPTIONS(3)  = 1e-2;    % the termination tolerance for F(x);
%  OPTIONS(2)  = 1e-1;    % the termination tolerance for x;
%  OPTIONS(3)  = 1e-1;    % the termination tolerance for F(x);
if 0
  % The HTML documentation software is unable to track dependencies to
  % functions evaluated with feval. Here is a
  % trick to get the html documentation right (pab 2005)
  M.s = 1; 
  X = tr_m2x(M);
  M = tr_x2m(X,known);
end

switch whatEst

case 'P'      % Estimate P

  P=init.P;
  X0 = tr_p2x(P,1);

  try 
    % For Matlab 5.3 and higher ???
    [X,OPTIONS] = fminsearch(@(x)f_smctp(x, Fobs,whatEst,method,known,whatKnown,init),X0,OPTIONS);
  catch
    % For Matlab 5.2 and lower ???
    [X,OPTIONS] = fmins(@(x)f_smctp(x, Fobs,whatEst,method,known,whatKnown,init),X0,OPTIONS);
  end

  Pest = tr_x2p(X,1);
  Est.P=Pest;

case 'MeanStd' % Estimate Mean and Std

  init.MeanStd(:,2) = log(init.MeanStd(:,2));
  X0 = init.MeanStd(:);

  try 
    % For Matlab 5.3 and higher ???
    [X,OPTIONS] = fminsearch(@(x)f_smctp(x, Fobs,whatEst,method,known,whatKnown,init),X0,OPTIONS);
  catch
    % For Matlab 5.2 and lower ???
    [X,OPTIONS] = fmins(@(x)f_smctp(x, Fobs,whatEst,method,known,whatKnown,init),X0,OPTIONS);
  end

  r = length(X)/2;
  MeanStd = reshape(X,r,2);
  MeanStd(:,2) = exp(MeanStd(:,2));
  Est.MeanStd = MeanStd;

case 'P,MeanStd' % Estimate P, Mean and Std

  P=init.P;
  X1 = tr_p2x(P,1);
  init.MeanStd(:,2) = log(init.MeanStd(:,2));
  X2 = init.MeanStd(:);
  X0 = [X1; X2];

  try 
    % For Matlab 5.3 and higher ???
    [X,OPTIONS] = fminsearch(@(x)f_smctp(x, Fobs,whatEst,method,known,whatKnown,init),X0,OPTIONS);
  catch
    % For Matlab 5.2 and lower ???
    [X,OPTIONS] = fmins(@(x)f_smctp(x, Fobs,whatEst,method,known,whatKnown,init),X0,OPTIONS);
  end

  r=(-1+sqrt(1+4*length(X)))/2;
  X1 = X(1:r*(r-1));
  X2 = X(r*(r-1)+1:end);

  Pest = tr_x2p(X1,1);
  Est.P=Pest;

  MeanStd = reshape(X2,r,2);
  MeanStd(:,2) = exp(MeanStd(:,2));
  Est.MeanStd = MeanStd;

%  Fest = smctp(Pest,known.F);


case {'CalcF','P,CalcF'} % Estimate P, Model parameters

  % transform model to vector
  if whatEst(1) == 'P' % 'P,CalcF'
    P=init.P;
    X0 = tr_p2x(P,1);
    r = length(init.P);
  else
    X0 = [];
    r = length(known.P);
  end

  for i = 1:r
    X = feval(known.trModel2X,init.M{i});
    nM(i) = length(X);
    X0 = [X0; X];
  end

  known.nM = nM;
    
  try 
    % For Matlab 5.3 and higher ???
    [X,OPTIONS] = fminsearch(@(x)f_smctp(x, Fobs,whatEst,method,known,whatKnown,init),X0,OPTIONS);
  catch
    % For Matlab 5.2 and lower ???
    [X,OPTIONS] = fmins(@(x)f_smctp(x, Fobs,whatEst,method,known,whatKnown,init),X0,OPTIONS);
  end

  % transform vector to model
  if whatEst(1) == 'P' % 'P,CalcF'
    r = length(init.P);
    X1 = X(1:r*(r-1));
    X2 = X(r*(r-1)+1:end);
    P = tr_x2p(X1,1);
    Est.P = P;
  else
    P = known.P;
    X2 = X;
  end
  r = length(P);

  % transform vector to model
  k1=1;
  for i = 1:r
    k2 = k1+nM(i)-1;
    M{i} = feval(known.trX2Model,X2(k1:k2),known);
    k1=k2+1;
  end

  Est.M = M;
  

case {'SimF','P,SimF'} % Estimate P, Model parameters

  r = length(known.P);

  % 1. Initial estimate

  M = init.M;

  % 2. Simulate each subload

  F = cell(2,1);
  for i = 1:r
    F{i} = feval(known.simFun,known.param,M{i},known.T,known.T0)/(known.T/2);
  end

%  while ~slut

    % 3. Uppdatera skattning

    % transform model to vector
    X0 = [];
    for i = 1:r
      X = feval(known.trModel2X,init.M{i});
      nM(i) = length(X);
      X0 = [X0; X];
    end

    known.nM = nM;
    known.F = F;

    % Estimate
    try 
      % For Matlab 5.3 and higher ???
      [X,OPTIONS] = fminsearch(@(x)f_smctp(x, Fobs,whatEst,method,known,whatKnown,init),X0,OPTIONS);
    catch
      % For Matlab 5.2 and lower ???
      [X,OPTIONS] = fmins(@(x)f_smctp(x, Fobs,whatEst,method,known,whatKnown,init),X0,OPTIONS);
    end

    % transform vector to model
    k1=1;
    for i = 1:r
      k2 = k1+nM(i)-1;
      M{i} = feval(known.trX2Model,X(k1:k2),known);
      k1=k2+1;
    end

    % 4. Simulate each subload and update

    for i = 1:r
      F1 = feval(known.simFun,M{i},known.T,known.T0)/(known.T/2);
      F{i} = known.theta*F{i} + (1-known.theta)*F1;
    end

    Est.M = M;
%  end

otherwise

  % This should not happen
  error(['Unexpected whatEst: ' whatEst '.'])

end % switch

% Calculate Estimated SMCTP-model


[y,F,P,FF] = f_smctp(X,Fobs,whatEst,method,known,whatKnown,init);

%Fest = smctp(P,FF);
Fest.P = P;
Fest.F = FF;

fprintf(1,'\n');