function [y,F,P,FF] = f_smctp(X,Fobs,whatEst,method,known,whatKnown,init)
%F_SMCTP  Auxiliary function used by ESTSMCTP
%
% CALL:  [y,F,P,FF] = Fsmtp(X,Fobs,whatEst,method,known,whatKnown,init)

%fprintf(1,'X = [%f',X(1));
%fprintf(1,', %f',X(2:length(X)));
%fprintf(1,']\n');


switch whatEst

case 'P'      % Estimate P

  P = tr_x2p(X,1);
  r = length(P);

  FF = known.F;

case 'MeanStd' % Estimate Mean and Std

  P = known.P;

  r = length(X)/2;
  MeanStd = reshape(X,r,2);
  MeanStd(:,2) = exp(MeanStd(:,2));

%  subplot(2,2,1), plotcmat(F1,3)
%  subplot(2,2,2), plotcmat(F2,3)

  n = length(known.F{1,1});
  param = [-1 1 n];
  FF = cell(r,2);
  for i = 1:r
    FF{i,1} = scalemat(param,known.F{i,1},MeanStd(i,1),MeanStd(i,2),param);
  end

%  subplot(2,2,3), plotcmat(F1,3)
%  subplot(2,2,4), plotcmat(F2,3)
%  drawnow

case 'P,MeanStd' % Estimate P, Mean and Std

  r=(-1+sqrt(1+4*length(X)))/2;
  X1 = X(1:r*(r-1));
  X2 = X(r*(r-1)+1:end);

  P = tr_x2p(X1,1);

  MeanStd = reshape(X2,r,2);
  MeanStd(:,2) = exp(MeanStd(:,2));

  n= length(known.F{1,1});
  param = [-1 1 n];
  FF = cell(r,2);
  for i = 1:r
    FF{i,1} = scalemat(param,known.F{i,1},MeanStd(i,1),MeanStd(i,2),param);
  end


case {'P,CalcF','CalcF'} % Estimate Model parameters

  if whatEst(1) == 'P' % 'P,CalcF'
    r = length(init.P);
    X1 = X(1:r*(r-1));
    X = X(r*(r-1)+1:end);
    P = tr_x2p(X1,1);
  else
    P = known.P;
  end
  r = length(P);

  % transform vector to model
  FF = cell(r,2);
  k1=1;
  for i = 1:r
    k2 = k1+known.nM(i)-1;
    M = feval(known.trX2Model,X(k1:k2),known);
    FF{i,1} = feval(known.Ffun,known.param,M);
    k1=k2+1;
  end

case {'SimF'} % Estimate Model parameters

  P = known.P;
  r = length(P);

  % transform vector to model
  FF = cell(r,2);
  k1=1;
  for i = 1:r
    k2 = k1+known.nM(i)-1;
    M = feval(known.trX2Model,X(k1:k2),known);
    FF{i,1} = feval(known.simFun,known.param,M,known.T,known.T0);
    k1=k2+1;
  end


otherwise

  % This should not happen
  error(['Unexpected whatEst: ' whatEst '.'])

end

if known.NOsubzero ~= 0
  for i = 1:r % Set to zero below and on the sub-diagonals
    FF{i,1} = flipud(triu(flipud(FF{i,1})',1+known.NOsubzero)');
  end
end

if known.SideInfo == 0
  F = smctp2rfm(P,FF);
elseif known.SideInfo == 11
  [F,Fsid] = smctp2arfm(P,FF,1,1);
  n = length(Fsid{1,1});
  F = zeros(n*r,n*r);
  FobsSid = Fobs;
  Fobs = zeros(n*r,n*r);
  for z = 1:r
    for w = 1:r
      Fsid{z,w}(1,1) = Fsid{z,w}(1,1) + Fsid{z,w}(n,n); 
      Fsid{z,w}(n,n) = 0; 
      FobsSid{z,w}(1,1) = FobsSid{z,w}(1,1) + FobsSid{z,w}(n,n); 
      FobsSid{z,w}(n,n) = 0; 
      F(1+n*(z-1):n*z,1+n*(w-1):n*w) = Fsid{z,w};
      Fobs(1+n*(z-1):n*z,1+n*(w-1):n*w) = FobsSid{z,w};
    end
  end
else
  error(['Method for SideInfo = ' num2str(known.SideInfo) ' not implemented']);
end

%F = F * sum(sum(Fobs));

if strcmp(method,'ML') == 1
  y = -loglcmat(Fobs,F);
elseif strcmp(method,'chi2') == 1
  y = chi2cmat(Fobs,F);
elseif strcmp(method,'HD') == 1
  y = hdcmat(Fobs,F);
elseif strcmp(method,'KL') == 1
  y = klcmat(Fobs,F);
else
  fprintf(1,['Method ' method ' not implemented']);
end

%fprintf(1,'y=%f\n',y);
fprintf(1,'*');
