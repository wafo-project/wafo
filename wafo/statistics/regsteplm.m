function [Q, I, B, BB] = regsteplm(y,x,varargin)
%REGSTEPLM Stepwise predictor subset selection for Linear Model regression
%
%  CALL:  [Q, I, B, BB] = regsteplm(y,X,options)
%
%          Q = Criterion value as a function of the number of parameters
%          I = Index numbers of the included columns in final model
%          B = Coefficients for the final model (Y = X(:,I)*B).
%         BB = Column p of BB is the best B of parameter size p. 
%      y     = Column vector of observed/response values
%      X     = matrix of regressors, with the last column filled with the constant value 1
%    options = struct defining performance of subset selection. 
%      .criterion : The criterion for subset selection. Valid values are:
%	                 'HT'   Hypothesis Test (using alpha)
%	                 'AIC'  Akaike's Information Criterion 
%	                 'BIC'  Bayesian Information Criterion
%	                 'CMV'  Cross Model Validation (inner criterion RSS) (default)
%      .how       : defines how new models are included. Chose between 
%	                 'all'      : All Subsets
%	                 'forward'  : Forward Inclusion
%	                 'backward' : Backward Elimination
%      .pmax      : Maximum number of parameters in model
%      .alpha     : Confidence level (default 0.05)
%
%	  REGSTEPLM selects a good subset of regressors in a multiple linear 
%	  regression model. The returned Q is the criterion as a function of 
%	  the number of parameters; it might be interpreted as an 
%	  estimate of the prediction standard deviation. For the method
%	  'HT' the reported Q is instead the successive p-values for
%	  inclusion or elimination.
%
%	  The last column of the prediction matrix x must be an intercept 
%	  column, i.e., all elements are ones. This column is never excluded 
%	  in the search for a good model. If it is not present it is added.  
%
%	  For the method 'HT' the optional input argument "alpha" is the 
%	  p-value reference used for inclusion or deletion.
% 
%  Example
%  x=(1:10)';  % Covariate
%  y= x+rand(10,1);
%  [Q,ix,B,BB] = regsteplm(y,[x x.^2,cos(x),sin(x),exp(x)]);
%  plot(Q),xlabel('Number of parameters'),ylabel('Q');
%  % chose model with 2 parameters
%  iy = find(BB(:,2));
%  b = reglm(y,x(:,iy(1:end-1)));
%  [y1,ylo,yup] = b.predict(); 
%  plot(x, y,'o',x,y1,x,[ylo,yup],'r','LineWidth',2);
%
%  close all;
%
%  See also reglm and polyfit.


%       Copyright (c) 1994, 2007 Anders Holtsberg, Per A. Brodtkorb
%
% This program is free software; you can redistribute it and/or modify
%     it under the terms of the GNU Lesser General Public License as published by
%     the Free Software Foundation; either version 3 of the License, or
%     (at your option) any later version.
% 
%     This program is distributed in the hope that it will be useful,
%     but WITHOUT ANY WARRANTY; without even the implied warranty of
%     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%     GNU General Public License for more details.
% 
%     You should have received a copy of the GNU Lesser General Public License
%     along with this program.  If not, see <http://www.gnu.org/licenses/>.




%	  This function is not highly optimized for speed but rather for
%	  flexibility. It would be faster if 'all subsets' were in a 
%	  separate routine and 'forward' and 'backward' were in another
%	  routine, especially for CMV.


options = struct('criterion','CMV','how','forward','pmax',[],'alpha',0.05);
if nargin==1 && strcmpi(y,'defaults')
  Q = options;
  return
end
%error(nargchk(1,inf,nargin))
narginchk(1,inf)
options = parseoptions(options,varargin{:});

validHow = {'all','forward','backward'};
ix = find(strncmpi(options.how,validHow,1));
if any(ix)
  options.how = validHow{ix};
else
  error('Option argument error: Invalid ''how'' option!')
end
validCrit = {'HT','AIC','BIC','CMV'};
ix = find(strncmpi(options.criterion,validCrit,1));
if any(ix)
  options.criterion = validCrit{ix};
else
  error('Option argument error: Invalid criterion option!')
end


n = length(y);
nc = size(x,2);

if isempty(options.pmax)
   pmax = nc;
else
   pmax = options.pmax;
end

if strcmpi(options.how,'backward')
   pmax = nc;
end

if any(x(:,nc)~=1)
   fprintf('   An intercept column added')
   x = [x ones(n,1)];
   nc = nc + 1;
   pmax = pmax + 1;
end
if nc<2, disp('only one model'), return, end 


      
Qsml = NaN*ones(pmax,1);

XX = x'*x;
XY = x'*y; 
YY = y'*y;  

% === If all subsets then set up an all-subsets-indicator-matrix ====

if strcmpi(options.how,'all')
   C = [];
   for i = 1:nc-1
      d = max(1,size(C,1));
      C = [zeros(d,1) C; ones(d,1) C];
      H = C * ones(size(C,2),1);
      J = find(H<pmax);
      C = C(J,:);
   end
   H = H(J);
   [S,I] = sort(H);
   C = C(I,:);
   C = [C ones(size(C,1),1)];
   AllSubsets = C;
   AllSubsetsH = sum(C')';
end

% === This is for CMV ===============================================

if strcmp(options.criterion,'CMV')
   dataloopend = n + 1;
   Qcmv = zeros(pmax,1);
   XXs = XX;
   XYs = XY;
   YYs = YY;
else
   dataloopend = 1;
end

for idata = 1:dataloopend
   
   if strcmp(options.criterion,'CMV')
      fprintf('\n %3.0f:', idata*(idata ~= dataloopend))
      if idata == dataloopend
         XX = XXs;
         XY = XYs;
         YY = YYs;
      else
         xi = x(idata,:);
         yi = y(idata);
         XX = XXs - xi'*xi;
         XY = XYs - xi'*yi;
         YY = YYs - yi^2; 
      end
   end

% === Now we begin to loop over model sizes =========================

   if strcmp(options.how,'backward')
      p = nc;
      loopto = 1;
      C = ones(1,nc);
   else
      p = 1;
      loopto = pmax;
      C = [zeros(1,nc-1) 1];
   end
   BB = zeros(nc,pmax);

   fprintf('  ');

   while 1;

% === The whole loop over the models of size p  ======================
     
      fprintf('%2.0f ', p)
      qmin = 1e99;
      for k = 1:size(C,1)
         Jk = find(C(k,:));
         q = YY - XY(Jk)'*(XX(Jk,Jk)\XY(Jk));               
         if q < qmin
            qmin = q;
            Cbest = C(k,:);
         end
      end
      Qsml(p) = qmin;
      I = find(Cbest);
      BB(I,p)  = XX(I,I)\XY(I);
   
% === And a piece for CMV only ======================================

      if strcmp(options.criterion,'CMV') && idata < n + 1
         Jk = find(Cbest);
         q = yi - xi(Jk)*(XX(Jk,Jk)\XY(Jk));               
         Qcmv(p) = Qcmv(p) + q^2;
      end

% === Next parameter size ===========================================
   
      if p == loopto, break, end

      if strcmp(options.how,'forward')
         p = p + 1;
         C = [];
         for i = 1:nc-1
            if Cbest(i)==0
               Cnew = Cbest;
               Cnew(i) = 1;
               C = [C; Cnew];
            end
         end

      elseif strcmp(options.how,'backward')
         p = p - 1;
         C = [];
         for i = 1:nc-1
            if Cbest(i)==1
               Cnew = Cbest;
               Cnew(i) = 0;
               C = [C; Cnew];
            end
         end
   
      elseif strcmp(options.how,'all')
         p = p + 1;
         C = AllSubsets((AllSubsetsH==p),:);
      end
   
   end
end

% === Finished at last ===============================================

switch options.criterion
  
  case 'CMV'
   Q = sqrt(Qcmv/n);
  case 'AIC'
   Q = sqrt(Qsml/n).*exp((1:pmax)'/n);
  case 'BIC'
   Q = sqrt(Qsml/n).*exp((1:pmax)'/n*log(n)/2);
  case 'HT'
   Qsml = sqrt(Qsml/n);
   Fvals = (Qsml(1:pmax-1).^2 - Qsml(2:pmax).^2) ...
            ./ (Qsml(2:pmax).^2 ./ (n-(2:pmax))' );
   Q = 1-cdff(Fvals,1,(n-(2:pmax))');
   %i = find(Q>alpha);
   if strcmpi(options.how,'backward')
      i = [1; Q<alpha];
      i = find(i);
      i = i(length(i));
   else
      i = [Q>alpha; 1];
      i = find(i);
      i = i(1);
   end
end
if ~strcmpi(options.criterion,'HT')
   [S,i] = sort(Q);
   i = i(1);
end

B = BB(:,i);
I = find(B);

fprintf('\n')
