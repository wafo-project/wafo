function [low,up] = findciproflog(Lp,alpha)
%FINDCIPROFLOG Find Confidence Interval from proflog function
%
% CALL [low,up] = findciproflog(Lp,alpha)
%
% Example
%  R = rndweib(1,3,0,100,1);
%  phat = fitweib(R);
%  [x,xlo,xup] = invweib(1/990,phat,'lowertail',false);
%  [Lp0,CI0] = proflog(phat,'i',2,'x',x,'link',@lnkweib,'plot',1);
%  alpha1 = 0.1:0.05:0.5;
%  [xlow,xup]= findciproflog(Lp0,alpha1);
%
% See also proglog

%
%     This program is free software; you can redistribute it and/or modify
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


  LLmax = max(Lp.data);
  LLrange     = 0.5 * invchi2(alpha, 1,'lowertail',false);
  cross_level = LLmax - LLrange;
  n = numel(alpha);
  low = zeros(size(alpha));
  up = low;
  for ix = 1:n
    CI = find1ci(Lp,cross_level(ix));
    low(ix) = CI(1);
    up(ix)  = CI(2);
  end
end  % main

function CI = find1ci(Lp,cross_level)

  ind = findcross(Lp.data,cross_level);
  
  switch length(ind)
  case 0
    CI  = [Lp.workspace.pmin, Lp.workspace.pmax];
    xlabtxt = strtok(Lp.labels{1},'(');
    warning('WAFO:FINDCIPROFLOG','Upper bound for %s is larger!',xlabtxt);
    warning('WAFO:FINDCIPROFLOG','Lower bound for %s is smaller!',xlabtxt);
  case 1
    x0 = ecross(Lp.args,Lp.data,ind,cross_level);
    xlabtxt = strtok(Lp.labels{1},'(');
    isUpcrossing = Lp.data(ind)<Lp.data(ind+1);
    if isUpcrossing
      CI = [x0 Lp.workspace.pmax];
      warning('WAFO:FINDCIPROFLOG','Upper bound for %s is larger!',xlabtxt);
    else
      CI = [Lp.workspace.pmin x0];
      warning('WAFO:FINDCIPROFLOG','Lower bound for %s is smaller!',xlabtxt);
    end
  case 2
    CI = ecross(Lp.args,Lp.data,ind,cross_level);
  otherwise
    warning('Number of crossings too large')
    CI = ecross(Lp.args,Lp.data,ind([1,end]),cross_level);
  end  
end  % find1ci