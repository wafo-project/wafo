function [hc,hi,ht,u1] = hill(data,k1,k2,Nmin)
% HILL Return the hill estimator for the tail index
%
% 

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



options = [];
sd = -sort(-data(:));

sd = sd(sd>0);
Ns = length(sd);

if nargin<2 || isempty(k1)
  k1 = floor(Ns/4)+1;
end
if nargin<3 || isempty(k2)
  k2 = floor(Ns/2)+1;
end
if nargin<4 || isempty(Nmin)
  Nmin = 0;
end
Nmin = min(Nmin,k1);


k = (1:Ns).';
logsd = log(sd);
clogsd = cumsum(log(sd))./k;
Hk = clogsd-logsd;
Hk(1) = NaN;

options.k = k(Nmin+1:Ns);
options.u = sd(options.k);

options.theta = log(options.k)/log(Ns);
titleTxt = 'Hill plot';
[p,s,mu] = polyfit(k(k1:k2),1./Hk(k1:k2),1);
[iHk,delta] = polyval(p,k(Nmin+1:k1),s,mu);
ei = 1./Hk(Nmin+1:k1)-iHk;
%plot(options.k,ei,options.k,3*[-delta delta],'r')
ix = find(ei<-3*delta | 3*delta<ei,1,'last');
if isempty(ix)
  ix = 1;
end
u1 = sd(ix+Nmin);
k1 = options.k;
tailindex = 1./Hk(k1);
hc = createwdata('data',tailindex,'args',sd(k1),...
'dataCI',[],'title',titleTxt,'labels',{'Threshold','Tail index'},...
    'workspace',options,'note',titleTxt);

hi = createwdata('data',tailindex,'args',k1,...
'title',titleTxt,'labels',{'Order statistics','Tail index'},...
    'workspace',options);

 ht = createwdata('data',tailindex,'args', options.theta  ,...
'title',titleTxt,'labels',{'Theta','Tail index'},...
    'workspace',options);

 
  
if ~isoctave
  hc = wdata(hc);
  hi = wdata(hi);
  ht = wdata(ht);
end
