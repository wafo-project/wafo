function h= nplot(ef,id,offset)
%NPLOT Normal probability plot of effects
%
% CALL:  nplot(ef,id)
%
%  ef = vector of average response, main effects and interaction effects
%       in the order given in id.
%  id = identification vector of main and interaction effects.
%       (default is standard order as given from the output of YATES)
%
% Two problems arise in the assessment of effects from unreplicated
% factorials: 1) occassionally real and meaningful high-order
% interactions occur, and 2) it is necessary to allow for selection.
% However, NPLOT by which effects are plotted on normal probability paper
% often provides an effective way of overcoming both difficulties.
%
% In the example below 3 effects with effects close to zero fit
% reasonably on a straight line. Those corresponding to A, B and AC do
% not. Conclusion: these effects are not easily explained as chance of
% occurences.
%
% Example:
%   y = [60 72 54 68 52 83 45 80]; % Responses to design D.
%   [ef,id] = yates(y);
%   nplot(ef,id) 
%
% See also  yates, fitmodel, getmodel, alias, cdr, sudg, ffd ,plotresponse, nplot

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


% Reference
% Daniel, C. (1959)
% Use of half-normal plot in interpreting factorial two-level
% experiments, Technimetrics, Vol 1, No 149. 
% Daniel, C. (1976)
% Applications of statistics to Industrial Experimentation.

%error(nargchk(1,3,nargin))
narginchk(1,3)
sz = size(ef);
n  = length(ef);
if prod(sz)==n,
  ef = ef(:);
else
  n = sz(1);
end


if nargin<2||isempty(id)
  % Make a list of all possible main effects and interaction effects
%------------------------------------------------------------------
  [y,id]= yates(ef);
end
% Secret option
if nargin<3||isempty(offset),offset = range(ef(2:end,:))/(3*n);end
[eff, ind] = sort(ef(2:end,:));

id = id(ind,:);
F = (0.5:n-1.5)'/(n-1);
nscore = invnorm(F);  % Normal score
ind = (diff(eff)==0);
ind = find(([0;ind] | [ind; 0])==0);
h1 = plot(eff(ind),nscore(ind),'.');
h2 = text(eff(ind)+offset,nscore(ind),id(ind,:),...
    'HorizontalAlignment','left','VerticalAlignment','middle');
eff(ind) = [];

if ~isempty(eff)
  nscore(ind) = [];
  id(ind,:) = [];
  iy = 1;
  ne = length(eff);
  hold on
  for ix=2:ne
    iy = iy+1;
    if ix==ne ||eff(ix) ~= eff(ix+1),
      tmp =  mean(nscore(ix-iy+1:ix));
      h1 = [h1 ; plot(eff(ix),tmp,'k.')];
      %h1 = text(eff(ix),tmp,num2str(iy));
      tmp2 = cellstr(id(ix-iy+1:ix,:)).';
      [tmp2{2,1:iy-1}] = deal(', ');
      tmp2{2,iy} = ' ';
      h2 = [h2;text(eff(ix)+offset/2,tmp,cat(2,[num2str(iy) ' '],tmp2{:}),...
    'HorizontalAlignment','left','VerticalAlignment','middle')];
      
      iy = 0;
    end
  end
  hold off
end
if nargout>0,
  h = [h1;h2];
end
xlabel('Effect')
ylabel('Nscore')



