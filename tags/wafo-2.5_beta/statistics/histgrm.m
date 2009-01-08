function binwidth = histgrm(x,N,odd,scale,lintype)
%HISTGRM Plot histogram
% 
% CALL:  binwidth = histgrm(x,N,odd,scale)
%
%        binwidth = the width of each bin
%               x = the data
%               N = approximate number of bins wanted 
%                   (default depending on length(x))  
%             odd = placement of bins (0 or 1) (default 0)
%           scale = argument for scaling (default 0)
%                   scale = 1 yields the area 1 under the histogram
%         lintype : specify color and lintype, see PLOT for possibilities.
%
% Example:
%   R=rndgumb(2,2,1,100);
%   histgrm(R,20,0,1)
%   hold on
%   x=linspace(-3,16,200);
%   plot(x,pdfgumb(x,2,2),'r')
%   hold off
%   



% References: 
%  Holtsberg, Anders (1999)
%  Stixbox. A statistics toolbox for Matlab and Octave. 
%  Lund University
%  http://www.maths.lth.se/matstat/stixbox

% Tested on: Matlab 5.3
% History: 
% Revised by es 31.03.2000 Added default in help + lintype  
% Revised by jr 22.12.1999

if nargin < 2, N = []; end
if nargin < 3, odd = []; end
if nargin < 4, scale = []; end
if nargin < 5, lintype = []; end

if isempty(N),       N = ceil(4*sqrt(sqrt(length(x))));end
if isempty(odd),     odd = 0;end
if isempty(scale),   scale = 0;end
if isempty(lintype), lintype = 'b-';end


mn = min(x);
mx = max(x);
d = (mx - mn)/N*2;
e = floor(log(d)/log(10));
m = floor(d/10^e);
if m > 5
   m = 5;
elseif m > 2
   m = 2;
end
d = m * 10^e;
mn = (floor(mn/d)-1)*d - odd*d/2;
mx = (ceil(mx/d)+1)*d + odd*d/2;
limits = mn:d:mx;

f = zeros(1,length(limits)-1);
for i = 1:length(limits)-1
   f(i) = sum(x>=limits(i) & x<limits(i+1));
end

xx = [limits; limits; limits];
xx = xx(:);
xx = xx(2:length(xx)-1);
yy = [f*0; f; f];
yy = [yy(:); 0];
if scale, yy = yy/length(x)/d; end

H = ishold;
plot(xx,yy,lintype)
hold on
plot(limits,limits*0)
if ~H 
  hold off
end

if nargout > 0
   binwidth = d;
end


