function tro = trmak(tr,ma,sa,note)
%TRMAK Put together a transformation object.
%
% CALL:  tro = trmak(tr,ma,sa,note);
%
%  tro   = Transformation struct with the fields:
%          tro.form = string identifying the transformation object,
%                     'pp', or 'table'.
%          tro.tr   = tr (see tr below)
%          tro.mean = ma
%          tro.std  = sa
%          tro.note = memorandum note
%          tro.date = createion date and time.
%  tr    = Two column table [x g(x)] or a pp-object giving the transformation. 
%  ma,sa = mean and standard deviation of the process.
%
% Examples:
%  tro = trmak; % Make an empty transformation object. 
%  sa = 2;  skew = .20; kurt = (4*skew/3).^2+3; ma =0;
%  tro1 = trmak(hermitetr([],[sa skew,kurt,ma]),ma,sa);
% 
% See also  trunmak

% Tested on: Matlab 5.3
% History:
% By pab 02.04.2001


% TODO % This is not complete. Implement this structure in all XXX2tr programs.

% Legal field names:
nam1 = { 'form','tr','mean','std','note','date'}; 


if nargin==0,  % make an empty object
  tr =[];ma =[]; sa =[];note =[];
else
  %error(nargchk(3,4,nargin))
  narginchk(3,4)
  if nargin<4, note = []; end
  switch class(tr),
   case 'struct',
    nam = fieldnames(tr);
     
    nam2 = unique(cat(1,nam(:),nam1(:)));
    if length(nam2)==length(nam), % tr is already a transformation struct.
      tro = tr;
      return
    elseif isfield(tr,'form') % We are dealing with a PP-object
      tro.form  = tr.form;
    else
      error('TR input is not a valid transformation function')
    end
   case 'double', 
    mvrs = version;
    ix   = find(mvrs=='.');
    mvrs = str2num(mvrs(1:ix(2)-1));
    trsz = size(tr);
    if (tr(1)==10) && prod(trsz)==length(tr) && mvrs<=5.2,
      tro.form ='pp';  % Old pp-form object. Matlab version 5.2 and below
    elseif trsz(2)==2,
      tro.form  = 'table';
    else
      error('TR input is not a valid transform function')
    end
  end
end
tro.tr    = tr;    % Transformation
tro.mean  = ma; % mean
tro.std   = sa; % standard deviation
tro.note  = note;
tro.date  = datestr(now);



