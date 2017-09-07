function [tr,ma,sa, form,not,dat] = trunmak(tro)
%TRUNMAK Split a transformation object into its pieces.
%
% CALL:  [tr,ma,sa,form,note,date] = trunmak(tro);
%
%  tro   = Transformation object
%  tr    = Two column table [x g(x)] or a pp-object giving the transformation. 
%  ma,sa = mean and standard deviation of the process.
%  form  = string identifying which form the transformation is
%          given as: 'pp' or 'table'.
%  note  = memorandum note
%  date  = date of creation
%
%
% See also  trmak


% TODO % This is not complete.

%error(nargchk(1,1,nargin))
narginchk(1,1)
switch class(tro),
 case 'struct',
  %OK
 otherwise
  error('This is not a transformation object!')
end

form = tro.form; 
tr  = tro.tr;   % Transformation
ma  = tro.mean; % mean
sa  = tro.std;  % standard deviation
not = tro.note;
dat = tro.date;



