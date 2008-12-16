function [filename] = writecov(cov,t,nr)
% WRITECOV Write the covariance R(t) and its derivatives to disk
% 
%  CALL filename = writecov(cov,t,nr)
% 
% filename = cellarray of filenames made.
%      cov = [R,R',R'',R''',R''''] matrix of covariance R(t) and
%            its derivatives as returned from spec2cov2, size N X M
%      t   = vector of times, length N. (default [])
%      nr  = number or derivatives to write. (default M-1)
%
% WRITECOV write the covariance R(t) and its derivatives to to ascii files  
% named Cd*.in on the hard disk. Each Cd0.in, Cd1.in,... contains either 
% 1 column with on of R(t), R'(t),R''(t),... only or 2 columns with
% t in 1'st column and one of R(t), R'(t),R''(t),... in the second column. 
% These files are used by the fortran programs: 
%   cov2mmpdfreg, cov2mmtpdf, cov2tccpdf, cov2tpdf and cov2thpdf.
%
% See also spec2mmpdfreg, spec2mmtpdf, spec2tccpdf, spec2tpdf and spec2thpdf.

% History: tested on matlab 7
% Revised pab august 2007
% -Added t as seperate input.
% -each Cd0.in, Cd1.in,... contains either 1 column with Rx^i(t) or 2 columns with
% t and Rx^i(t)
% Revised pab may 2007
% -each Cd0.in, Cd1.in,... now only contains 2 columns with t, and Rx^i(t)

n = size(cov,1);
m = size(cov,2);
formatStr = '%12.10f %12.10E \n';

if nargin<3 || isempty(nr)
  nr = m-1;
end
if nargin<2 || isempty(t)
  t = [];
  formatStr = '%12.10E \n';
elseif length(t)~=n
  error('Length of t must equal size(cov,1)!')
end

if m<nr+1
  error(['You must supply at least nr=', int2str(nr), ' derivatives.'])
end

filename = cell(1,nr+1);

for k=0:nr
  filename{k+1} =sprintf('Cd%d.in',k);
end
cleanup(filename{:})


for k=0:nr
  %   covar=[cov(:,1), cov(:,k+2), zeros(n,3)];
  %   fid=fopen(filename{k+1},'wt');
  %   fprintf(fid,'%12.10f %12.10E %4.2f %4.2f %4.2f\n',covar');
  covar=[t, cov(:,k+1)];
  fid=fopen(filename{k+1},'wt');
  fprintf(fid,formatStr,covar.');
  fclose(fid);
end
