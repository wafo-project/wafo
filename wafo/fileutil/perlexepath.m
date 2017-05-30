function p = perlexepath
%PERLEXEPATH Returns the path to the PERL executable
% 
%  CALL:  Str = perlexepath 
%  
%    Str  =  a string with the path to the directory containing 
%            the perl.exe function including a terminating slash 
%            or backslash, depending on the system on wich Matlab is run.
%            The executable files must be located in a directory found
%            by the command getenv('PATH').
%
%  The global variable PERLEXEPATH is also set.
% 
%  See also: perl

% History
% revised pab 8 Oct 2003 
%  -change PERLEXEPATH from global to persistent varaible  
% by pab 14.12.2000
% based on perl.m by Peter Acklam

%global PERLEXEPATH
persistent PERLEXEPATH
if ~isempty(PERLEXEPATH)
  p =  PERLEXEPATH;
  return
end
p = '';
if isunix, 
  PERLEXEPATH = p;
  return, 
end

if ispc
  path = getenv('PATH');
  k = find( path == pathsep );              % find path separators
  k = [ 0 k length(path)+1 ];               % find directory boundaries
  ndirs = length(k)-1;                      % number of directories
  for i = 1:ndirs
    dir = path( k(i)+1 : k(i+1)-1 );       % get i'th directory
    perlexe = fullfile(dir, 'perl.exe');   % path to perl.exe
    if exist(perlexe, 'file')
      PERLEXEPATH = fullfile(dir,'');
      p = PERLEXEPATH;
      return
    end
  end
  error('Perl binary not found. Is Perl installed properly?');
end
error('Unknown platform.');

