function [Triv,Math,Contr,Exim] = fnames
% FNAMES Creates string matrices with standard function names
%
% CALL: [triv,math,contr,exim] = fnames  
%
%  triv  = string matrices of some general language commands 
%  math  = elementary math operators,
%  contr = control flow commands 
%  exim  = "export-import" commands.
%
% See also: findref

% History
% revised pab 19.10.2000
%  - updated header info
%  Kirill K. Pankratov, kirill@plume.mit.edu
%  02/06/95
  
%w    = what('elmat');
%w    = strrep(w.m,'.m',''); % remove .m endings
%Triv = strvcat(w{:});
%
%w    = what('elfun');
%w    = strrep(w.m,'.m',''); % remove .m endings
%Math = strvcat(w{:}); 
%  
%w     = what('lang');
%w     = strrep(w.m,'.m',''); % remove .m endings
%     or
%w     = iskeyword;  
%Contr = strvcat(w{:}); 
  
%w     = what('iofun');
%w     = strrep(w.m,'.m',''); % remove .m endings
%Exim  = strvcat(w{:}); 
  
% "Trivia" = Elementary matrices and matrix manipulation.    
Triv = ['ones      '; 'zeros     '; 'rand      '; 'randn     ';
        'find      '; 'max       '; 'min       '; 'all       ';
        'any       '; 'format    '; 'sum       '; 'cumsum    ';
        'prod      '; 'cumprod   '; 'length    '; 'size      ';
        'abs       '; 'log       '; 'reshape   '; 'sort      ';
        'i         '; 'j         '; 'inf       '; 'nan       ';
        'eye       '; 'diag      '; 'fliplr    '; 'flipud    ';
        'function  '; 'nargchk   '; 'nargoutchk'; 'varargin  ';
	'varargout '; 'nargin    '; 'nargout   '; 'script    ';
	'repmat    '; 'eps       '; 'ans       '; 'isnan     ';
	'isinf     '; 'isfinite  '; 'numel     '; 'ndims     ';
	'disp      '; 'isempty   '; 'islogical '; 'logical   ';
        'realmax   '; 'realmin   ';];
 % "Math" - some elementary math functions
Math = ['abs     '; 'acos    '; 'acosh   '; 'acot    '; 'acoth   '; 
        'acsc    '; 'acsch   '; 'asec    '; 'asech   '; 'asin    ';
        'asinh   '; 'atan    '; 'atan2   '; 'atanh   '; 'ceil    ';
        'conj    '; 'cos     '; 'cosh    '; 'cot     '; 'coth    ';
	'mod     '; 'csc     '; 'csch    '; 'exp     '; 'fix     '; 
	'floor   '; 'log2    '; 'imag    '; 'log     '; 'log10   '; 
	'real    '; 'rem     '; 'pow2    ';
        'round   '; 'sec     '; 'sech    '; 'sign    '; 'sin     ';
        'sinh    '; 'sqrt    '; 'tan     '; 'tanh    '; 'pi      ';
        'angle   '; 'wrap    '; 'unwrap  '; 'isreal  '; 
	'complex '; 'cplxpair';];

 % "Controls" = matlab keywords
Contr = [ 'break     '; 'case      ';  'catch     ';
	  'continue  '; 'else      ';  'elseif    ';
	  'end       '; 'for       ';  'function  ';
	  'global    '; 'if        ';  'otherwise ';
	  'persistent'; 'return    ';  'switch    ';
	  'try       '; 'while     ';];   

 % "Export-import"
Exim = ['load   '; 'fprintf'; 'fread  '; 'fscanf ';
        'save   '; 'fwrite '; 'eval   '; 'feval  '];

