function data=loaddata(filename)
% LOADDATA Loads a matrix from a text file.
%          
%  CALL: data=loaddata(filename)
%
%        data     = the output data matrix.
%        filename = a string containing the file to read.
%
%  This routine is used to allow m-functions to load data files
%  with their own name. ("fun.m" can't use "load fun.dat;v=fun;", but 
%  can use "v=loaddata('fun.dat');")
% Example:% 
%
% x=loaddata('sea.dat'); size(x)

% varname=filename;
% i=findstr(varname,'.');
% if ~isempty(i)
%   varname=varname(1:max(i)-1);
% end
% i=findstr(varname,'\'); % PC systems
% if ~isempty(i)
%   varname=varname(max(i)+1:length(varname));
% end
% i=findstr(varname,'/'); % Unix systems
% if ~isempty(i)
%   varname=varname(max(i)+1:length(varname));
% end

data = load(filename,'-ascii');
%eval(['data=' varname ';']);
%eval(['clear ' varname]);
