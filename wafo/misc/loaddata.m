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
%
% Example:
%
%   x = loaddata('sea.dat'); 
%   assert(size(x), [9524, 2])
%   assert(x(1:3,2)', [ -1.2004945, -1.0904945, -0.79049454], 1e-6)
%   assert(x(1:3,1)', [0.05, 0.3, 0.55], 1e-6)

data = load(filename,'-ascii');
