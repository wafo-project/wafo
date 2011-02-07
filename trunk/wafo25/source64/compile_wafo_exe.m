%% COMPILE_WAFO_EXE Script that compiles WAFO exe-files.
%
% COMPILE_WAFO_EXE compiles WAFO Fortran files and copies them 
% to the correct location. In order to use this you must have 
% a working fortran compiler installed in the path and wafo installed.
%
% Example
% cd(fullfile(waforoot,'source'))
% compile_wafo_exe
% 




try
  %SOURCE_DIR = fullfile(waforoot,'source');
  OUT_DIR    = wafoexepath;
catch
  error('WAFO must be in the PATH! Run initwafo(''full'') to set correct PATH.')
end

%p0 = pwd;
%cd(SOURCE_DIR)



if ~exist(OUT_DIR,'dir')
  [success,msg,msgid] = mkdir(outdir);
  if ~success
    error(msg)
  end
end


%% Choose fortran compiler
disp('A fortran compiler must be installed on the system prior to running this script. ')
disp(' ')
disp('Free fortran compiler may be downloaded from http://www.g95.org/.')
disp(' ')
disp('Choose a Fortran-compiler! (default ''g95'') ')
disp(' ''g95''   : G95 Fortran Compiler')
disp(' ''ifort'' : Intel Fortran Compiler')
disp(' ''df''    : Compaq visial fortran 6.1')
disp(' ')

FC = input('Type in your selection here (in quotes and press enter): ');
if isempty(FC)
  FC = 'g95';
  % FC = 'gfortran'
  % FC = 'ifort'; % Intel Fortran Compiler
  % FC = 'df';    % Compaq visial fortran 6.1
end

% Set your compiler parameters  
OUTPAR = '-o ';
if strcmp(FC,'g95')
  PARAM = '-O3 -ffixed-form'; % -ftrace=full'; 
elseif strcmp(FC,'df') %  compaq visial fortran 6.1
  PARAM = '/fast /fixed /transform_loops'; % /traceback /fast  
  OUTPAR = '/exe:';
else
  PARAM = '-O -fixed';
  % PARAM  = '-automatic' ; % DIGTAL (alpha)
  % PARAM = '-O -fixed -Bdynamic'; % SOLARIS (sol2)
end

%% FORTRAN FILE COMPILATION
%The files written in Fortran is installed by the following commands in the matlab prompt:

try
 miss = cell(1,0);
 Nmiss = 0;
  f77lib  = cell(1,9);
  f77file = cell(1,9);
  [f77file{1:2}] = deal('down2cc.f','simduff.f');
  [f77lib{1:2}]  = deal({''});
  
  f77file{3} = 'cov2mmpdfreg.f'; % old minmax
  f77lib{3}  = {'dsvdc.f','mregmodule.f'};
  
  [f77file{4:9}] = deal('cov2acdf.f','cov2mmpdf.f','cov2tccpdf.f',...
    'cov2thpdf.f','cov2mmtpdf.f','rindd.f');
  [f77lib{4:9}]  = deal({'intmodule.f', 'jacobmod.f', 'rind70.f'});
 
  FC
  
 
  ext = '.exe';
  foutfile = strrep(f77file,'.f',ext);
  h = fwaitbar(0,[],'this may take a while');
  Nf = numel(f77file);
  for ix=1:Nf
    outfile = sprintf('%s%s',OUTPAR,fullfile(OUT_DIR,foutfile{ix}));
    infile = sprintf('%s ', f77lib{ix}{:},f77file{ix});
    fwaitbar(ix/Nf,h,sprintf('Compiling %s',foutfile{ix}))
    errorCode = dos(sprintf('%s %s %s %s',FC,PARAM,outfile,infile));
    if errorCode
      Nmiss = Nmiss+1;
      miss{Nmiss} = outfile;
      warning('WAFO:COMPILE_WAFO_EXE','Unable to compile %s',outfile)
    end
  end
  close(h)

 
  
  %clean-mod:
	delete('*.mod')
  if isempty(miss)
    disp('Successfully compiled all .exe-files!')
  else
    disp(sprintf( 'Failed to compile: %s\n',miss{:}))
  end
catch
  disp('Unexpected termination of compilation of fortran-files!')
end


%cd(p0) % move to the starting directory
