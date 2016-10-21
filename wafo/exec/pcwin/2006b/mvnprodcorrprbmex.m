
%MVNPRODCORRPRBMEX Computes multivariate normal probability 
%                with product correlation structure.
%
%  CALL [value,error,inform]=mvnprodcorrprbmex(rho,A,B,abseps,releps,useBreakPoints);
%
%     RHO    REAL, array of coefficients defining the correlation
%            coefficient by:
%                correlation(I,J) =  RHO(I)%RHO(J) for J/=I
%            where 
%                1 <= RHO(I) <= 1
%     A		 REAL, array of lower integration limits.
%     B		 REAL, array of upper integration limits.
%	       NOTE: any values greater the 10, are considered as
%                   infinite values.
%     ABSEPS REAL absolute error tolerance.
%     RELEPS REAL relative error tolerance.
%     USEBREAKPOINTS = 1 If extra integration points should be used
%                        around possible singularities
%                      0 If no extra
%  
%     ERROR  REAL estimated absolute error, with 99% confidence level.
%     VALUE  REAL estimated value for the integral
%     INFORM INTEGER, termination status parameter:
%            if INFORM = 0, normal completion with ERROR < EPS;
%            if INFORM = 1, completion with ERROR > EPS and MAXPTS 
%                           function vaules used; increase MAXPTS to 
%                           decrease ERROR;
%
%
% This file was successfully compiled for matlab 5.3
% using Compaq Visual Fortran 6.1, and Windows 2000.
% The example here uses Fortran77 source.
% First, you will need to modify your mexopts.bat file.
% To find it, issue the command prefdir(1) from the Matlab command line,
% the directory it answers with will contain your mexopts.bat file.
% Open it for editing. The first section will look like:
%
%rem %%%%%%%%%***********************************************************
%rem General parameters
%rem ********************************************************************
%set MATLAB=%MATLAB%
%set DF_ROOT=C:\Program Files\Microsoft Visual Studio
%set VCDir=%DF_ROOT%\VC98
%set MSDevDir=%DF_ROOT%\Common\msdev98
%set DFDir=%DF_ROOT%\DF98
%set PATH=%MSDevDir%\bin;%DFDir%\BIN;%VCDir%\BIN;%PATH%
%set INCLUDE=%DFDir%\INCLUDE;%DFDir%\IMSL\INCLUDE;%INCLUDE%
%set LIB=%DFDir%\LIB;%VCDir%\LIB
%
% then you are ready to compile this file at the matlab prompt using the
% following command:
%  mex -O mvnprodcorrprbmex.f
