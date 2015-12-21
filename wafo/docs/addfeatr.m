% ADDFEATR  How To Add New Features to WAFO
%
% NEW Functions 
%~~~~~~~~~~~~~~~
% General conventions:
% Each line of code in the functions must not exceed 76 columns
% to ease the readability of the code. Function names must be written in lower 
% case letters without any underscores. Also local variable names should 
% preferably be in lower case letters. Global variables should be written in 
% upper case letters in order to avoid confusion with local variables. 
%
% Help header:
% Generally a function included in WAFO should be well commented and include 
% the following information in the given order in the help header:
%
% Header: (H1 line) Name of the function (in uppercase letters) followed by a
%         concise and definitive description appropriate for retrieval purposes by
%          the MATLAB \verb+lookfor+ command. Not to exceed one line.  
% Call  : Specifying the calling prototype, i.e. output argument(s),
%         function name and input argument(s). Function name in lower case letters. 
% Outputs/inputs: A list describing the output and input arguments in
%         the order they appear in the calling prototype. The description should
%         also indicate the data type of the argument. An optional input argument is
%         indicated by specifying the default value enclosed in parenthesis. Each
%         line should be aligned so that all the equal-signs are in the same column.
% Description: Detailed explanation of what the function does, the
%         assumptions made, the limitations and how the algorithm works.  This
%         should be as detailed as possible without exceeding one page.
% Side effects: If any like changing/using global variables, 
%         changing properties of figure windows etc.
%         should also be notified in the help header
% Example(s): An example of how to use the function in practice is desirable.
% See also  Comma separated list of related functions or functions which this 
%           routine calls.
% Secret help:
% Immediately after the help header the following information should be given: 
%
% References: A complete reference from which the user can obtain further information
% Tested on : Specifies on which matlab version the function has been tested.
% History   : Revision log of the function in chronological order. 
%
% See jonswap for an example on how it looks like in the WAFO toolbox
%
%
%  NEW DATASETS:
% ~~~~~~~~~~~~~~
%
%  Signals: 
%
%  [1] Choose a filename  (preferably with less than eight characters)  e.g. "gullfaks" 
%
%  [2] Create an ascii data file named  "gullfaks.dat", which contains your
%      signal in vertical format with sampled values in columns.
%      The documentation m-file explains what each column is a measure of.
%      Sampling times if given should be in the first column.
%
%  [3] Create an documentation file named "gullfaks.m" with the 
%      following structure: (all in lower case letters if not explicitly stated otherwise)
%
% Title         : name of the m file in upper case letters followed by a 
%                 descriptrive text. Not to exceed one line.
% Call          : indicating how to access the data
% Size          : The size of data.
% Sampling rate :
% Device        : Measuring device used to sample the data.
% Source        : Indication of the original source of the data
% Format        : Specifying the format of the data
% Description   : and other relevant facts of the data , \ie{} where and when the 
%                 data were measured. An indication of the quality of the data should 
%                 also be included.  Descriptive measures of the data like significant 
%                 wave height, Hm0 , and peak period, Tp, may also be given here.  
%                 Also restrictions on the use of the data if any should be given here.
% See also      : a comma separated list of related files or functions
%
%  Example:
%
% GULLFAKS  surface elevation measured at Gullfaks C 30.12.1990
%              
%  CALL:  xn = load('gullfaks.dat')
%  
% Size             :    2560 X 2
% Sampling Rate    :    2.5 Hz
% Device           :    EMI laser
% Source           :    STATOIL 
% Format           :    ascii, c1: time c2: surface elevation
% Description      :
%   The wave data was measured 24th December 1989 at the Gullfaks C platform  
%   in the North Sea from 17.00 to  21.20. The  water depth of  218 m is regarded 
%   as deep water  for the most  important wave components.
%   There are two EMI laser sensors named 219 and 220. This data  set  is obtained from
%   sensor 219, which is  located in the
%   Northwest corner  approximately two platform leg diameters away from   
%   the closest  leg.  
%   Thus the wave elevation is not expected  to be significantly 
%   affected by  diffraction effects for incoming waves in the western sector.   
%   The wind direction for this period is from the southwest.
%   Some  difficulties in calibration of the instruments have been reported
%   resulting in several consecutive measured values being equal or nearly equal 
%   in the observed data set.
% 
%   Hm0 = 6.8m, Tm02 = 8s, Tp = 10.5
% 
% See also  gullfaks.jpg
%
%
% NEW  DEMONSTRATIONS / PUBLISHED PAPERS:
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%  To add new demonstration scripts to the DEMOS directory:
%
%  [1] Choose a name for your demo and a short prefix for files
%      implementing your demo. For example, "isope99" and "is".
%
%  [2] Create the new subdirectory accordingly after the chosen name:
%
%         WAFO/DEMOS/isope99. 
%
%  [3] Create the m-files which implement your demo:
%      isdemo    - starts a user controlled demonstration with choices
%      isinit    - sets up data structures, globals etc...
%      isfig1,2..- called from isdemo. 
%      isintro   - help file describing the purpose of the demo
%      iscleanup - clears all globals created by the demo and restores the 
%                  original workspace
%
%  [4]  Specialized tools/functions not available in the official 
%       release of WAFO needed for generating the figures should be put into 
%       WAFO/DEMOS/isope99/private subdirectory.
%
%      These are the rules for the scripts:
%
%      I.   One m-file creates one complete figure, not a series
%              of figures or a subplot of a figure.
%      II.  If some variables in one script is needed by another script
%              they  must be saved as global variables.
%      III. No pause or prints statements allowed.
%      IV.  As far as possible try to use the tools in WAFO
%
%      Inspection of existing demonstrations will help in following these rules.
%
% 
% 
%  To add new demonstration scripts to the PAPER directory is 
%  similar to the DEMOS except that the scripts must recreate figures 
%  in  published articles or technical reports.
%
% 


% history
% by pab 14.12.1999
% 

more on
help addfeatr
more off
