function h=hos(A,kernel)
%HOS Oversmoothing Parameter.
%
% CALL:  h = hos(data,kernel)
%
%   h      = one dimensional maximum smoothing value for smoothing parameter
%            given the data and kernel.  size 1 x D
%   data   = data matrix, size N x D (D = # dimensions )
%   kernel = 'epanechnikov'  - Epanechnikov kernel.
%            'biweight'      - Bi-weight kernel.
%            'triweight'     - Tri-weight kernel.  
%            'triangluar'    - Triangular kernel.
%            'gaussian'      - Gaussian kernel
%            'rectangular'   - Rectanguler kernel. 
%            'laplace'       - Laplace kernel.
%            'logistic'      - Logistic kernel.
%  
%  Note that only the first 4 letters of the kernel name is needed.
% 
% The oversmoothing or maximal smoothing principle relies on the fact
% that there is a simple upper bound for the AMISE-optimal bandwidth for
% estimation of densities with a fixed value of a particular scale
% measure. While HOS will give too large bandwidth for optimal estimation 
% of a general density it provides an excellent starting point for
% subjective choice of bandwidth. A sensible strategy is to plot an
% estimate with bandwidth HOS and then sucessively look at plots based on 
% convenient fractions of HOS to see what features are present in the
% data for various amount of smoothing. The relation to HNS is given by:
% 
%           HOS = HNS/0.93
%
%  Example: 
%    % data = rndnorm(0, 1,20,1)
%  data = [-0.0233845632050972   0.9070186193622006;...
%           0.6529594866766634   1.3689145060433903;...
%           0.4477857310723146  -0.6311953712037597;...
%          -1.9256785038579962   0.5886257667993168;...
%          -0.5290011931824666  -0.3602090880229930];
%  assert(hos(data,'epan'), [1.73513679136905, 1.43948322577017]/0.93, 1e-10);
%  assert(hos(data,'biwe'), [2.05555487703312, 1.70530460760076]/0.93, 1e-10);
%  assert(hos(data,'triw'), [2.33418149081877, 1.93645545333964]/0.93, 1e-10);
%  assert(hos(data,'tria'), [1.90615281623682, 1.58135947458212]/0.93, 1e-10);
%  assert(hos(data,'gaus'), [0.783780547013426, 0.650230549961770]/0.93, 1e-10);
%  assert(hos(data,'rect'), [1.36382287194830, 1.13143825711994]/0.93, 1e-10);
%  assert(hos(data,'lapl'), [0.579817701798895, 0.481021358025369]/0.93, 1e-10);
%  assert(hos(data,'logi'), [0.438140924596874, 0.363485181466877]/0.93, 1e-10);
%  
% See also  hste, hbcv, hboot, hldpi, hlscv, hscv, hstt, kde, kdefun

% Reference:  
%  B. W. Silverman (1986) 
% 'Density estimation for statistics and data analysis'  
%  Chapman and Hall, pp 43-48 

%  Wand,M.P. and Jones, M.C. (1986) 
% 'Kernel smoothing'
%  Chapman and Hall, pp 60--63


%Tested on: matlab 5.3
% History:
% revised pab feb2005
% -updated example + see also line  
% revised pab 21.09.99
% 
% updated string comparisons
% from kdetools

h=hns(A,kernel)/0.93;