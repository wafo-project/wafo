% Paper REC in WAFO Toolbox. 
% Version 2.1.1   15-Sep-2005 
% 
%      Brodtkorb, P.A., Myrhaug, D. and Rue, H (1999). Joint distribution of
%      wave height and wave crest velocity from reconstructed data.
%      Proceedings of 9th ISOPE Conference, Brest, Vol III, pp. 66-73 
%
% 
% recdemo       - Show Rec figures with point and click interface 
% recfig        - Callback implementing functions of RecDemo 
% recfig1       - Location of Gullfaks C and Statfjord A platforms in The North Sea 
% recfig10      - Probability of exceeding H: Model (dash); data (dots) 
% recfig11      - Probability of exceeding V: Model (dash); data (dots) 
% recfig12      - The conditional probability of exceeding V given H: Model (dash); data (dots) 
% recfig2       - 10 minutes mean values of wind  (dash) and direction (solid) 
% recfig3       - Example of the reconstructed (dots) and original data set (pluses) 
% recfig4       - Estimated spectral density (solid) and 95% confidence intervals (dots) 
% recfig5       - Transfer function, g, versus the crossing level u 
% recfig6       - Variability of simulated e(g(u)-u) (circles) 
% recfig7       - Estimated Weibull parameters versus h=H/Hrms (circles) 
% recfig8       - Conditional mean and standard deviation of V given H: 
% recfig9       - Joint distribution of  V and H: 
% recinit       - Setup all global variables 
% recintro      - Info about RECDEMO: A statistical procedure for reconstruction of 1D signals. 
% recleanup     - Clears global variables defined and used by the RECDEMO 
