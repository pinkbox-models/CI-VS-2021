----------------------------------------------------------------------------------
 Theoretical relationship between two measures of spike synchrony: Correlation index and vector strength
 Matlab implementation 
----------------------------------------------------------------------------------

%%% Versions %%% 

+ Ver. 0.99 (October 17, 2021): Initial release of the code on GitHub. 

%%% Authors %%% 

Dominik Kessler (University of Oldenburg) dominik.kessler@uni-oldenburg.de
Go Ashida (University of Oldenburg) go.ashida@uni-oldenburg.de

%%% Aims %%% 

The functions and scripts of this repository can be used to ...
 1. simulate von Mises (vM) spike trains,
 2. calculate vector strength (VS), correlation index (CI) and shuffled autocorrelogram (SAC) from spike trains,
 3. calculate or estimate theoretical VS and CI via the concetration parameter kappa of the vM distribution,
 4. generate similar plots as in Figure 2, Figure 3A, and Figure 6 of Kessler et al.

%%% Contents %%% 

++ test scripts ++
testVM.m : script to produce similar plots as Figure 2 and Figure 3A of Kessler et al.
           by generating Poissonian spike trains with a von Mises distribution; 
           computing VS and CI; plotting rasters, phase histograms and SACs; and 
           comparing empirical VS and CI values to their theoretical relation.
testSACdecay.m  : script to reproduce Figure 6 of Kessler et al. about the effect of the data length on the SAC.

++ internally called functions ++
genPhaseLock.m  : function to generate Poissonian spike trains with a von Mises spiking distribution.
calcPhaseHist.m : function to calculate the phase histogram and vector strength (VS) of spike trains.
calcSAC.m    : function to calculate the shuffled autocorrelogram (SAC) and correlation index (CI) for spike trains.
estimateCI.m : function to estimate CI from VS with their theoretical relation for the vM distribution
estimateVS.m : function to estimate VS from CI with their theoretical relation for the vM distribution


%%% Reference %%% 

Kessler D, Carr CE, Kretzberg J, Ashida G (2021)
"Theoretical relationship between two measures of spike synchrony: Correlation index and vector strength"
(Submitted)

%%% Copyright/License %%% 

CC0 1.0 Universal (Public Domain Dedication).
See LICENSE.txt for the statement.   

