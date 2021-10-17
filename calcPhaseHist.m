function [PH,PHtv,VS] = calcPhaseHist(SPin,T1,T2,NB,FQ)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculating phase histogram and vector strength shuffled of spike trains 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Inputs
%   SPin: 1-D cell array of spike time vectors [ms] 
%   T1: starting time [ms] for the analysis 
%   T2: ending time [ms] for the analysis 
%   NB: number of bins used for the phase histogram
%   FQ: frequency [Hz] 
%
%  +++ Notes +++
%   Spike times between T1 and T2 will be used for the analysis. 
%
% Outputs
%   PH: calculated period histogram values
%   PHtv: time vector for period histogram [cycle]
%   VS: vector strength
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% preprocessing to remove spike times out of the analysis time window
Nreps = length(SPin);
SP = cell(1,Nreps);
for c = 1:Nreps
  v = SPin{c}; 
  SP{c} = v( v>=T1 & v<=T2 );  
end

% vector strength
VSnn = 0; 
VScc = 0; 
VSss = 0; 
for c = 1:Nreps
  VSnn = VSnn + length(SP{c});
  VScc = VScc + sum( cos( 2 * pi * FQ * SP{c} / 1000 ) );
  VSss = VSss + sum( sin( 2 * pi * FQ * SP{c} / 1000 ) );
end
VS = sqrt( VScc*VScc + VSss*VSss ) / VSnn; 

% period histogram
PHtv = (0:NB-1)/NB;
PH = zeros(1,NB); 
for c = 1:Nreps
 PHidx = floor( mod((FQ*SP{c}/1000),1) * NB) + 1;  % phase index vector
 for i=1:length(PHidx)
   PH(PHidx(i)) = PH(PHidx(i)) + 1;
 end
end

end % [eof]
