function [SAC,SACtv,CI,CN,Nsp] = calcSAC(SPin,BW,T1,T2,TL)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculating shuffled autocorrelogram of spike trains 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Inputs
%   SPin: 1-D cell array of spike time vectors [ms] 
%   BW: SAC time bin width [ms] 
%   T1: starting time [ms] for the analysis 
%   T2: ending time [ms] for the analysis 
%   TL: maximum time difference [ms] for SAC calculation
%
%  +++ Notes +++
%   The total number of spikes should be below ~20000. 
%   Spike times between T1 and T2 will be used for constructing the SAC. 
%   The values of resulting time vector will be between -TL and +TL.  
%
% Outputs
%   SAC: calculated SAC values
%   SACtv: time vector [ms] for SAC 
%   CI: correlation index (SAC at time zero) 
%   CN: normalization factor = (#reps)*(#reps-1)*(rate^2)*(bin width)*(data length) 
%   Nsp: number of spikes 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Release: 17 October 2021
%  Dominik Kessler (dominik.kessler@uni-oldenburg.de)
%  Go Ashida (go.ashida@uni-oldenburg.de)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% preprocessing to remove spike times out of the analysis time window
Nreps = length(SPin);
SP = cell(1,Nreps);
for c = 1:Nreps
  v = SPin{c}; 
  SP{c} = v( v>=T1 & v<=T2 );  
end

% SAC time/data vectors
Nbins = ceil(TL/BW); 
SACtv = (-Nbins:Nbins) * BW;  % total number of bins = 2*Nbins + 1
SAC = zeros(1, length(SACtv)); 
Nsp = 0;  % number of spikes

% main loop 
for c = 1:length(SP)
  Nsp = Nsp + length(SP{c});
  for k = [ 1:c-1, c+1:length(SP) ]  % loop across trials
    if( ~isempty(SP{c}) && ~isempty(SP{k}) )
      d = SP{c} - SP{k}'; % matrix containing spike time differences
      v = d( d>-TL-BW/2 & d<TL+BW/2 );  % data within SAC time range
      
      % histogram counts
      for j = 1:length(v)
        i = round(v(j)/BW) + Nbins + 1;
        SAC(i) = SAC(i) + 1;
      end
    end
  end
end

% normalization
T = T2-T1;  % [ms] data length
R = Nsp / (Nreps*T);  % [spikes/ms] average spike rate
CN = Nreps * (Nreps-1) * R^2 * BW * T;  % normalization factor (Louage norm)

% final results 
SAC = SAC/CN; 
CI = SAC(Nbins+1); 

end % end of function
