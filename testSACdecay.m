%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Script to test the effect of limited data length on SAC 
% (cf. Figure 6 of Kessler et al.). 
% Four different data lengths (20, 30, 40, and 50 ms) are tested.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% 1. Set parameters for spike train generation 
DT = 0.002;        % [ms] time step
T = 100;           % [ms] total time length
N = round(T/DT);   % number of time steps
tv = (0:N-1) * DT; % [ms] time vector

M = 2000; % number of trials
F = 500;  % [Hz] stimulus frequency
R = 0.8;  % target vector strength value
L = 200;  % [spikes/sec] mean rate 
P = pi;   % initial phase

%% 2. Generating inhomogeneous poisson spikes
disp('making spikes');
[SP,K] = genPhaseLock(M,N,F,R,L,P,DT); 
SPin = cell(1,M);
for c = 1:M
  SPin{c} = tv( logical(SP(c,:)) );
end

%% 3. Choose analysis parameters
NB = 41;   % number of bins (for phase histogram)
TL = 10;   % [ms] maximum time delay for SAC
BW = 0.05; % [ms] SAC bin width

% data length = 50ms
T1_D50 = 30; % [ms] analysis start time
T2_D50 = 80; % [ms] analysis end time

% data length = 40ms
T1_D40 = 30; % [ms] analysis start time
T2_D40 = 70; % [ms] analysis end time

% data length = 30ms
T1_D30 = 30; % [ms] analysis start time
T2_D30 = 60; % [ms] analysis end time

% data length = 20ms
T1_D20 = 30; % [ms] analysis start time
T2_D20 = 50; % [ms] analysis end time

%% 4. Calculating VS and SAC for different data lengths
disp('calculating PH');
[PH_D50, PHtv_D50, VS_D50] = calcPhaseHist(SPin, T1_D50, T2_D50, NB, F);
[PH_D40, PHtv_D40, VS_D40] = calcPhaseHist(SPin, T1_D40, T2_D40, NB, F);
[PH_D30, PHtv_D30, VS_D30] = calcPhaseHist(SPin, T1_D30, T2_D30, NB, F);
[PH_D20, PHtv_D20, VS_D20] = calcPhaseHist(SPin, T1_D20, T2_D20, NB, F);

disp('calculating SAC');
[SAC_D50, SACtv_D50, CI_D50] = calcSAC(SPin, BW, T1_D50, T2_D50, TL);
[SAC_D40, SACtv_D40, CI_D40] = calcSAC(SPin, BW, T1_D40, T2_D40, TL);
[SAC_D30, SACtv_D30, CI_D30] = calcSAC(SPin, BW, T1_D30, T2_D30, TL);
[SAC_D20, SACtv_D20, CI_D20] = calcSAC(SPin, BW, T1_D20, T2_D20, TL);

%% 5. Calculate theoretical values
VSthr = besseli(1,K) ./ besseli(0,K);
CIthr = besseli(0,2*K) ./ besseli(0,K).^2; 
SACest = besseli(0,2*K*cos(pi*F*SACtv_D50/1000)) ./ besseli(0,K).^2; 

%% 6. Compute decay factors for different data lengths
% 50ms
SACdecay_D50 = ones(1,length(SACtv_D50));
SACdecay_D50(SACtv_D50<0) = 1 + SACtv_D50(SACtv_D50<0)/(T2_D50-T1_D50);
SACdecay_D50(SACtv_D50>0) = 1 - SACtv_D50(SACtv_D50>0)/(T2_D50-T1_D50);

% 40ms
SACdecay_D40 = ones(1,length(SACtv_D40));
SACdecay_D40(SACtv_D40<0) = 1 + SACtv_D40(SACtv_D40<0)/(T2_D40-T1_D40);
SACdecay_D40(SACtv_D40>0) = 1 - SACtv_D40(SACtv_D40>0)/(T2_D40-T1_D40);

% 30ms
SACdecay_D30 = ones(1,length(SACtv_D30));
SACdecay_D30(SACtv_D30<0) = 1 + SACtv_D30(SACtv_D30<0)/(T2_D30-T1_D30);
SACdecay_D30(SACtv_D30>0) = 1 - SACtv_D30(SACtv_D30>0)/(T2_D30-T1_D30);

% 20ms
SACdecay_D20 = ones(1,length(SACtv_D20));
SACdecay_D20(SACtv_D20<0) = 1 + SACtv_D20(SACtv_D20<0)/(T2_D20-T1_D20);
SACdecay_D20(SACtv_D20>0) = 1 - SACtv_D20(SACtv_D20>0)/(T2_D20-T1_D20);

%% 7. Plotting (cf. Figure 6 in Kessler et al.)
figure(3);
set(gcf, 'PaperPositionMode', 'auto');
set(gcf, 'Position', [450 60 1120 420]);

% subplot1: SAC curve (50ms and 40ms)
subplot(1,2,1); cla; hold on;
plot(SACtv_D50, SACest, 'k-'); % SAC for infinite data length
plot(SACtv_D50, SAC_D50,'b-'); % calculated SAC (50ms)
plot(SACtv_D50, max(SACest).*SACdecay_D50, 'b:'); % theo decay line (50ms)
plot(SACtv_D40, SAC_D40,'g-'); % calculated SAC (40ms)
plot(SACtv_D40, max(SACest).*SACdecay_D40, 'g:'); % theo decay line (40ms) 

text(-8, max(SAC_D50)*1.24, sprintf('CI_t_h_e_o = %.4f',CIthr));
text(-8, max(SAC_D50)*1.16, sprintf('CI_5_0_m_s = %.4f',CI_D50));
text(-8, max(SAC_D40)*1.08, sprintf('CI_4_0_m_s = %.4f',CI_D40));
box off;
set(gca,'TickDir','out');
yticks([0 1 2 3]);
yticklabels({'0','1','2','3'});
ylim([0 max(SAC_D50)*1.24]);
ylabel('normalized #coincidences');
xlabel('delay (ms)');

% subplot 2: SAC curve (30ms and 20ms)
subplot(1,2,2); cla; hold on;
plot(SACtv_D30, SACest, 'k-'); % SAC for infinite data length
plot(SACtv_D30, SAC_D30, 'r-'); % calculated SAC (30ms)
plot(SACtv_D30, max(SACest).*SACdecay_D30, 'r:'); % theo decay line (30ms)
plot(SACtv_D20, SAC_D20,'y-'); % calculated SAC (20ms)
plot(SACtv_D20, max(SACest).*SACdecay_D20,'y:'); % theo decay line (20ms) 

text(-8, max(SAC_D30)*1.24, sprintf('CI_t_h_e_o = %.4f',CIthr));
text(-8, max(SAC_D30)*1.16, sprintf('CI_3_0_m_s = %.4f',CI_D30));
text(-8, max(SAC_D20)*1.08, sprintf('CI_2_0_m_s = %.4f',CI_D20));
set(gca,'TickDir','out');
box off;
yticks([0 1 2 3]);
yticklabels({'0','1','2','3'});
ylim([0 max(SAC_D50)*1.24]);
xlabel('delay (ms)');
