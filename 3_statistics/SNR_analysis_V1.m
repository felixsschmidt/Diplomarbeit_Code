% SNR evaluation
%% Set up
clear, clc, close all
format compact
cd 'C:\Users\Felix\Desktop\Close Loop\CMC\mat'
currentDir = pwd;
addpath(genpath(fullfile(currentDir, 'data')));
addpath(genpath(fullfile(currentDir, 'code')));
addpath(genpath(fullfile(currentDir, 'kakearney-boundedline-pkg-50f7e4b')));
sb_figDir = 'C:\Users\Felix\Desktop\Close Loop\CMC\mat\data\figures';
sb_resultDir = 'C:\Users\Felix\Desktop\Close Loop\CMC\mat\data\results';

%% SNR Evaluation
load 'SNR estimates.mat'
%load 'CMC_values_in_both_conditions.mat'
% descriptive results
M_SNR = [mean(SNR_high) mean(SNR_low)];
SD_SNR = [std(SNR_high) std(SNR_low)];
[SNR_high SNR_low];

% inferntial results
[H,P,CI,STATS] = ttest(rau(SNR_high), rau(SNR_low), 'tail', 'both');

% plot CMC-SNR correlation
% plot results
figure
plot(cmcValues(:,1), SNR_high, 'o', 'MarkerSize', 4)
title('Correlation of CMC and SNR in high MEP trials')
xlabel('SNR'), ylabel('Coherence'), hold on

plot(cmcValues(:,2), SNR_low, 'o', 'MarkerSize',4 )
title('Correlation of CMC and SNR in high and low MEP trials')
xlabel('SNR'), ylabel('Coherence')
legend('highMEP', 'lowMEP')

[Cmc_SNR_r', Cmc_SNR_p']
Phs_SNR_r
[v i] = min(Phs_SNR_r(:,5)), Phs_SNR_r(i,1)
struktur_PhsSNR_r
