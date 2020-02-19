% statistical analysis of CMC in pre vs. post stim condition
% 1.  Descriptives + plotting
% 2.  Inferential statistics
% 2.1 Wilcoxon signed-rank-test on peak CMC differences
% 2.2 Wilcoxon signed-rank-test on freq-wise CMC differences
% 3   AMT/RMT pre post comparison
% 4   MVC pre post comparison

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

%% 1. Descriptives
load 'CMC_values_in_pre_and_post_condition.mat';

freqWiseDiff = 1;

% Plot mean and sem in both conditions
% create target vector [subNum x cmc pre x cmc post x foi]
n = length(Coh_values_nonStim);
cmcValues_pp = zeros(n,4);
for frq=1:n
    cmcValues_pp(frq,1) = Coh_values_nonStim{frq}.subNum;
    cmcValues_pp(frq,2) = Coh_values_nonStim{frq}.CohPreStim;
    cmcValues_pp(frq,3) = Coh_values_nonStim{frq}.CohPostStim;
    cmcValues_pp(frq,4) = Coh_values_nonStim{frq}.cmcFrq;

end
% descriptive values
Mpre = mean(cmcValues_pp(:,2));
SEMpre = std(cmcValues_pp(:,2)) / sqrt(n);
Mpost = mean(cmcValues_pp(:,3));
SEMpost = std(cmcValues_pp(:,3)) / sqrt(n);

M = [Mpre Mpost]; 
SEM = [SEMpre SEMpost];

figure(1), 
if freqWiseDiff ~= 1
    subplot(1,2,2)    
    errorbar(M, SEM, 'ko', 'LineWidth', 1.5, 'MarkerSize', 3)
    xlim([0 3]), ylim([0 0.25]), ylabel('Coherence')
    %xlabel('Mean ± SEM of individual peak Frequency of Coherence')
    set(gca,'XTick',1:2);
    set(gca,'XTickLabel',{'PreStim', 'PostStim' });
    set(gca,'box','off','ZTick',[],'ZColor','w')
end
% Continuous Plot
% 1 PreStim trials: freq specific Coh values across subjects
frqLoopLength = 0:(length(Coh_values_nonStim{1}.CohspctrmPre)-1);
CohVecPre = zeros(11, 40);
for vp = 1:n
    
    for frq = frqLoopLength
            CohVecPre(vp,frq+8) = Coh_values_nonStim{vp}.CohspctrmPre(frq+1);
    end
    
end

% for every freq 8:40, calculate mean + sd
M_SD_pre = zeros(2, length(CohVecPre)); % vector containing mean + sem values
for frq = 1:length(CohVecPre)
    M_SD_pre(1,frq) = mean(CohVecPre(:,frq));
    M_SD_pre(2,frq) = std(CohVecPre(:,frq))/sqrt(n);
    
end

% 2 PostStim trials: freq specific Coh values across subjects
CohVecPost = zeros(11, 40);
for vp = 1:n
    
    for frq = frqLoopLength
            CohVecPost(vp,frq+8) = Coh_values_nonStim{vp}.CohspctrmPost(frq+1);
    end
    
end

% for every freq 8:40, calculate mean + sd
M_SD_post = zeros(2, length(CohVecPost)); % vector containing mean + sem values
for frq = 1:length(CohVecPost)
    M_SD_post(1,frq) = mean(CohVecPost(:,frq));
    M_SD_post(2,frq) = std(CohVecPost(:,frq))/sqrt(n);
    
end

mPre = M_SD_pre(1,:); % Mean PreStim trials
ePre = M_SD_pre(2,:); % SEM PreStim trials
mPost = M_SD_post(1,:);  % Mean PostStim trials
ePost = M_SD_post(2,:);  % SEM PostStim trials
x = (1:size(mPre,2));  % x-axis values

if freqWiseDiff ~= 1
    subplot(1,2,1)
end

[h1 hp] = boundedline(x, mPre, ePre, 'b', x, mPost, ePost, 'r',...
    'cmap', [0 0 0; 0.6 0.2 0.2],'alpha', 'transparency', 0.17);
hold on
ax = gca;
ax.FontSize = 12; 
legend('PreStim trials', 'PostStim trials', 'Location', 'northwest')
xlim([8 40]), xlabel('Frequency (Hz)', 'FontSize', 15), 
ylabel('Corticomuscular Coherence', 'FontSize', 15)
%title('Mean +- SEM of Coherence in both conditions')
set(h1,'LineWidth',1.5)
set(gca,'box','off','ZTick',[],'ZColor','w')
set(gcf,'color','w');

%saveas(figure(1), [sb_figDir '/matlab_fig/CMC spectra_PrePostStim.png'])

%% 2. Inferential statistics
%% 2.1 Wilcoxon signed-rank-test
% test for normality
cmcValues_pp(:,2:3), allValues = ans(:);
[H, pValue, W] = swtest(allValues, 0.05);

if H == 0
    fprintf('CMC values are not normally distributed \n')
end

if freqWiseDiff ~= 1
    % right sided: 
    % H0: mu2(CMC_Pre) <= mu1(CMC_Pre)
    % H1: mu2(CMC_Pre) >  mu1(CMC_Pre) 
    [p, h] = signrank(cmcValues_pp(:,2),cmcValues_pp(:,3),'tail','right');
    if h == 0
        fprintf('Wilcoxon: Non-significant - H0 cannot be rejected \n')
        fprintf('p-value: %f \n', p) 
    else
        fprintf('Significant - H1 can be accepted \n')
        fprintf('p-value: %f \n', p) 
    end

    save(fullfile(sb_resultDir, 'Wilcoxon on peak CMC comparison between pre- and postStim.mat'),...
        'p', 'h');
end

%% 2.2 Wilcoxon signed-rank-test on freq-wise CMC differences
freq = 8:40;
decisionPP = zeros(length(freq), 4); decisionPP(:,1) = freq;
% compute p-values 
% two sided:
% H0: mu2(CMC_Pre(?)) >= mu1(CMC_Post(?))
% H1: mu2(CMC_Pre(?)) <=  mu1(CMC_Post(?)) 
% [frq x p-val x decision x test statistic value T x alpha level applied at this p ]
struktur = {'Frequenz', 'P-Wert', 'Entscheidung', 'T-Wert', 'alpha'};
for frq = 8:length(freq)+7
    
    [p, h, T] = signrank(CohVecPre(:,frq),CohVecPost(:,frq), 'tail', 'left');
    decisionPP(frq-7,2) = p;
    decisionPP(frq-7,3) = h;
    decisionPP(frq-7,4) = T.signedrank;
end
    
% check for significance with Simes-Hochberg procedure
sorted_dec = sort(decisionPP(:,2), 'descend');
alpha = 0.05;
for ii = 0:length(sorted_dec)-1
    decisionPP(ii+1,5) = alpha/(ii+1);
        if sorted_dec(ii+1) <= alpha/(ii+1)
            decisionPP(ii+1,3) = 1;
        else
            decisionPP(ii+1,3) = 0;
        end       
end

fprintf('%d significant results out of %d comparisons \n',...
    sum(decisionPP(:,3) > 0), length(freq)) 
fprintf('CAVE: very little p-values for frequencies 27, 28, 29 \n')

save(fullfile(sb_resultDir, 'Wilcoxon on freq-wise CMC comparisons between pre- and postStim.mat'),...
    'decisionPP', 'struktur');

%% standardized effect size (Cohen, 1988 [pp. 48])
% compute d for the frequency, at which the largest mean CMC differences
% occur between conditions

% find out largest M1-M2 difference
M1_M2 = [8:40]'; % vector of mean CMC differences
for frq = 8:length(M1_M2)+7 % 8:length(freq)+7
    M1_M2(frq-7,2) = mean(CohVecPre(:,frq)) - mean(CohVecPost(:,frq));      
end
[v i] = max(abs(M1_M2(:,2)));
fprintf('largest mean difference between pre and post CMC at %d Hz \n',i+7)

%powDiffValues = allVP_pow_v(:,1,1 + 5) - allVP_pow_v(:,1,1 + 5); % 13 Hz
Mdiff = abs(M1_M2(i,2));
SDdiff = std(CohVecPre(:,i) - CohVecPost(:,i));
%mean(allVP_pow_v(:,1,1 + 7)) - mean(allVP_pow_v(:,2,1 + 7))
%std(allVP_pow_v(:,1,frq + 7)), std(allVP_pow_v(:,2,frq + 7))
%corrcoef(allVP_pow_v(:,1,1 + 7), allVP_pow_v(:,2,1 + 7))

dz = Mdiff/SDdiff;

%% 3 RMT pre post comparison
% values were just copied out of excel sheet

preRMT = [57 39 39 79 46 39 42 48 54 52 68]';
postRMT = [62 44 42 72 42 42 45 49 58 57 68]'; prepostR = [preRMT; postRMT];

% Shapiro-Wilk-test for normality
[H, pValue, W] = swtest(preRMT- postRMT, 0.05); % input: data, alpha
if H == 1
    fprintf('RMTs are normally distributed on 5 %% alpha level \n')
end

% H0: Pairwise differences had a mean smaller or equal to zero
% HA: Pairwise differences had a mean larger than zero 
[h1,p1,CI,STATS] = ttest(preRMT, postRMT, 'tail', 'right');

if h1 == 0
    fprintf('Differences in motor thresholds is not significant \n')
else
    fprintf('Differences in motor thresholds is significant \n')
end

figure(2),subplot(1,2,1)
errorbar([mean(preRMT) mean(postRMT)], [std(preRMT) std(postRMT)], 'k.',...
    'LineWidth', 1.5, 'MarkerSize', 15)
xlim([0.5 2.8]),  ylim([25 75])
set(gca, 'XTick', [1:2],'XTickLabel',{'preStim RMT', 'postStim RMT'})
xlabel('Resting Motor Thresholds', 'Position', [1.5 21.5]),
ylabel('Stimulation intensity (% MSO)')
set(gca,'box','off') %,'ZTick',[],'ZColor','w')
set(gcf,'color','w')

% extras
%preAMT = [42 30 32 64 37 26 28 37 41 44 45]';
%postAMT = [40 30 34 52 37 30 29 32 47 44 59]'; prepostA = [preAMT postAMT];

%% 4 Pre- post-stim MVC level
% values were just copied out of excel sheet
preMVC = [55.3 59 37 49 43.5 81.4 58.7 72.3 59 42.8 56.3];
postMVC = [55.4 64.8 27 42.8 34.9 81.1 54.4 69.1 60.6 41.1 60.1];
prePostMVC=[preMVC postMVC]';

subplot(1,2,2),
errorbar([median(preMVC) median(postMVC)], [std(preMVC) std(postMVC)], 'k.',...
    'LineWidth', 1.5, 'MarkerSize', 15)
xlim([0.5 2.5]), ylim([20 80])
set(gca, 'XTick', [1:2],'XTickLabel',{'preStim MVC', 'postStim MVC'})
xlabel('Maximal Voluntary Contraction', 'Position', [1.5 15.5])
ylabel('Force level (N)')
%title('MVC force output (N) prior and after TMS intervention')
set(gca,'box','off') %,'ZTick',[],'ZColor','w')
set(gcf,'color','w')

% Shapiro-Wilk-test for normality
[H, pValue, STAT] = swtest([preMVC postMVC]', 0.1); % input: data, alpha
if H == 1
    fprintf('MVC levels are normally distributed on 10 %% alpha level \n')
end

zPrePostMVC = atanh(prePostMVC);
% H0: Pairwise differences had a median smaller or equal to zero
% HA: Pairwise differences had a median larger than zero 
[p2, h2, T] = signrank(preMVC, postMVC, 'tail', 'right');

if h2 == 0
    fprintf('Difference in MVC force levels were not significant \n')
else
    fprintf('Difference in MVC force levels were significant \n')
end

save(fullfile(sb_resultDir, 'PrePost threshold and MVC comparisons.mat'),...
    'h1', 'p1', 'h2', 'p2');

%saveas(figure(2), [sb_figDir '/matlab_fig/RMT and MVC distribution.jpg'])

