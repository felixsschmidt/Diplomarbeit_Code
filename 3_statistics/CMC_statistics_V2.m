% statistical analysis of CMC
% 1.1 Descriptives
% 1.2 Continous plot
% 2.1 Pre-test to inferential statistics
% 2.2/3 Inferential statistics (Wilcoxon)

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

%% 1.1 Transform data and calculate + plot descriptive statistics
load('CMC_values_in_both_conditions.mat', 'Coh_values');
% CMC has only to be RAU transformed where it is compared with a different
% measurement (e.g. in correlation). Then, both value pairs need to have
% the same distribution

% create target vector [subNum x cmc high x cmc low]
cmcValues = zeros(length(Coh_values),3);
for ii=1:length(Coh_values)
    cmcValues(ii,1) = Coh_values{ii}.subNum;
    cmcValues(ii,2) = Coh_values{ii}.CohHighMEP;
    cmcValues(ii,3) = Coh_values{ii}.CohLowMEP;
end

%% 1.2 Plot continuous CMC in both conditions
% for every vpn, put cohspctrm [1:40] into one row of vector CohVec
% first 7 values are 0, as coh was only calculated for 8:40 Hz
frqLoopLength = 0:(length(Coh_values{1}.CohspctrmHigh)-1);
% 1 high MEP trials
CohVecH = zeros(11, 40);
for dd = 1:length(Coh_values)
    
    for ii = frqLoopLength
            CohVecH(dd,ii+8) = Coh_values{dd}.CohspctrmHigh(ii+1);
    end
    
end

% for every freq 8:40, calculate mean + sd
M_SD_high = zeros(2, length(CohVecH)); % vector containing mean + sd values
% std(x)/sqrt(length(x))
for ii = 1:length(CohVecH)
    M_SD_high(1,ii) = mean(CohVecH(:,ii));
    M_SD_high(2,ii) = std(CohVecH(:,ii))/sqrt(length(CohVecH(:,ii)));
    
end

%%%%%%%%

% 2 low MEP trials
CohVecL = zeros(11, 40);
for dd = 1:length(Coh_values)
    
    for ii = frqLoopLength
            CohVecL(dd,ii+8) = Coh_values{dd}.CohspctrmLow(ii+1);
    end
    
end

% for every freq 8:40, calculate mean + sd
M_SD_low = zeros(2, length(CohVecL)); % vector containing mean + sd values

for ii = 1:length(CohVecL)
    M_SD_low(1,ii) = mean(CohVecL(:,ii));
    M_SD_low(2,ii) = std(CohVecL(:,ii))/sqrt(length(CohVecL(:,ii)));
    
end

mH = M_SD_high(1,:); % Mean high MEP trials
eH = M_SD_high(2,:); % SEM high MEP trials
mL = M_SD_low(1,:);  % Mean low MEP trials
eL = M_SD_low(2,:);  % SEM high MEP trials
x = (1:size(mH,2));  % x-axis values

figure(1)
[h1 hp] = boundedline(x, mH, eH, 'b', x, mL, eL, 'r', 'alpha', 'transparency', 0.17);
hold on
ax = gca;
ax.FontSize = 12; 
legend('high MEP trials', 'low MEP trials', 'Location', 'northwest')
xlim([8 40]), xlabel('Frequency (Hz)', 'FontSize', 15), ylabel('Corticomuscular Coherence', 'FontSize', 15)
%title('Mean ± SEM of Coherence in both conditions')
set(h1,'LineWidth',1.5)
set(gcf,'color','w');
%saveas(figure(1), [sb_figDir '/matlab_fig/CMC spectra.png'])
%% 1.3 Descriptive plots
% Bar plot of individual CMC values in both conditions
barVec = zeros(length(cmcValues),2);
for ii=1:length(cmcValues)
    barVec(ii,1) = cmcValues(ii,2);
    barVec(ii,2) = cmcValues(ii,3);
end

max(barVec(:,2)) - min(barVec(:,2))

figure(2), subplot(2,2,1), b=bar(barVec); b(2).FaceColor = [0.8500 0.3250 0.0980];
b(1).FaceColor =[0 0.4470 0.7410];
legend('high MEP trials', 'low MEP trials', 'Location', 'northwest')
ylabel('Coherence'), xlabel('Subjects')
set(gca,'box','off','ZTick',[],'ZColor','w')

% Mean + SEM plot for both conditions
Mhigh = mean(cmcValues(:,2)); Mdhigh = median(cmcValues(:,2));
Mlow = mean(cmcValues(:,3)); Mdlow = median(cmcValues(:,3));
SDhigh = std(cmcValues(:,2));
SDlow = std(cmcValues(:,3));
SEM_high = SDhigh/sqrt(length(cmcValues(:,2))); % std(x)/sqrt(length(x))
SEM_low = SDlow/sqrt(length(cmcValues(:,3)));
%
Means = [Mhigh Mlow]; Mds = [Mdhigh Mdlow];
Stds = [SDhigh SDlow];
SEMs = [SEM_high SEM_low];

subplot(2,2,2), errorbar(Means, SEMs, 'ko', 'LineWidth', 1.5, 'MarkerSize', 2.2)
hold on
ylim([0 .3]), xlim([0 3]), ylabel('Coherence'),
set(gca,'XTick',1:2);
set(gca,'XTickLabel',{'High MEP trials', 'Low MEP trials' });
set(gca,'box','off','ZTick',[],'ZColor','w')
%title('Mean + SEM of CMC in both conditions')
% xlim([.6 2.4]

% Scatter plot of frequency distribution of CMC per subject
% plot distribution of frequency, in which subjects' CMC value peaks
frqOfCMC = zeros(11,2); frqOfCMC(:,1) = 1:11; %[subj x freq]
for ii = 1:length(frqOfCMC(:,1)) 
    frqOfCMC(ii,2) = Coh_values{ii}.cmcFrq;
end
subplot(2,2,3)
p= plot(frqOfCMC(:,1), frqOfCMC(:,2), '.', 'MarkerSize', 15);
xlim([0.5 11.5]), ylim([10 35]), xticks([1:11])
xlabel('Subjects'), ylabel({'Peak Frequency ' , 'of Coherence (Hz)'})
set(p,'color', [0 0 0]), set(gca,'box','off','ZTick',[],'ZColor','w')
%title('Subject specific peak frequency of Corticomuscular Coherence', 'Position', [6 36])

% Channel distribubtion of indiv peak CMC
Pk_e = zeros(length(Coh_values),1);
for vp = 1:length(Coh_values)    
    Pk_e(vp) = Coh_values{vp}.EOI;  
    if Pk_e(vp) == 17
        Pk_e(vp) = 16.2;
    end
end

subplot(2,2,4)
q=histogram(Pk_e,'BinWidth', 0.1, 'FaceColor', [0.9550, 0.0780, 0.001840]); 
xlim([16 16.3]), ylim([0 11]), yticks(1:11)
set(gca,'XTickLabel',{'', 'e16', '', 'e17', '', 'e18'})
xlabel('Electrode of highest Corticomuscular Coherence '), ylabel('Number of Subjects')
%title('Count of subjects with peak Coherence in electrode ''e16'' or ''e17'' ',...
%    'Position', [16.15 11.4])
set(gca,'box','off','ZTick',[],'ZColor','w')
set(gcf,'color','w');

%saveas(figure(2), [sb_figDir '/matlab_fig/Descriptive data visualization.jpg'])

%% Additional descriptives
find(cmcValues(:,2:3) > 0.2); % cmc values for both conditions exceeding 0.2
hh = find(cmcValues(:,2) > 0.2); % cmc values for high MEP trials exceeding 0.2
ll = find(cmcValues(:,3) > 0.2); % cmc values for low MEP trials exceeding 0.2
mean(frqOfCMC(:,2)); std(frqOfCMC(:,2)); range(frqOfCMC(:,2));
%% TOPO
wantTopo = 0;
if wantTopo == 1

    load DRCMR_EEG_Head.mat
    cfg = [];
    cfg.layout = lay;
    cfg.box = 'no';
    topo = ft_layoutplot(cfg);
    set(gcf,'color','w');
end

% save manually
%% 2. Inferential statistics
%% 2.1 Assess if cmc values are normally distributed
% put all values in one vector
allValues=zeros(length(cmcValues(:,2))*2,1);
allValues(1:length(cmcValues(:,2))) = cmcValues(:,2);
allValues(12:end) = cmcValues(:,3);
% descriptive: 
subplot(2,1,1), histCMC = histogram(allValues, 8); % Histogram
hold on, axis tight, xlabel('CMC'), ylabel('number of observations')
title('Histogramm of Coherence distribution') 
subplot(2,1,2), QQp = qqplot(allValues);           % QQ-Plot

% inferential: Shapiro-Wilk-test of normality
[H, pValue, W] = swtest(allValues, 0.05); % input: data, alpha
if H == 1
    fprintf('CMC values are not normally distributed \n')
end

% pValue: probability of observing the given
% result by chance given that the null hypothesis is true
% if pValue < alpha -> significant

%% 2.2 Wilcoxon signed-rank test

% has to be implemented: alpha-adjustment 

% question: are central tendencies in both dependent samples different?
 
% right sided: 
% H0: mu2(CMC_highMEP) <= mu1(CMC_lowMEP)
% H1: mu2(CMC_highMEP) >  mu1(CMC_lowMEP) 
[P,H,STATS] = signrank(cmcValues(:,2),cmcValues(:,3),'tail','right');

% right tailed test assumes that diffence of cmcValues(:,2) - cmcValues(:,3)
% (cmc high MEP vs cmc low MEP)comes from a distribution whose median is > 0

if H == 0
    fprintf('Non-significant - H0 cannot be rejected \n')
    fprintf('p-value: %f \n', P) 
else
    fprintf('Significant - H1 can be accepted \n')
    fprintf('p-value: %f \n', P) 
end

save(fullfile(sb_resultDir, 'Wilcoxon on peak CMC differences.mat'),...
    'P', 'H', 'STATS');

%% standardized effect size (Cohen, 1988 [pp. 48])
cmcDiffValues = cmcValues(:,2) - cmcValues(:,3);
Mdiff = mean(cmcDiffValues);
SDdiff = std(cmcDiffValues);
%SDdiff = std(cmcValues(:,2)) - std(cmcValues(:,3));

dz = Mdiff/std(cmcDiffValues);

%% Cohen's d effect size
M1=mean(cmcValues(:,2));M2=mean(cmcValues(:,3));
SD1=std(cmcValues(:,2));SD2=std(cmcValues(:,3));
[M1 M2; SD1 SD2];
corrcoef(cmcValues(:,2), cmcValues(:,3))

Mdiff/ mean([SD1 SD2]);
%% 2.3 Wilcoxon manual
N = 11;
wxVec = zeros(length(cmcValues), 8);

% [1.ID 2.CMC_h 3.CMC_l 4.Diff 5.Sign 6.Rank]

for ii = 1:length(cmcValues)
    wxVec(ii,1) = cmcValues(ii,1);
    wxVec(ii,2) = cmcValues(ii,2);
    wxVec(ii,3) = cmcValues(ii,3);
    wxVec(ii,4) = abs(cmcValues(ii,2)-cmcValues(ii,3));
    if (cmcValues(ii,2)-cmcValues(ii,3)) < 0
        wxVec(ii,5) = -1;
    elseif (cmcValues(ii,2)-cmcValues(ii,3)) > 0
        wxVec(ii,5) = 1;
    else
        wxVec(ii,5) = 9999;
    end
end

% create digit positions in gst(:,2), gst(:,1) will be nonsense
gst = zeros(length(wxVec), 2);
gst(:,1) = wxVec(:,4);
for ii = 1:length(cmcValues)
    [val id] = min(gst(:,1));
    gst(id,2) = ii;
    gst(id,1) = 99+ii;
end

wxVec(:,6) = gst(:,2);
 
% [1.ID 2.CMC_h 3.CMC_l 4.Diff 5.Sign 6.Rank 7.neg Rank 8.pos Rank]
for ii = 1:length(cmcValues)
    if wxVec(ii,5) < 0
       wxVec(ii,7) = wxVec(ii,6);
    else
        wxVec(ii,8) = wxVec(ii,6);
    end
end

T_minus = sum(wxVec(:,7)); T_plus = sum(wxVec(:,8));

W = min([T_minus T_plus]); % if W <= Wcrit -> significant
Wcrit = 13; % for N=11 and alpha=.05, taken from table

if W > Wcrit
    fprintf('non-significant differences between CMC values \n')
end

