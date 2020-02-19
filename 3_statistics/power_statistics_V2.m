% statistical analysis of power
% 1. Data transformation + plotting
% 2. Inferential statistics
%% Set up
clear, clc, close all
format compact
cd 'C:\Users\Felix\Desktop\Close Loop\CMC\mat'
currentDir = pwd;
addpath(genpath(fullfile(currentDir, 'data')));
addpath(genpath(fullfile(currentDir, 'code')));
addpath(genpath(fullfile(currentDir, 'kakearney-boundedline-pkg-50f7e4b')));
load DRCMR_EEG_Head.mat
sb_figDir = 'C:\Users\Felix\Desktop\Close Loop\CMC\mat\data\figures\';
sb_resultDir = 'C:\Users\Felix\Desktop\Close Loop\CMC\mat\data\results';

% check if descriptive correlation plots are desired
wantCorrPlots = 1;
load(fullfile(sb_resultDir, 'Correlation Analysis on frequency-wise Phase vs. CMC differnces.mat'),...
    'Corco')
CorcoPhs=Corco;

%% 1. Transformation of values and figure
load('pow_values_in_both_conditions.mat')
loopLength = length(power_values{1}.highMEP_pow_val);
H_MStdpowerV = zeros(3,loopLength);
ghostV = zeros(1,11);

%%%%%% high MEP trials %%%%%%
% form: [M x SEM x Freq]'
for ii = 1:loopLength
         
    for jj = 1:11         
        ghostV(ii,jj) = power_values{jj}.highMEP_pow_val(1,ii);
    end
    
    % first row: grand mean
    H_MStdpowerV(1,ii) =   mean(ghostV(ii,:));
    % second row: SEM   
    H_MStdpowerV(2,ii) =   std(ghostV(ii,:))/sqrt(loopLength);   
    
end
% third row: frequencies
H_MStdpowerV(3,:) = power_values{1}.highMEP_pow_val(2,:); 


%%%%%% low MEP trials %%%%%%
L_MStdpowerV = zeros(3,loopLength);
ghostV = zeros(1,11); % stores power values of freq x temporarily

% high MEP trials
% form: [M x SEM x Freq]'
for ii = 1:loopLength
         
    for jj = 1:11         
        ghostV(ii,jj) = power_values{jj}.lowMEP_pow_val(1,ii);
    end
    
    % first row: grand mean
    L_MStdpowerV(1,ii) =   mean(ghostV(ii,:));
    % second row: SEM   
    L_MStdpowerV(2,ii) =   std(ghostV(ii,:))/sqrt(loopLength);   
    
end
% third row: frequencies
L_MStdpowerV(3,:) = power_values{1}.lowMEP_pow_val(2,:);  

% plot
mH = H_MStdpowerV(1,:);
eH = H_MStdpowerV(2,:);
mL = L_MStdpowerV(1,:);
eL = L_MStdpowerV(2,:);
x = H_MStdpowerV(3,:);

figure(1)
[h1 hp] = boundedline(x, mH, eH, 'b', x, mL, eL, 'r', 'alpha', 'transparency', 0.17);
hold on
h = legend('high MEP trials', 'low MEP trials');
rect = [.2, .8, .2 .07 ];
set(h, 'Position', rect);
ax = gca;
ax.FontSize = 11.5; 
axis tight, xlabel('Frequency (Hz)', 'FontSize', 14), ylabel('Power (\muV^2)', 'FontSize', 14) % arbitrary
%title('Power spectra from sensorimotor electrodes in high vs. low MEP trials')
set(h1,'LineWidth',1.3)
set(gca,'box','off','ZTick',[],'ZColor','w')
set(gcf,'color','w');
%saveas(figure(1), [sb_figDir '/matlab_fig/Power spectra.png'])

% unresolved issues:
% power umrechnen: y-Axen values sind sehr klein

%% 2. Topo Power
%% 2.1 Prepping data
load Topo_power_Val.mat

frqOfint = 8:40; % vector with analyzed frequencies
twLength = zeros(1,length(frqOfint));
ind = zeros(1,length(Hpow.time));
power_valuesH = zeros(2,length(Hpow.freq),3);

for ii = 1:length(Hpow.freq)
    
    for jj = 1:63
       % vector with time window (tw) length (3 cycles) of respect freq
       twLength(ii) = 3/frqOfint(ii); 
       % output indices from time variable,
       % where tw center + tw/2 dont exceed stimulus onset  
       ind = (Hpow.time + (twLength(ii)/2)) <= (-0.019999); 
       % pick center point from tw, where it is closest to stimulus onset
       [val idx] = min(abs(Hpow.time(ind)));
       % output power value for that TF-sample closest to trigger
       power_valuesH(1,ii,jj) = Hpow.powspctrm(jj,ii,idx); % safe power value
       power_valuesH(2,ii,jj) = frqOfint(ii); % corresponding frequency
    end
   
end

power_valuesL = zeros(2,length(Hpow.freq),3);
for ii = 1:length(Lpow.freq)
    
    for jj = 1:63
       % vector with time window (tw) length (3 cycles) of respect freq
       twLength(ii) = 3/frqOfint(ii); 
       % output indices from time variable,
       % where tw center + tw/2 dont exceed stimulus onset  
       ind = (Lpow.time + (twLength(ii)/2)) <= (-0.019999); 
       % pick center point from tw, where it is closest to stimulus onset
       [val idx] = min(abs(Lpow.time(ind)));
       % output power value for that TF-sample closest to trigger
       power_valuesL(1,ii,jj) = Lpow.powspctrm(jj,ii,idx); % safe power value
       power_valuesL(2,ii,jj) = frqOfint(ii); % corresponding frequency
    end
   
end

%% 2.2 Plotting

% select time-freq block closest to stimulus onset
    frqOfint = 8:40; % vector with analyzed frequencies
    
    % vector with time window (tw) length (3 cycles)
    % of lowest freq (8 Hz), as this will have the largest tw of all
    % frequencies and therefore the most negative center time.
    % It's spanning the max time frame right before SO.
    
    twLength_8Hz = 3/frqOfint(1); 
    twLength_13Hz = 3/frqOfint(6); 
    % output indices from time variable,
    % where tw center + tw/2 dont exceed stimulus onset 
    % irrelevant of condition
    ind = (Lpow.time + (twLength_13Hz/2)) <= (-0.019999); 
    % pick center point from tw, where it is closest to stimulus onset
    Mintim = min(abs(Lpow.time(ind))); Mintim = -Mintim;
    Maxtim = Mintim + (twLength_13Hz/2);
        
% subtract the two power spectra and then divide them by their sum 
% - this normalizes the difference by the common activity
    cfg = [];
    cfg.parameter = 'powspctrm';
    cfg.operation = '(x1-x2)/(x1+x2)';
    tfr_difference = ft_math(cfg, Hpow, Lpow);

% plot difference
% not baseline corrected
    figure(2)
    cfg = [];
    cfg.xlim = [Mintim Maxtim]; % TF-block closest to SO [-0.22 -0.0325]
    cfg.ylim = [13 35]; % freq dimension
    cfg.zlim = 'maxabs';
    cfg.layout = lay;
    cfg.highlight = 'on';
    cfg.highlightchannel = {'e16'; 'e17'; 'e18'};
    cfg.highlightsymbol = '.';
    cfg.highlightsize = 20;
    cfg.highlightcolor = [1 0 0];
    cfg.comment = 'no';
    ft_topoplotTFR(cfg, tfr_difference);
    ax = gca;
    ax.FontSize = 20; 
    c = colorbar; y = ylabel(c, 'Power (\muV^2)', 'FontSize',20);
    set(c,'Location', 'southoutside') % xposition yposition width height
    set(gcf,'color','w');
    
%saveas(figure(2), [sb_figDir '/matlab_fig/Topoplot Power differences_V2.png'])
 
    
%% singleplot at e11
single = 0;
if single == 1
    cfg = [];
    cfg.parameter = 'powspctrm';
    cfg.xlim = [Mintim Maxtim]; % time
    cfg.zlim = 'maxabs';
    cfg.ylim = [13 35]; % freq dimension
    cfg.channel = {'e11'};
    cfg.colorbar = 'yes';
    cfg.colorbartext = '(\muV^2)';
    figure; ft_singleplotTFR(cfg, tfr_difference);
end

% [power_valuesH(:,15,11), power_valuesL(:,15,11)]

%% 2. Inferential statistics
% create vector with power values in both conditions for each freq
n = 11; % number of subjects
m = length(power_values{1}.highMEP_pow_val(1,:)); % number of frequencies
k = 2; % number of conditions
allVP_pow_v = zeros(n,k,m); % has dimensions: (vp, cond, freq)

for vp = 1:n
    
    % save power values for high MEP trials
    for frq = 1:m
        allVP_pow_v(vp,1,frq) = power_values{vp}.highMEP_pow_val(1,frq);
    end
    
    % save power values for low MEP trials
    for frq = 1:m
        allVP_pow_v(vp,2,frq) = power_values{vp}.lowMEP_pow_val(1,frq);
    end
    
end
    
%save(fullfile(sb_resultDir, 'PowerValues_allVP_allFrq_allCond.mat'),'allVP_pow_v')
%load(fullfile(sb_resultDir, 'PowerValues_allVP_allFrq_allCond.mat'),'allVP_pow_v')

%% 2.1 Wilcoxon signrank test
% test if frequency-wise power in beta band [13:35] Hz differs between conditions 
% Simes-Hochberg alpha adjustment: highest p-value is handled the most
% liberally; sequential testing of all comparisons

% two sided wilcoxon sign rank tests: 
% H0i: mu2(pow(?)_highMEP) == mu1(pow(?)_lowMEP)
% H1i: mu2(pow(?)_highMEP) ~=  mu1(pow(?)_lowMEP) 

% vector to save freq that has been analysed, p-value and decision value
% [freq x p-value x decision (0='no') x T x alpha level applied at this p]
decisionS = zeros(length(13:35), 3); decisionS(:,1) = 13:35;
alpha = 0.05;

% calculate p-values,
% test checks significance with Bonferroni-Holm alpha adjustment
% automatically.
% test statistic has to be >= the critical value
for frq = 1:length(13:35)
    
    [P,H,STATS] = signrank(allVP_pow_v(:,1,frq + 5), allVP_pow_v(:,2,frq + 5),...
        'tail','left', 'alpha', alpha); 
    decisionS(frq,2) = P;
    decisionS(frq,3) = H;
    decisionS(frq,4) = STATS.signedrank;
end

% check for significance with Simes-Hochberg procedure
sorted_dec = sort(decisionS(:,2), 'descend');
for ii = 0:length(sorted_dec)-1
    decisionS(ii+1,5) = alpha/(ii+1);
        if sorted_dec(ii+1) <= alpha/(ii+1)
            decisionS(ii+1,3) = 1;            
        else
            decisionS(ii+1,3) = 0;
        end       
end

struktur = {'Frequenz', 'P-Wert', 'Entscheidung', 'T-Wert', 'alpha'};
fprintf('%d significant results out of %d comparisons \n',...
    sum(decisionS(:,3) > 0), length(13:35) )   

save(fullfile(sb_resultDir, 'Wilcoxon on frequency-wise power differences between conditions.mat'),...
    'decisionS', 'struktur'); % beta range 13-35 Hz

%% standardized effect size (Cohen, 1988 [pp. 48])
% compute d for the frequency, at which the largest mean power differences
% occur between conditions

% find out largest M1-M2 difference
M1_M2 = [13:35]'; % vector of mean power differences
for frq = 1:length(13:35)
    M1_M2(frq,2) = mean(allVP_pow_v(:,1,frq + 5)) - mean(allVP_pow_v(:,2,frq + 5));      
end
[v i] = max(M1_M2(:,2));
[decisionS(:,1:2), M1_M2(:,2)];

%powDiffValues = allVP_pow_v(:,1,1 + 5) - allVP_pow_v(:,1,1 + 5); % 13 Hz
Mdiff = M1_M2(i,2);
SDdiff = std(allVP_pow_v(:,1,frq + 7) - allVP_pow_v(:,2,frq + 7));
%mean(allVP_pow_v(:,1,1 + 7)) - mean(allVP_pow_v(:,2,1 + 7))
%std(allVP_pow_v(:,1,frq + 7)), std(allVP_pow_v(:,2,frq + 7))
%corrcoef(allVP_pow_v(:,1,1 + 7), allVP_pow_v(:,2,1 + 7))

dz = Mdiff/SDdiff;

%% 2.2 Correlation analysis of freq wise power differences between conditions 
% with CMC values in corresponding frequencies

% correlate differences in beta power between conditions 
% with differences in CMC between conditions for [13:35] Hz
load('CMC_values_in_both_conditions.mat', 'Coh_values');

% get differences between power values and between CMC values
m = length(13:35);
allVP_powDiffCMC_v = zeros(n,2,m); % [abs powDiff x abs CMCDiff] 
for vp = 1:n    
    for frq = 1:m
        allVP_powDiffCMC_v(vp,1,frq) = allVP_pow_v(vp,1,frq+5) - allVP_pow_v(vp,2,frq+5);
        allVP_powDiffCMC_v(vp,2,frq) = rau(Coh_values{vp}.CohspctrmHigh(1,frq+5)) - ...
            rau(Coh_values{vp}.CohspctrmLow(1,frq+5)); 
    end
end

%%%%% OBS: %%%%%
% CMC values have been RAU transformed before their difference was taken

% RAU-transform CMC data for linear correlation with power values
% allVP_powDiffCMC_v(:,2,:) = rau(allVP_powDiffCMC_v(:,2,:));

% get frequency-wise correlation coefficients (un-normalized, uncorrected)
% [freq x corrcoef x p-value x crit t value x decision]
Corco = zeros(m,4); Corco(:,1) = [13:35];n=length(allVP_powDiffCMC_v(:,1,1));
for frq = 1:m
    [r p] = corrcoef(allVP_powDiffCMC_v(:,1,frq), allVP_powDiffCMC_v(:,2,frq));
    Corco(frq,2) = r(1,2);  
    Corco(frq,3) = p(1,2);
    % calculate value of t-statistic (Zar, 2010 [pp.383]):
    Corco(frq,4) = r(1,2) / ((1-r(1,2)^2)/(n-2));
end

% test for significance with Simes-Hochberg adjustment
sorted_CorPval = sort(Corco(:,3), 'descend');
alpha = 0.05;
for ii = 0:length(sorted_CorPval)-1
        if sorted_CorPval(ii+1) <= alpha/(ii+1)
            Corco(ii+1,5) = 1;
        else
            Corco(ii+1,5) = 0;
        end       
end
fprintf(' %d significant results out of %d comparisons \n',...
    sum(Corco(:,5) > 0), length(Corco(:,1)))

save(fullfile(sb_resultDir, 'Correlation Analysis of frequency-wise Power vs. CMC differnces.mat'),...
    'Corco'); % beta range 13-15 Hz
              % not Fisher-z-transformed

%figure(2), plot(allVP_powDiffCMC_v(:,1,14), allVP_powDiffCMC_v(:,2,14), '.',...
%    'MarkerSize', 10)
%% Further descriptives for Corr analysis for power AND phase
% power
 CorcoPow=Corco;
 LnegCorr=length(find(CorcoPow(:,2) < 0));
 LposCorr=length(find(CorcoPow(:,2) > 0));

if wantCorrPlots == 1
    % power
    % figure,histogram(Corco(:,2)) unused
    figure(3),subplot(1,2,1)
    b=bar(CorcoPow(:,2),0.8); b.FaceColor =[1 1 1];
    set(gca, 'XTick', [1:2:23],'XTickLabel',[1:2:23]+12)
    ax = gca;
    ax.FontSize = 13;
    ylabel('Pearson''s Linear Correlation (\rho)', 'FontSize', 15),
    xlabel('Frequencies (Hz)', 'FontSize', 15)
    set(gca,'box','off','ZTick',[],'ZColor','w'), %axis tight
    %title('A', 'Position', [-0.1 0.45])

    % phase
    subplot(1,2,2)
    a=bar(CorcoPhs(:,2),0.8); a.FaceColor =[1 1 1];
    set(gca, 'XTick', [1:2:23],'XTickLabel',[1:2:23]+12)
    ax = gca;
    ax.FontSize = 13;
    ylabel('Angular-Linear Correlation (\rho_A_L)', 'FontSize', 15),
    xlabel('Frequencies (Hz)', 'FontSize', 15)
    set(gca,'box','off','ZTick',[],'ZColor','w')
    set(gcf,'color','w');
end

%saveas(figure(3), [sb_figDir '/matlab_fig/correlation plots.png'])

%% %%%%%%%%% UN-USED %%%%%%%%%%

% correlate high MEP trial power with CMC fluctuations in same conditions
% correlate low MEP trial power with CMC fluctuations in same conditions
allVP_powDiffCMC_v = allVP_pow_v; % [pow High x pow Low x CMC High x CMC Low]
for vp = 1:n
    
    for frq = 1:m
        allVP_powCMC_v(vp,3,frq) = Coh_values{vp}.CohspctrmHigh(1,frq);
        allVP_powCMC_v(vp,4,frq) = Coh_values{vp}.CohspctrmLow(1,frq);
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Correlation coefficient transformation
% Fisher führt eine asymptotische Normalisierung durch, wobei die Verteilung
% der Korrelationskoeffizienten approximativ in eine Normalverteilung überführt wird
