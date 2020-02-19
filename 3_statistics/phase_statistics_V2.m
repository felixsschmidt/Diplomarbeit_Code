% statistical analysis of phase
% 1. Data transformation 
% 2. Inferential statistics
% 3. Plotting
%% Set up
clear, clc, close all
format compact
cd 'C:\Users\Felix\Desktop\Close Loop\CMC\mat'
currentDir = pwd;
addpath(genpath(fullfile(currentDir, 'data')));
addpath(genpath(fullfile(currentDir, 'code')));
sb_figDir = 'C:\Users\Felix\Desktop\Close Loop\CMC\mat\data\figures\';
sb_resultDir = 'C:\Users\Felix\Desktop\Close Loop\CMC\mat\data\results';
%% 1. Transformation of values
load('phase_values_in_both_conditions.mat')
freqOfint = 13:35;
N = length(phase_values);

% get mean complex magnitudes of fourier coefficients in both conditions
MeanMagnitude = zeros(N,length(freqOfint)); 
for vp = 1:N     
    for freq = 1:length(freqOfint)
        % vec with indiv freq-specific mean values
        MeanMagnitude(vp,freq) =  mean(abs(phase_values{vp}.fourierCoef_H(freq+5,:)));
    end    
end

%%%%%%% (1) HIGH MEP %%%%%%%
% [rows:sub ID x columns:phase values_h]
MeanPhaseVal_h = zeros(N, length(freqOfint)); 

for vp = 1:N     
    for freq = 1:length(freqOfint)
        % vec with indiv freq-specific mean values
        MeanPhaseVal_h(vp,freq) =  phase_values{vp}.Mean_hMEP_phase_val(freq+5);
    end    
end

% [rows:freq x columns:freq-spec mean across subjects]
MeanAcrossVP_h = zeros(length(freqOfint),2);
MeanAcrossVP_h(:,1) = freqOfint;

for freq = 1:length(freqOfint)
    % vec with group freq-specific mean values
    MeanAcrossVP_h(freq,2) = circ_mean(MeanPhaseVal_h(:,freq));
end


%%%%%%% (1) LOW MEP %%%%%%%
% [rows:sub ID x columns:phase values_l]
MeanPhaseVal_l = zeros(N, length(freqOfint)); 

for vp = 1:N     
    for freq = 1:length(freqOfint)
        % vec with indiv freq-specific mean values
        MeanPhaseVal_l(vp,freq) =  phase_values{vp}.Mean_lMEP_phase_val(freq+5);
    end    
end

% [rows:freq x columns:freq-spec mean across subjects]
MeanAcrossVP_l = zeros(length(freqOfint),2);
MeanAcrossVP_l(:,1) = freqOfint;

for freq = 1:length(freqOfint)
    % vec with group freq-specific mean values
    MeanAcrossVP_l(freq,2) = circ_mean(MeanPhaseVal_l(:,freq));
end

%% 2. Inferential statistics

% alpha level is 0.05
% alpha level adjustment after Simes-Hochberg (2.4):
% highest p-value is handled the most liberally
% sequential testing of all comparisons, no break condition as in 
% Bonferroni-Holm (smalles p-values is handled most strictly)

% Hypotheses:
% H0i: mu2(phs(?)_highMEP) == mu1(phs(?)_lowMEP)
% H1i: mu2(phs(?)_highMEP) ~=  mu1(phs(?)_lowMEP) 

% structure of MeanPhaseVal_X: [rows:sub ID x columns:phase values]
%% 2.1 Parametric Watson-Williams multi-sample test for equal mean direction
wantOtherTests = 0;
if wantOtherTests == 1
% assumptions: data follow a von Mises distribution, concentration
% parameter kappa are equal in different samples

% circ_ktest: test for concentration parameter to be equal
    Phs_decisionS = zeros(length(freqOfint), 2); Phs_decisionS(:,1) = freqOfint;

    for freq = 1:length(freqOfint)

        [p table] = circ_wwtest(MeanPhaseVal_h(:,freq), MeanPhaseVal_l(:,freq));
        Phs_decisionS(freq,2) = p;

    end
%save(fullfile(sb_resultDir, 'Mean Phase values across subjects for both conditions.mat'), ...
%    'MeanPhaseVal_h', 'MeanPhaseVal_l');

% test outputs departures from assumptions: test might not be applicable
end
%% 2.2 nonparametric equal median test (Fisher NI, 1995)
if wantOtherTests == 1
    % H0: populations have equal medians
    Phs_decisionS21 = zeros(length(freqOfint), 4); Phs_decisionS21(:,1) = freqOfint;
    for freq = 1:length(freqOfint)

        [pval med P] = circ_cmtest(MeanPhaseVal_h(:,freq), MeanPhaseVal_l(:,freq));
        Phs_decisionS21(freq,2) = pval;
        Phs_decisionS21(freq,3) = med;
        Phs_decisionS21(freq,4) = P;

    end

% nothing is significant

end

%% 2.3 Non-parametric two-sample WATSON's U^2 test
if wantOtherTests == 1
    % # after the Zar implementation
    % https://de.mathworks.com/matlabcentral/fileexchange/43543-pierremegevand-watsons_u2

    % H0(i): two samples came from the same population/from two populations
    % having the same direction
    % H1(i): two samples do not come from the same population/from two populations
    % having the same direction

    % calculate critical values
    Phs_decisionS23 = zeros(length(freqOfint), 3); Phs_decisionS23(:,1) = freqOfint;
    U2_critV = 0.1857;
    for freq = 1:length(freqOfint)

         U2 = watsons_U2(MeanPhaseVal_h(:,freq), MeanPhaseVal_l(:,freq));
         Phs_decisionS23(freq,2) = U2; % outputs critical values 

    end

    % critical value for U2(alpha(0.05),n1(11),n2(11)) = 0.1857
    % Zar JH (1999). Biostatistical Analysis. 4th edition. Prentice Hill.

    % check significance with Bonferroni-Holm alpha adjustment
    if max(Phs_decisionS23(:,2)) < U2_critV
        fprintf(' 0 significant results \n')
    end
end
%% 2.4 Non-parametric two-sample WATSON's U^2 test
freqOfint = 13:35;
% # after the Zar implementation with p-value calculation
Phs_decisionS24 = zeros(length(freqOfint), 4); Phs_decisionS24(:,1) = freqOfint;
% calculate p-values
for freq = 1:length(freqOfint)
     [pval U2] = watsons_U2_approx_p(MeanPhaseVal_h(:,freq), MeanPhaseVal_l(:,freq));
     Phs_decisionS24(freq,2) = pval; % outputs p-values 
     Phs_decisionS24(freq,3) = U2;         
end

% check for significance with Simes-Hochberg procedure
sortedPhs_dec = sort(Phs_decisionS24(:,2), 'descend');
alpha = 0.05;
for ii = 0:length(sortedPhs_dec)-1
        if sortedPhs_dec(ii+1) <= alpha/(ii+1)
            Phs_decisionS24(ii+1,4) = 1;
        else
            Phs_decisionS24(ii+1,4) = 0;
        end       
end

fprintf(' %d significant results out of %d comparisons \n',...
    sum(Phs_decisionS24(:,4) > 0), length(freqOfint))   

save(fullfile(sb_resultDir, 'Watsons U^2 on frequency-wise phase differences between conditions.mat'),...
    'Phs_decisionS24'); % freq range 13-35 Hz
%% effect size
% difficult for U^2
%% 2.5 Phase-Coherence Correlation Analysis (Circular-Linear)
% Frage: Macht es Sinn den frequenzspez. Phasenunterschied zwischen den 
% beiden Bedingungen mit den korrespondierenden Kohärenzunterschieden zw.
% den Bedingungen zu korrelieren, wo es ja keine signifikante Phasenunter-
% schiede gibt? Macht es nicht vielleicht vielmehr Sinn, Phase in Bed. A
% mit CMC in Bed. A zu vergleichen (äquivalent in Bedingung B):
% Forschungsfrage dann: Hängen Phase und Kohärenz zusammen? Geht weniger
% auf die Unterschiede.
% Das gleiche Überlgegung gilt für Power-Analyse.

% Man könnte auch Phase mit CMC in ganzem 1s preTMS interval korrelieren,
% nur dann wären die Messzeitpunkte für Phase und CMC respektive sehr
% verschieden.
load('CMC_values_in_both_conditions.mat', 'Coh_values');

n = 11; % number of subjects
m = length(13:35); % number of frequencies [13-35 Hz]
k = 2; % number of conditions

% get differences between phase values and between CMC values
allVP_phsDiffCMC_v = zeros(n,2,m); % [abs phsDiff x abs CMCDiff] 
for vp = 1:n    
    for frq = 1:m
        allVP_phsDiffCMC_v(vp,1,frq) = MeanPhaseVal_h(vp,frq) - MeanPhaseVal_l(vp,frq);
        % should the phase vector contain the abs diff between phases?
        % I don't think so. Otherwise you wouldn't get the true difference
        allVP_phsDiffCMC_v(vp,2,frq) = rau(Coh_values{vp}.CohspctrmHigh(1,frq+5)) - ...
            rau(Coh_values{vp}.CohspctrmLow(1,frq+5)); 
    end
end

%%%%% OBS: %%%%%
% CMC values have been RAU transformed before their difference was taken

% RAU-transform CMC data for linear correlation with phase values
%allVP_phsDiffCMC_v(:,2,:) = rau(allVP_phsDiffCMC_v(:,2,:));

% get frequency-wise circ-linear correlation coefficients 
%  [freq x corrcoef x p-value x crit chi^2 value x decision]
Corco = zeros(m,5); Corco(:,1) = [13:35];n=length(allVP_phsDiffCMC_v(:,1,frq));
for frq = 1:m
    [rho pval] = circ_corrcl(allVP_phsDiffCMC_v(:,1,frq), allVP_phsDiffCMC_v(:,2,frq));
    Corco(frq,2) = rho;  
    Corco(frq,3) = pval;
    % n*r^2 can be compared with X^2(2):
    Corco(frq,4) = n*rho^2;
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

save(fullfile(sb_resultDir, 'Correlation Analysis on frequency-wise Phase vs. CMC differnces.mat'),...
    'Corco'); % beta range 13-35Hz
              % Phase values were not transformed
              
    
%% 3. Plotting
%% 3.0 scatterplot mean phases at spec freq + corresponding mean phases across subj

figure(1), subplot(1,2,1), 
for freq = 1:length(freqOfint)
    
    for vp=1:N
        p=polarplot( MeanPhaseVal_h(vp,freq), 1, 'ko'); hold on
        set(p,'MarkerSize',10)
    end
    p = polarplot( [0 MeanAcrossVP_h(freq,2)], [0 1], 'r' );
    set(p, 'linewidth', 1.2)
end
%title({'Individual frequency specific phase angles with respective average phases',
%'across subjects for these frequencies in high MEP trials'})
Ax = gca; % current axes
Ax.RGrid = 'off';
Ax.RTickLabel = []; 

subplot(1,2,2)
for freq = 1:length(freqOfint)
    
    for vp=1:N
        p=polarplot( MeanPhaseVal_l(vp,freq), 1, 'ko'); hold on
        set(p,'MarkerSize',10)
    end
    p = polarplot( [0 MeanAcrossVP_l(freq,2)], [0 1], 'r' );
    set(p, 'linewidth', 1.2)
end
%title({'Individual frequency specific phase angles with respective average phases',
%'across subjects for these frequencies in low MEP trials'})
set(gcf,'color','w');
Ax = gca; % current axes
Ax.RGrid = 'off';
Ax.RTickLabel = []; 

%saveas(figure(1), [sb_figDir '/matlab_fig/phase differences in both conditions.png'])

%circ_plot(MeanPhaseVal_h(vp,freq), 'pretty','bo',true,'linewidth',2,'color','r')
%% 3.1 plots mean phases at spec freq + corresponding mean phases across subj
wantOtherFig = 0;

if wantOtherFig == 1

    %figure(1), 
    subplot(2,2,3),
    for freq = 1:length(freqOfint)

        for vp=1:N
            p=polarplot([0 MeanPhaseVal_h(vp,freq)], [0 1], 'Color', [0.7 0.7 0.7],...
              'LineWidth',1); hold on
        end
        p = polarplot( [0 MeanAcrossVP_h(freq,2)], [0 1],'Color', [1 0 0],...
            'LineWidth',1.8);
    end
    %title({'Individual frequency specific phase angles with respective average phases',
    %'across subjects for these frequencies in high MEP trials'})
    Ax = gca; % current axes
    %Ax.ThetaGrid = 'off';
    Ax.RGrid = 'off';
    Ax.RTickLabel = []; 
    %Ax.ThetaTickLabel = [];

     subplot(2,2,4) 
    for freq = 1:length(freqOfint)

        for vp=1:N
            p=polarplot([0 MeanPhaseVal_l(vp,freq)], [0 1], 'Color', [0.7 0.7 0.7],...
              'LineWidth',1); hold on
        end
        p = polarplot( [0 MeanAcrossVP_l(freq,2)], [0 1],'Color', [1 0 0],...
            'LineWidth',1.5);
    end
    %title({'Individual frequency specific phase angles with respective average phases',
    %'across subjects for these frequencies in low MEP trials'})
    set(gcf,'color','w');

    Ax = gca; % current axes
    %Ax.ThetaGrid = 'off';
    Ax.RGrid = 'off';
    Ax.RTickLabel = []; 
    %Ax.ThetaTickLabel = [];
end
%% 3.2 plots freq spec mean phases across subjects

needed = 0; 
if needed == 1
    figure(2), subplot(1,2,1)
    for freq = 1:length(freqOfint)

        p = polarplot( [0 MeanAcrossVP_h(freq,2)], [0 1],...
            'LineWidth',1.5); 
        hold on
    end

    subplot(1,2,2)
    for freq = 1:length(freqOfint)

        p = polarplot( [0 MeanAcrossVP_l(freq,2)], [0 1],...
            'LineWidth',1.5); 
        hold on
    end
end




