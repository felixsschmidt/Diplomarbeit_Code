%%      Analysis Pipeline for Corticomuscular Coherence
% 1.    Define 1s preTMS period where CMC is computed
% 2.    Define cell for channelcombination
% 3.1.1 Compute coherence over all trials
% 3.1.2 Define individual EOI and FOI
% 3.2   Compute Coherence for trials, corresponding to MEPs in both
%       conditions
%% Set up
clear, clc, close all
cd '/mnt/projects/CMCloop/data/random_stimulation_experiment/FT_Analysis'
setup_path; 
addpath(genpath(fullfile(projectdir, 'data'))); 
% Folder information
iSub = 11;                           % adjust
sb_preStimDir = fullfile(baseDir,indiDir,IDs{iSub},preprocDir, preStimDir);
sb_resultDir = fullfile(baseDir,indiDir,IDs{iSub}, resultDir);
sb_groupDir = fullfile(baseDir, groupDir);
%% 1. Define 1s preTMS period [-1020 -.02] ms
load(fullfile(sb_preStimDir, [num2str(IDs{iSub}) '_preTMS_EEG_EMG.mat']));

cfg = [];
cfg.toilim = [-1.02 -0.02];
data_preTMS = ft_redefinetrial(cfg, data_preproc);

% optional: plot to confirm
disp_to_comp = 0;
if disp_to_comp == 1
    
    figure
    subplot(2,1,1)
    plot(data_preTMS.time{1}, data_preTMS.trial{1}(17,:));
    axis tight;
    legend(data_preTMS.label(17));
    subplot(2,1,2)
    plot(data_preproc.time{1},data_preproc.trial{1}(17,:));
    axis tight;
    legend(data_preproc.label(17));
end

%% 2. define Nx2 cell for channelcmb
% combine all EEG channel with EMGright
ccomb = cell(length(data_preTMS.label)-1,2);
for ii=1:length(data_preTMS.label)-1
    ccomb{ii,1} = char(lay.label(ii));
    ccomb{ii,2} = 'EMGright';       
end

%% 3.1.1 Compute cross- & power & cmc spectra over all trials
% cross- & power spectra
cfg = [];
cfg.output = 'powandcsd';
cfg.method = 'mtmfft'; % freq analysis or tf-analysis using 'mtmconvol'?
cfg.foi = [8:1:40];
cfg.tapsmofrq = 3; % smoothing parameter
cfg.keeptrials = 'yes';
cfg.channel = ccomb;
cfg.channelcmb = ccomb;
%cfg.tapsmofrq=ones(length(cfg.foi));
%cfg.t_ftimwin = 3./cfg.foi;
%cfg.taper = 'hanning'; % default 'dpss'

freq = ft_freqanalysis(cfg, data_preTMS);

% compute coherence
cfg = [];
cfg.method = 'coh';
cfg.channelcmb = ccomb;
fd = ft_connectivityanalysis(cfg, freq);

% display the coherence
cfg = [];
cfg.parameter = 'cohspctrm';
cfg.xlim = [8 40];
cfg.refchannel = 'EMGright';
cfg.layout = lay;
cfg.showlabels = 'yes';
figure, ft_multiplotER(cfg,fd), close figure 1 

%% 3.1.2 Define EOI & FOI
% determine freq of max value of CMC 
% in sensorimotor electrodes (e17, e18)

% find electrode of interest (eoi)
[cmcpk16 foi16] = max(fd.cohspctrm(16,5:28));
[cmcpk17 foi17] = max(fd.cohspctrm(17,5:28));
[cmcpk18 foi18] = max(fd.cohspctrm(18,5:28));

% EXPL.: 
% fd.cohspctrm contains cmc values for frq limits 8-40
% whereas size(fd.cohspctrm,2) -> 33.
% e.g.: fd.cohspctrm(1) -> value for 8Hz; fd.cohspctrm(33) -> value for 40Hz
% max cmc value shell be calculated in beta for 12-35 Hz.
% e.g.: fd.cohspctrm(1+4) -> cmc at 8+4= 12Hz
%       fd.cohspctrm(33-5) -> cmc at 40-5= 35Hz
% therefore max(fd.cohspctrm(e,5:28)) calculates the max of cmc at
% electrode e between frequencies 12-35. 
% fd.cohspctrm(e, 5:28) has different indices again: 
% size(fd.cohspctrm(e, 5:28),2) -> 24. 1st element corresponds to 12Hz
% and fd.cohspctrm(e, 5); 24th element corresponds to 35Hz and fd.cohspctrm(28) 
% due to index transformation, +11 has to be added to foi[16 17 18] to
% output the correct peak freq. of cmc.

% TABLE
% cmc frqs.                 8 9 ... 12 13 14 ... 33 34 35 ...  39 40
% fd.cohspctrm indc.        1 2 ...  5  6  7 ... 26 27 28 ...  32 33   (+7)
% fd.cohspctrm(5:28) indc.           1  2  3 ... 22 23 24              (+4)


if (cmcpk17 > cmcpk18) && (cmcpk17 > cmcpk16) % if cmc pk lies in e17
    eoi = 17;
    foi = foi17+11;
elseif (cmcpk18 > cmcpk16) && (cmcpk18 > cmcpk17) % if cmc pk lies in e18
    eoi = 18;
    foi = foi18+11;
elseif (cmcpk16 > cmcpk17) && (cmcpk16 > cmcpk18) % if cmc pk lies in e18
    eoi = 16;
    foi = foi16+11;
       
end

% plot coherence in selected channel
CLall = 1-(1-95/100)^(1/(numel(freq.powspctrm(:,1,1))-1));
xval = 8:40; 
CLall(1,1:length(xval)) = CLall; 

% assign first 7 values to 0, values from 8-40 will be cmc values from
% fd.cohspctrm -> necessary for plotting 
transVec = zeros(1, 40);
for ii = 0:(length(fd.cohspctrm(eoi,:))-1)
    transVec(ii+8) = fd.cohspctrm(eoi,ii+1);
end

figure, plot(transVec)
axis tight
xlim([8 40])
hold on
plot(xval, CLall, '-- k') 
title('Corticomuscular Coherence over trials'), xlabel('Frequency')

fprintf('CMC peak at %d Hz, strongest in electrode %d \n', foi, eoi)


% topoplot
wantTopo=0;
if wantTopo == 1
    cfg = [];
    cfg.parameter = 'cohspctrm';
    cfg.xlim = [1 50];
    cfg.zlim = [0 0.08];
    cfg.refchannel = 'EMGright';
    cfg.layout = lay;
    figure; ft_topoplotER(cfg, fd), close figure 3
end

%% 3.2 Get cross- & power & cmc spectra over trials in conditions
load(fullfile(sb_resultDir, [num2str(IDs{iSub}) '_MEP_trlNum.mat']))
trlVec = [trlsHigh; trlsLow];
store_freq = cell(2,2);

for ii = 1:length(trlVec(:,1))
    
    % cross- & power-spectra
    cfg = [];
    cfg.output = 'powandcsd';
    cfg.method = 'mtmfft'; 
    cfg.foi = [8:1:40];
    %cfg.t_ftimwin = 3./cfg.foi;
    %cfg.taper = 'hanning'; % default 'dpss'
    cfg.tapsmofrq = 3; % smoothing parameter
    cfg.keeptrials = 'yes';
    cfg.channel = ccomb;
    cfg.channelcmb = ccomb;
    cfg.trials = trlVec(ii,:);
    freq1 = ft_freqanalysis(cfg, data_preTMS);
               
    % coherence
    cfg = [];
    cfg.method = 'coh';
    cfg.channelcmb = ccomb;
    fd1 = ft_connectivityanalysis(cfg, freq1);
    
    % store output structs in memory cell
    % one row stores data from one condition
    % first row contains CMC for high MEPs, second row from low MEPs
    store_freq{ii,1} = freq1;
    store_freq{ii,2} = fd1;
      
end

% transform output memory cell to output structure for overview
output_val.frqH=store_freq{1,1}; output_val.fdH=store_freq{1,2};
output_val.fdL=store_freq{2,2};  output_val.frqL=store_freq{2,1}; 

% calculate CMC peak value for high and low MEP condition
% implemented as such, maximum could lay at a different frequency
pkH = output_val.fdH.cohspctrm(eoi,foi-7); % trnslate back to 
pkL = output_val.fdL.cohspctrm(eoi,foi-7);

% save coh values for both conditions in vector 
% [subject ID x subject Num x cmcFrq x EOI x CMC_h x CMC_l x Cohspctrm]
if exist('CMC_values_in_both_conditions.mat')
    load(fullfile(sb_groupDir, 'CMC_values_in_both_conditions.mat'), 'Coh_values')
end
Coh_values{iSub} = struct('subNum', iSub, 'subID', IDs{iSub}, 'cmcFrq', foi, 'EOI', eoi, 'CohHighMEP', pkH, 'CohLowMEP', pkL, 'Cohspctrm', output_val.fdH, ...
    'CohspctrmHigh', output_val.fdH.cohspctrm(eoi,:), 'CohspctrmLow', output_val.fdL.cohspctrm(eoi,:));
save(fullfile(sb_groupDir, 'CMC_values_in_both_conditions.mat'), 'Coh_values')

% plot
CL = 1-(1-95/100)^(1/(numel(freq1.powspctrm(:,1,1))-1)); 
xval = 8:40; 
CLval(1,1:length(xval)) = CL;

% assign first 7 values to 0, values from 8-40 will be cmc values from
% output_val.fdH.cohspctrm/output_val.fdL.cohspctrm  -> necessary for plotting 
transVecH = zeros(1, 40);
transVecL = zeros(1, 40);
  
for ii = 0:(length(output_val.fdH.cohspctrm(eoi,:))-1)
    transVecH(ii+8) = output_val.fdH.cohspctrm(eoi,ii+1);
end

for ii = 0:(length(output_val.fdL.cohspctrm(eoi,:))-1)
    transVecL(ii+8) = output_val.fdL.cohspctrm(eoi,ii+1);
end

% plot 
figure, plot(transVecH)
axis tight
hold on
plot(transVecL)
plot(xval, CLval, '-- k')
xlim([8 40])
legend('highMEPs', 'lowMEPs', '95% CL')
title('Corticomuscular Coherence'), xlabel('Frequency')

fprintf('CMC peak at %d Hz, strongest in electrode %d \n', foi, eoi)
fprintf('CMC peak in high MEP is %f and %f in low MEP \n', pkH, pkL)
%% FT plots
wantFTplots = 0;
if wantFTplots == 1
    cfg = [];
    cfg.parameter = 'cohspctrm';
    cfg.xlim = [1 60];
    cfg.refchannel = 'EMGright';
    cfg.layout = lay;
    cfg.showlabels = 'yes';
    cfg.channel = lay.label(eoi); % choose channel
    ft_singleplotER(cfg, output_val.fdH)

    cfg = [];
    cfg.parameter = 'cohspctrm';
    cfg.xlim = [1 60];
    cfg.refchannel = 'EMGright';
    cfg.layout = lay;
    cfg.showlabels = 'yes';
    cfg.channel = lay.label(eoi); % choose channel
    ft_singleplotER(cfg, output_val.fdL)
end

%% convolution based analysis of cmc
conv = 0;
if conv ~= 0
    cfg = [];
    cfg.output = 'fourier';
    cfg.method = 'mtmconvol';
    cfg.foi = [5:40];
    cfg.t_ftimwin = 3./cfg.foi;
    cfg.toi = (cfg.t_ftimwin/2).*-1;
    cfg.tapsmofrq = 2;
    cfg.taper = 'hanning';
    cfg.keeptrials = 'yes';
    cfg.channel = ccomb;
    cfg.channelcmb = ccomb;

    freq = ft_freqanalysis(cfg, data_preTMS);

    % compute coherence
    cfg = [];
    cfg.method = 'coh';
    cfg.channelcmb = ccomb;
    fd = ft_connectivityanalysis(cfg, freq);

    % display the coherence
    cfg = [];
    cfg.parameter = 'cohspctrm';
    cfg.xlim = [5 40];
    cfg.refchannel = 'EMGright';
    cfg.layout = lay;
    cfg.showlabels = 'yes';
    figure, ft_multiplotER(cfg,fd), close figure 1 
end
