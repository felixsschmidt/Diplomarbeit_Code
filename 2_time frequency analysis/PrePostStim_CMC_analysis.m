% CMC analysis of pre and post intervention trials
%% Set up
clear, clc, close all
cd '/mnt/projects/CMCloop/data/random_stimulation_experiment/FT_Analysis'
setup_path; 
addpath(genpath(fullfile(projectdir, 'data'))); 
% Folder information
iSub = 11;                           % adjust
eoi = 16; foi = 20;                 % adjust
sb_preprocDir = fullfile(baseDir,indiDir,IDs{iSub},preprocDir);
sb_resultDir = fullfile(baseDir,indiDir,IDs{iSub}, resultDir);
sb_groupDir = fullfile(baseDir, groupDir);

%% 1. Load data and define Nx2 cell for channelcmb
load(fullfile(sb_preprocDir, [num2str(IDs{iSub}), '_preStim_EEG_EMG.mat']), 'data_preStim');
load(fullfile(sb_preprocDir, [num2str(IDs{iSub}), '_postStim_EEG_EMG.mat']), 'data_postStim');

% combine all EEG channel with EMGright
ccomb = cell(length(data_preStim.label)-1,2);
for ii=1:length(data_preStim.label)-1
    ccomb{ii,1} = char(lay.label(ii));
    ccomb{ii,2} = 'EMGright';       
end

%% 2. Corticomuscular Coherence Analysis pre and post intervention
%% 2.1 CMC Analysis in preStim condition
store_freq = cell(2,2);

for ii = 1:2   
    % cross- & power-spectra
    cfg = [];
    cfg.output = 'powandcsd';
    cfg.method = 'mtmfft'; 
    cfg.foi = [8:1:40];
    cfg.tapsmofrq = 3; % smoothing parameter
    cfg.keeptrials = 'yes';
    cfg.channel = ccomb;
    cfg.channelcmb = ccomb;
    cfg.trials = 'all';
    if ii == 1       
        freq = ft_freqanalysis(cfg, data_preStim);        
    else
        freq = ft_freqanalysis(cfg, data_postStim);
    end     
                  
    % coherence
    cfg = [];
    cfg.method = 'coh';
    cfg.channelcmb = ccomb;
    if ii == 1       
        fd = ft_connectivityanalysis(cfg, freq);
    else
        fd = ft_connectivityanalysis(cfg, freq);
    end    
    
    % store output structs in memory cell
    % one row stores data from one condition
    % first row contains CMC for high MEPs, second row from low MEPs
    store_freq{ii,1} = freq;
    store_freq{ii,2} = fd;
      
end

% transform output memory cell to output structure for overview
output_val.frqPre=store_freq{1,1}; output_val.fdPre=store_freq{1,2};
output_val.fdPost=store_freq{2,2};  output_val.frqPost=store_freq{2,1}; 

% calculate CMC values for each condition in freq that has been used for
% comparison between high and low trials
pkPre = output_val.fdPre.cohspctrm(eoi,foi-7); 
pkPost = output_val.fdPost.cohspctrm(eoi,foi-7);

% save
% [subject ID x subject Num x cmcFrq x EOI x CMC_pre x CMC_post]
if exist('CMC_values_in_pre_and_post_condition.mat')
    load(fullfile(sb_groupDir, 'CMC_values_in_pre_and_post_condition.mat'), 'Coh_values_nonStim')
end
Coh_values_nonStim{iSub} = struct('subNum', iSub, 'subID', IDs{iSub}, 'cmcFrq', foi, 'EOI', eoi, 'CohPreStim', pkPre, 'CohPostStim', pkPost,...
    'CohspctrmPre', output_val.fdPre.cohspctrm(eoi,:), 'CohspctrmPost', output_val.fdPost.cohspctrm(eoi,:));
save(fullfile(sb_groupDir, 'CMC_values_in_pre_and_post_condition.mat'), 'Coh_values_nonStim')


% plot
freq = [8:40];
figure(1), plot(freq,output_val.fdPre.cohspctrm(eoi,:)),
hold on, plot(freq,output_val.fdPost.cohspctrm(eoi,:)), axis tight
legend('PreStim', 'PostStim')
title('Corticomuscular Coherence'), xlabel('Frequency')

fprintf('CMC peak at %d Hz, strongest in electrode %d \n', foi, eoi)
fprintf('CMC peak in preStim is %f and %f in postStim \n', pkPre, pkPost)
