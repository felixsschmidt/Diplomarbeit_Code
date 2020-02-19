% Preprocessing of 1s pre and post intervention EEG and EMG data
% 1. Trial definition
% 2. EEG preprocessing
% 3. EEG artifact rejec
% 4. EMG preprocessing
% 5. Merge data

%% Set up
clear, clc, close all
cd '/mnt/projects/CMCloop/data/random_stimulation_experiment/FT_Analysis'
setup_path; 
addpath(genpath(fullfile(projectdir, 'data'))); 
% Folder information
iSub = 11;                           % adjust
sessionNum = 2;
sb_rawDir = fullfile(baseDir,indiDir,IDs{iSub},rawDir);
sb_preprocDir = fullfile(baseDir,indiDir,IDs{iSub},preprocDir);
sb_resultDir = fullfile(baseDir,indiDir,IDs{iSub}, resultDir);
sb_groupDir = fullfile(baseDir, groupDir);

%% 1. Trial definition and -selection, Import EEG+EMG data 
% include 3s around TMS artifact
cfg = [];
if sessionNum == 1
    cfg.dataset = fullfile(sb_rawDir, ['X85256_CMC_no_stim_pre' '.eeg']); % adjust
else
    cfg.dataset = fullfile(sb_rawDir, ['X85256_CMC_no_stim_post' '.eeg']); % adjust
end
cfg.trialfun = 'ft_trialfun_general';
cfg.trialdef.eventtype = 'Stimulus';
cfg.trialdef.eventvalue = 'A';
cfg.trialdef.prestim = 0.5; 
cfg.trialdef.poststim = 0.5; 
cfg.channel = {'all'}; 
cfg = ft_definetrial(cfg);
% save cfg structure for later; incorporates cfg.trl & cfg.dataset
save(fullfile(sb_preprocDir, 'trial_def_cfg_nonStim.mat'), 'cfg');  

dataRaw = ft_preprocessing(cfg);

%figure, plot(dataRaw.time{1}, dataRaw.trial{1} (17,:))
%%  2. Preprocessing EEG
%% 2.1 FILTERING 
cfg = [];
cfg.detrend = 'yes'; % remove linear trend
cfg.demean = 'yes'; % baseline correction
cfg.baselinewindow = 'all'; 

% High-pass filter
cfg.hpfilter = 'yes'; % apply hp filter
cfg.hpfreq = 1; % hp filter freq 1 Hz
cfg.hpfilttype = 'firws'; % windowed FIR filter 

% Low-pass filter
cfg.lpfilter = 'yes'; 
cfg.lpfreq = 100; % lp filter freq 100 Hz
cfg.lpfilttype = 'firws'; 

% Notch Filter
notchFreq = 50; % notch filter freq 50 Hz line noise
cfg.bsfilter = 'yes';
cfg.bsfreq = [notchFreq-1 notchFreq+1];
cfg.bsfilttype = 'firws';

% only on eeg channels
cfg.channel = {'all' '-EMGright' '-EMGleft' }; 

dataEEG_filtered = ft_preprocessing(cfg, dataRaw);

% update layout
dataEEG_filtered.label = lay.label;

% Downsample to 500 Hz
downsampFreq = 500; 
cfg = [];
cfg.resamplefs = downsampFreq;
dataEEG_filtered_rs = ft_resampledata(cfg, dataEEG_filtered);

%figure, plot(dataEEG_filtered_rs.time{1}, dataEEG_filtered_rs.trial{1} (17,:))

%% 2.2 INTERPOLATION 
% removing channels and fill them with average from neighbors
% should be done before re-referencing

% define elec structure
dataEEG_filtered_rs.elec.chanpos = lay.pos;
dataEEG_filtered_rs.elec.elecpos = lay.pos;
dataEEG_filtered_rs.elec.label = lay.label;

% define neighbours
cfg = [];
cfg.layout = lay;
cfg.method = 'triangulation';
cfg.channel = lay.label;
neighbours = ft_prepare_neighbours(cfg);
% 5 neighbour electrodes for e16, e17 & e18
neighbours(18).neighblabel = {'e6';'e15';'e17';'e29';'e30'};
neighbours(17).neighblabel = {'e31';'e32';'e18';'e7';'e6'}; 
neighbours(18).neighblabel = {'e32';'e33';'e19';'e7';'e17'};

%%% ONLY REMOVE BAD CHANNELS HERE %%%

% OBS: DONT REMOVE CHANNELS DURING VISUAL INSPECTION SUMMARY
% JUST FIND OUT WHAT CHANNELS ARE BAD AND ENTER IT MANUALLY IN NEXT LINE
% otherwise, channel order will be disrupted
cfg = [];
cfg.method = 'summary';
cfg.layout = eeglayout;
dataEEG_filtered_rs = ft_rejectvisual(cfg, dataEEG_filtered_rs);

% manually add the names of removed channels, e.g. {'e36'; 'e51'}
artif.badchannel  = input('write badchannels: ');

% interpolate
cfg = [];
cfg.elec = dataEEG_filtered_rs.elec;
cfg.badchannel    = artif.badchannel;
cfg.method         = 'weighted';
cfg.neighbours     = neighbours;
dataEEG_interpol = ft_channelrepair(cfg, dataEEG_filtered_rs);

%% 2.3 HJORTH MONTAGE REF
cfg = [];
cfg.method = 'hjorth';
cfg.elec = dataEEG_interpol.elec;
cfg.trials = 'all';
cfg.neighbours = neighbours;
dataEEG_reref = ft_scalpcurrentdensity(cfg, dataEEG_interpol);
 
%% 3. Visual artifact rejection EEG

%%% ONLY REMOVE BAD TRIALS HERE %%%
cfg = [];
cfg.method = 'summary';
cfg.layout = eeglayout;

if sessionNum == 1 
    dataEEG_vrejec_pre = ft_rejectvisual(cfg, dataEEG_reref);
    save(fullfile(sb_preprocDir , [num2str(IDs{iSub}) '_preproc_EEG_dataPre.mat']), 'dataEEG_vrejec_pre');
else
    dataEEG_vrejec_post = ft_rejectvisual(cfg, dataEEG_reref);
    save(fullfile(sb_preprocDir , [num2str(IDs{iSub}) '_preproc_EEG_dataPost.mat']), 'dataEEG_vrejec_post');
end

%% 4 EMG preprocessing + artifact rejection (rectified)
cfg = [];
cfg.demean = 'yes';
cfg.dftfilter = 'yes'; % notch filter
cfg.hpfilter = 'yes';
cfg.hpfreq = 10;
cfg.hpfilttype = 'firws';
cfg.rectify = 'yes';
cfg.channel = {'EMGright' 'EMGleft'}; % only EMG channels
cfg.reref = 'yes';                    % subtract reference (EMGleft) from active electrode (EMGright)
cfg.refchannel = {'EMGleft'};

dataEMG_filtered = ft_preprocessing(cfg, dataRaw);

% Downsample to 500 Hz
cfg = [];
cfg.resamplefs = downsampFreq;
dataEMG_filtered_rs = ft_resampledata(cfg, dataEMG_filtered);

% remove rejected EEG trials from EMG structure
cfg=[];
load(fullfile(sb_preprocDir, 'trial_def_cfg_nonStim.mat'), 'cfg'); % incorporates cfg.trl & cfg.dataset
if sessionNum ==1 
    cfg.artfctdef = dataEEG_vrejec_pre.cfg.artfctdef; 
else
    cfg.artfctdef = dataEEG_vrejec_post.cfg.artfctdef;
end
dataEMG_rejec = ft_rejectartifact(cfg, dataEMG_filtered_rs);

% remove EMGleft ref chan
cfg = [];
cfg.method = 'summary';
cfg.layout = eeglayout;

if sessionNum == 1 
    dataEMG_rejec_pre = ft_rejectvisual(cfg, dataEMG_rejec);
    save(fullfile(sb_preprocDir , [num2str(IDs{iSub}) '_preproc_EMG_data_pre.mat']), 'dataEMG_rejec_pre');
else
    dataEMG_rejec_post = ft_rejectvisual(cfg, dataEMG_rejec);
    save(fullfile(sb_preprocDir , [num2str(IDs{iSub}) '_preproc_EMG_data_post.mat']), 'dataEMG_rejec_post');
end

%% 5. Merge data
if sessionNum == 1
    cfg = [];
    cfg.keepsampleinfo = 'yes';
    dataEMG_rejec_pre.time = dataEEG_vrejec_pre.time; % adjust time
    data_preStim = ft_appenddata(cfg, dataEEG_vrejec_pre, dataEMG_rejec_pre);
    data_preStim=ft_struct2single(data_preStim); 
    save(fullfile(sb_preprocDir, [num2str(IDs{iSub}), '_preStim_EEG_EMG.mat']), 'data_preStim');
end


if sessionNum == 2
    cfg = [];
    cfg.keepsampleinfo = 'yes';
    dataEMG_rejec_post.time = dataEEG_vrejec_post.time; % adjust time
    data_postStim = ft_appenddata(cfg, dataEEG_vrejec_post, dataEMG_rejec_post);
    data_postStim=ft_struct2single(data_postStim);
    save(fullfile(sb_preprocDir, [num2str(IDs{iSub}), '_postStim_EEG_EMG.mat']), 'data_postStim');
end
