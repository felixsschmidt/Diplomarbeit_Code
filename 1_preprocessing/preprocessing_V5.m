%% Preprocesing Pipeline
% 1. Trial definition
% 2. TMS artifact rejection
% 3. Preproc EEG data
% 4. Artifact rejection EEG
% 5.1 Preproc and artifact rejection preTMS EMG data 
% 5.2 Preproc and artifact rejection postTMS EMG data
% 6. Append and sasve datasets 
%% Set up
clear, clc, close all
cd '/mnt/projects/CMCloop/data/random_stimulation_experiment/FT_Analysis/'
setup_path; 
addpath(genpath(fullfile(projectdir, 'data'))); 

% Folder information
iSub = 11;                                                          % define
sessionNum = 2;                                                    % define
sb_rawDir = fullfile(baseDir,indiDir,IDs{iSub},rawDir);
sb_preprocDir = fullfile(baseDir,indiDir,IDs{iSub},preprocDir);
sb_preStimDir = fullfile(baseDir,indiDir,IDs{iSub},preprocDir, preStimDir);
sb_postStimDir = fullfile(baseDir,indiDir,IDs{iSub},preprocDir, postStimDir);

%% 1. Trial definition and -selection, Import EEG+EMG data 
% include 3s around TMS artifact
cfg = [];
cfg.dataset = fullfile(sb_rawDir, ['X85256_CMC_random_session2' '.eeg']); % adjust 
cfg.trialfun = 'ft_trialfun_general';
cfg.trialdef.eventtype = 'Stimulus';
cfg.trialdef.eventvalue = 'A';
cfg.trialdef.prestim = 3; 
cfg.trialdef.poststim = 3; 
cfg.channel = {'all'}; 
cfg = ft_definetrial(cfg);

% save cfg structure for later; incorporates cfg.trl & cfg.dataset
save(fullfile(sb_preprocDir, 'trial_def_cfg.mat'), 'cfg');  

dataRaw = ft_preprocessing(cfg); % upload unfiltered raw data

%% 2. Detect and remove TMS artifacts
% it's done now on all the data

% specify structure for TMS artifact removal
cfg.trials = 1:length(dataRaw.trial); 
cfg.plottrial = 'no';
cfg.ampsd = 1;
cfg.refwin = [200 20]; 
cfg.interp = 'std'; % fill with noise approach
cfg.method = 'single';
cfg.pretime = 20; 
cfg.posttime = 15; % 20ms around TMS are interpolated with noise

data_TMSr = cleanTMS_V2(cfg, dataRaw);

% plot to check after artifact rejection
dispExCh = 0;   % plot example channel to check, if TMS artifact was removed?

if dispExCh == 1
    figure 
    plot(data_TMSr.time{1}, data_TMSr.trial{1}(17,:));
    axis tight;
    legend(data_TMSr.label(17)); 
end

%% 3. Preprocessing EEG
%%% 3.1 FILTERING %%%
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

dataEEG_filtered = ft_preprocessing(cfg, data_TMSr);

% update layout
dataEEG_filtered.label = lay.label;

% Downsample to 500 Hz
downsampFreq = 500; 
cfg = [];
cfg.resamplefs = downsampFreq;
dataEEG_filtered_rs = ft_resampledata(cfg, dataEEG_filtered);

%%% 3.2 INTERPOLATION %%%
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
 
%%% 3.3 RE-REFERENCE %%%
refAppr = 1;
if refAppr == 1
    
    %%% 3.3.1 HJORTH MONTAGE %%%
    cfg = [];
    cfg.method = 'hjorth';
    cfg.elec = dataEEG_interpol.elec;
    cfg.trials = 'all';
    cfg.neighbours = neighbours;
    dataEEG_reref = ft_scalpcurrentdensity(cfg, dataEEG_interpol);
    
else
    
    %%% 3.3.2 Common av reference %%%
    cfg = [];
    cfg.reref = 'yes';
    cfg.refchannel = 'all';
    dataEEG_reref = ft_preprocessing(cfg, dataEEG_interpol);
    
end

disp = 0;
if disp == 1
    cfg = [];
    cfg.viewmode = 'vertical';
    ft_databrowser(cfg, dataEEG_reref);
end

dispExCh = 0;
if dispExCh == 1
    figure 
    plot(dataEEG_filtered_rs.time{1}, dataEEG_filtered_rs.trial{1}(17,:));
    axis tight;
    legend(dataEEG_filtered_rs.label(17)); 
end

%% 4. Visual artifact rejection EEG

%%% ONLY REMOVE BAD TRIALS HERE %%%
cfg = [];
cfg.method = 'summary';
cfg.layout = eeglayout;


if sessionNum == 1 
    dataEEG_vrejec1 = ft_rejectvisual(cfg, dataEEG_reref);
    save(fullfile(sb_preStimDir, [num2str(IDs{iSub}) '_preproc_EEG_data1.mat']), 'dataEEG_vrejec1');
else
    dataEEG_vrejec2 = ft_rejectvisual(cfg, dataEEG_reref);
    save(fullfile(sb_preStimDir, [num2str(IDs{iSub}) '_preproc_EEG_data2.mat']), 'dataEEG_vrejec2');
end

%% 5.1 preTMS EMG preprocessing + artifact rejection (rectified)
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

dataEMG_filtered = ft_preprocessing(cfg, data_TMSr);

% Downsample to 500 Hz
cfg = [];
cfg.resamplefs = downsampFreq;
dataEMG_filtered_rs = ft_resampledata(cfg, dataEMG_filtered);

% remove rejected EEG trials from EMG structure
cfg=[];
load(fullfile(sb_preprocDir, 'trial_def_cfg.mat'), 'cfg'); % incorporates cfg.trl & cfg.dataset
if sessionNum ==1 
    cfg.artfctdef = dataEEG_vrejec1.cfg.artfctdef; 
else
    cfg.artfctdef = dataEEG_vrejec2.cfg.artfctdef;
end
dataEMG_rejec = ft_rejectartifact(cfg, dataEMG_filtered_rs);

% remove EMGleft ref chan
cfg = [];
cfg.method = 'summary';
cfg.layout = eeglayout;

if sessionNum == 1 
    dataEMG_rejec1 = ft_rejectvisual(cfg, dataEMG_rejec);
    save(fullfile(sb_preStimDir, [num2str(IDs{iSub}) '_preproc_EMG_data1.mat']), 'dataEMG_rejec1');
else
    dataEMG_rejec2 = ft_rejectvisual(cfg, dataEMG_rejec);
    save(fullfile(sb_preStimDir, [num2str(IDs{iSub}) '_preproc_EMG_data2.mat']), 'dataEMG_rejec2');
end

%% 5.2 postTMS EMG preprocessing + artifact rejection (unrectified)
cfg = [];
cfg.demean = 'yes';
cfg.dftfilter = 'yes'; % notch filter
cfg.hpfilter = 'yes';
cfg.hpfreq = 10;
cfg.hpfilttype = 'firws';
cfg.rectify = 'no';
cfg.channel = {'EMGright' 'EMGleft'}; % only EMG channels
cfg.reref = 'yes';                    % subtract reference (EMGleft) from active electrode (EMGright)
cfg.refchannel = {'EMGleft'};

UnRec_dataEMG_filtered = ft_preprocessing(cfg, data_TMSr);

% Downsample to 500 Hz
cfg = [];
cfg.resamplefs = downsampFreq;
UnRec_dataEMG_filtered_rs = ft_resampledata(cfg, UnRec_dataEMG_filtered);

% remove rejected EEG trials from EMG structure
cfg = [];
load(fullfile(sb_preprocDir, 'trial_def_cfg.mat'), 'cfg'); % incorporates cfg.trl & cfg.dataset
if sessionNum ==1 
    cfg.artfctdef = dataEEG_vrejec1.cfg.artfctdef; 
else
    cfg.artfctdef = dataEEG_vrejec2.cfg.artfctdef;
end
UnRec_dataEMG_rejec = ft_rejectartifact(cfg, UnRec_dataEMG_filtered_rs);

% remove EMGleft ref chan
cfg = [];
cfg.method = 'summary';
cfg.layout = eeglayout;

if sessionNum == 1 
    UnRec_dataEMG_rejec1 = ft_rejectvisual(cfg, UnRec_dataEMG_rejec);
    save(fullfile(sb_postStimDir, [num2str(IDs{iSub}) '_UnRec_preproc_EMG_data1.mat']), 'UnRec_dataEMG_rejec1');
else
    UnRec_dataEMG_rejec2 = ft_rejectvisual(cfg, UnRec_dataEMG_rejec);
    save(fullfile(sb_postStimDir,[num2str(IDs{iSub}) '_UnRec_preproc_EMG_data2.mat']), 'UnRec_dataEMG_rejec2');
end    

%% 6. Merge datasets
if sessionNum == 2
    merge_data
end

%%
wantPlot = 0; 
if wantPlot == 1
    cfg = [];
    cfg.viewmode = 'vertical';
    ft_databrowser(cfg, data_preproc)

    % plot to confirm
    figure % EEG data over C3
    subplot(2,1,1)
    plot(data_preproc.time{32}, data_preproc.trial{32}(17,:));
    axis tight;
    legend(data_preproc.label(17)); 

    subplot(2,1,2) % EMG data
    plot(data_preproc.time{32}, data_preproc.trial{32}(64,:));
    axis tight;
    legend(data_preproc.label(64)); 
   
end




