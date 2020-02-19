% phase analysis of 1s preTMS eeg

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
sb_figDir = fullfile(baseDir, figDir);

%% 1. Define 1s preTMS period [-1020 -.02] ms
load(fullfile(sb_preStimDir, [num2str(IDs{iSub}) '_preTMS_EEG_EMG.mat']));

cfg = [];
cfg.toilim = [-1.02 -0.02];
data_preTMS = ft_redefinetrial(cfg, data_preproc);

%% 2. Phase analysis for both conditions
%% 2.1 Convolution based, time resolved estimate of phase

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% frequencies (foi):       8-40 Hz
% time (toi):              1 s
% time window (t_ftimwin): freq. specific: 3 cycles/freq. 
% stepsize & time (toi):   window slides from -1.02 till -0.02 ms preTMS
%                          in steps of 50 ms
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

load(fullfile(sb_resultDir, [num2str(IDs{iSub}) '_MEP_trlNum.mat']))
trlVec = [trlsHigh; trlsLow];

for ii = 1:length(trlVec(:,2))
    
    cfg = [];
    cfg.channel = {'all' '-EMGright'};
    cfg.trials = trlVec(ii,:); 
    cfg.method = 'mtmconvol';
    cfg.output = 'fourier';
    cfg.taper = 'hanning';
    cfg.foi = 8:40;                              
    cfg.t_ftimwin = 3./cfg.foi; 
    cfg.toi = -1.02:0.05:-0.02; 
    if ii == 1                            
        Hfourier = ft_freqanalysis(cfg, data_preproc);
    else
        Lfourier = ft_freqanalysis(cfg, data_preproc);  
    end
    
end
 
%% 2.2 Determine mean phases across trials 
%      for each respective freq in high vs low MEP trials

% fourierspctrm: [56×63×33×21 double]
%%%%%% (1) IDX VEC  %%%%%%
% idx vector with freq specific idx of tw closest to stimulus onset (SO)
% this vector can be used across trials for the same freq to indicate
% tw of interest for that freq

% vec can be used as 4th argument for
% Hphase.fourierspctrm(trl,17,frq,idx_v(frq))

% suitable for high as well as low MEP trials

frqOfint = 8:40;
twLength = zeros(2,length(frqOfint)); twLength(2,:) = 8:40;
ind = zeros(1,length(Hfourier.time));
idx_v = zeros(2,length(Hfourier.freq));

for ii = 1:length(Hfourier.freq)
    
    twLength(1, ii) = 3/frqOfint(ii); % tw covering 3 cycles
    % output indices from time variable,
    % where tw center + tw/2 dont exceed stimulus onset  
    ind = (Hfourier.time + (twLength(1,ii)/2)) <= (-0.019999);
    % pick center point from tw, where it is closest to stimulus onset
    [val idx] = min(abs(Hfourier.time(ind)));
    idx_v(1,ii) = idx;          % save indx for respective freq
    idx_v(2,ii) = frqOfint(ii); % save corresponding freq
end

%%%%%% (2) %%%%%%
% fourierspctrm: [56×63×33×21 double]
% get phase values per freq over trials

% (2.1) HIGH MEP trials
numOftrls = length(Hfourier.fourierspctrm(:,1,1,1));
numOffreq = length(Hfourier.freq);
fourierCoefH = zeros(numOffreq, numOftrls); 
phase_valuesH = zeros(numOffreq, numOftrls);% [rows:freq x columns:trl]

for freq = 1:numOffreq % for every freq 8-40 Hz
    
    for trl = 1:numOftrls % for every trl
        fourierCoefH(freq, trl) = Hfourier.fourierspctrm(trl,17,freq,...
            idx_v(1,freq));
        phase_valuesH(freq, trl) = angle(Hfourier.fourierspctrm(trl,17,freq,...
            idx_v(1,freq)));
        
        
    end
    
end

% (2.2) LOW MEP trials
numOftrls = length(Lfourier.fourierspctrm(:,1,1,1));
numOffreq = length(Lfourier.freq);
fourierCoefL = zeros(numOffreq, numOftrls); 
phase_valuesL = zeros(numOffreq, numOftrls);% [rows:freq x columns:trl]

for freq = 1:numOffreq % for every freq 8-40 Hz
    
    for trl = 1:numOftrls % for every trl
        fourierCoefL(freq, trl) = Lfourier.fourierspctrm(trl,17,freq,...
            idx_v(1,freq));
        phase_valuesL(freq, trl) = angle(Lfourier.fourierspctrm(trl,17,freq,...
            idx_v(1,freq)));
        
        
    end
    
end


%% plot data
% mean out of phases per freq: every freq, 1 mean phase per subj
% mean phase val for each frq across high MEP trials 
MfrqThetaH = zeros(numOffreq,1);
figure(1), subplot(1,2,1)
for freq = 1:numOffreq
    
    MfrqThetaH(freq,1) = circ_mean(phase_valuesH(freq,:));
    p = polarplot( [0 MfrqThetaH(freq)], [0 abs(fourierCoefH(freq))],...
        'LineWidth',1.5); 
    hold on
end
title('Mean phases across frequencies for one subject in high MEP trials')

subplot(1,2,2)
% mean phase val for each frq across low MEP trials 
MfrqThetaL = zeros(numOffreq,1);
for freq = 1:numOffreq
    
    MfrqThetaL(freq,1) = circ_mean(phase_valuesL(freq,:));
    p = polarplot( [0 MfrqThetaL(freq)], [0 abs(fourierCoefL(freq))],...
        'LineWidth',1.5); 
    hold on
end
title('Mean phases across frequencies for one subject in low MEP trials')

%% save data

% save phase values for both conditions in vector
% [subject x ID x Mphase_h x Mphase_l x fourierCoef_h x fourierCoef_l]
if exist('phase_values_in_both_conditions.mat')
    load(fullfile(sb_groupDir, 'phase_values_in_both_conditions.mat'), 'phase_values')
end
phase_values{iSub} = struct('subNum', iSub, 'subID', IDs{iSub}, ...
    'Mean_hMEP_phase_val', MfrqThetaH, 'Mean_lMEP_phase_val', MfrqThetaL,...
    'fourierCoef_H', fourierCoefH, 'fourierCoefL', fourierCoefL);
save(fullfile(sb_groupDir, 'phase_values_in_both_conditions.mat'), 'phase_values' )


%%
%abs(z).*exp(i*theta)
%figure, polarplot( [0 theta(1)], [0 z(1)], 'r' ), hold on
%polarplot( [0 theta(2)], [0 z(2)], 'b' ), hold on
%polarplot( [0 theta(3)], [0 z(3)], 'k' )