% power analysis of 1s preTMS eeg

%% Set up
clear, clc, close all
cd '/mnt/projects/CMCloop/data/random_stimulation_experiment/FT_Analysis'
setup_path; 
addpath(genpath(fullfile(projectdir, 'data'))); 
% Folder information
iSub = 1;                           % adjust
sb_preStimDir = fullfile(baseDir,indiDir,IDs{iSub},preprocDir, preStimDir);
sb_resultDir = fullfile(baseDir,indiDir,IDs{iSub}, resultDir);
sb_groupDir = fullfile(baseDir, groupDir);
sb_figDir = fullfile(baseDir, figDir);

%% 1. Define 1s preTMS period [-1020 -.02] ms
load(fullfile(sb_preStimDir, [num2str(IDs{iSub}) '_preTMS_EEG_EMG.mat']));

cfg = [];
cfg.toilim = [-1.02 -0.02];
%cfg.toilim = [-1 1];
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

%% 2.  Power analysis for both conditions
%% 2.1. Convolution based, time resolved estimate of power 
load(fullfile(sb_resultDir, [num2str(IDs{iSub}) '_MEP_trlNum.mat']))
trlVec = [trlsHigh; trlsLow];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% frequencies (foi):       8-40 Hz
% time (toi):              1 s
% time window (t_ftimwin): freq. specific: 3 cycles/freq. 
% stepsize & time (toi):   window slides from -1.02 till -0.02 ms preTMS
%                          in steps of 50 ms
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for ii = 1:length(trlVec(:,2))
    
    cfg = [];
    cfg.channel = {'all' '-EMGright'};
    cfg.trials = trlVec(ii,:); 
    cfg.method = 'mtmconvol';
    cfg.output = 'pow';
    cfg.taper = 'hanning';
    cfg.foi = 8:40;                              
    cfg.t_ftimwin = 3./cfg.foi; 
    cfg.toi = -1.02:0.05:-0.02; 
    if ii == 1                            
        Hpow = ft_freqanalysis(cfg, data_preproc);
    else
        Lpow = ft_freqanalysis(cfg, data_preproc);  
    end
    
end
 
save(fullfile(sb_groupDir, 'Topo_power_Val.mat'), 'Hpow', 'Lpow')

%% 2.2 extract time-freq window closest to stim-onset for every freq
% dimension 'time' describes where time windows are centered

%%%%% (1) For high MEP trials for CP3 (16), C3 (17) & FC3 (18)%%%%%

frqOfint = 8:40; % vector with analyzed frequencies
twLength = zeros(1,length(frqOfint));
ind = zeros(1,length(Hpow.time));
power_valuesH = zeros(2,length(Hpow.freq),3);

for ii = 1:length(Hpow.freq)
    
    for jj = 16:18
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

% average power values for sensorimotor electrodes
% take mean from power values from each freq from all different electrodes
Mpow_valuesH = zeros(2,length(power_valuesH));
for ii = 1:length(power_valuesH)
    
    % first row: mean
    Mpow_valuesH(1,ii) = mean([power_valuesH(1,ii,16)...
        power_valuesH(1,ii,17) power_valuesH(1,ii,18)]);         
    % second row: frequencies
    Mpow_valuesH(2,:) = power_valuesH(2,:,18);   
    
end

%%%%% (2) For low MEP trials for CP3 (16), C3 (17) & FC3 (18) %%%%%

frqOfint = 8:40; % vector with analyzed frequencies
twLength = zeros(1,length(frqOfint));
ind = zeros(1,length(Hpow.time));
power_valuesL = zeros(2,length(Hpow.freq),3);

for ii = 1:length(Lpow.freq)
    
    for jj = 16:18
    
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

% average power values for sensorimotor electrodes
% take mean from power values from each freq from all different electrodes
Mpow_valuesL = zeros(2,length(power_valuesL));
for ii = 1:length(power_valuesL)
    
    % first row: mean
    Mpow_valuesL(1,ii) = mean([power_valuesL(1,ii,16)...
        power_valuesL(1,ii,17) power_valuesL(1,ii,18)]); 
    
    % third row: frequencies
    Mpow_valuesL(2,:) = power_valuesL(2,:,18);   
    
end

% save power values for both conditions in vector
% [subject x ID x power_h x power_l]
if exist('pow_values_in_both_conditions.mat')
    load(fullfile(sb_groupDir, 'pow_values_in_both_conditions.mat'), 'power_values')
end
power_values{iSub} = struct('subNum', iSub, 'subID', IDs{iSub}, ...
    'highMEP_pow_val', Mpow_valuesH, 'lowMEP_pow_val',Mpow_valuesL);
save(fullfile(sb_groupDir, 'pow_values_in_both_conditions.mat'), 'power_values' )

%%
%% plot
wantPlot = 1;

if wantPlot == 1
    
    figure(1) 
    plot(Mpow_valuesH(2,:),  Mpow_valuesH(1,:), 'Linewidth', 1.5); 
    hold on
    plot(Mpow_valuesL(2,:),  Mpow_valuesL(1,:), 'Linewidth', 1.5) 
    xlabel('Frequency (Hz)'), ylabel('Power (\muV^2)')
    legend('high MEP trials', 'low MEP trials')
    title('Power spectrum from sensorimotor electrodes in high vs. low MEP trials')
    axis tight 
    %saveas(figure(1),fullfile(sb_figDir, 'power_plot.jpg'))
end



%%

