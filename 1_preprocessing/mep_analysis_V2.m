%% MEP analysis pipeline
% 1. Define postTMS period to be analyzed (15-60 ms)
% 2. Compute MEPs, remove outliers
% 3. Divide MEPs by median split and save trial indices for cmc analysis
%% Set up
clear, clc, close all
cd '/mnt/projects/CMCloop/data/random_stimulation_experiment/FT_Analysis'
setup_path; 
addpath(genpath(fullfile(projectdir, 'data'))); 

% Folder information
iSub = 7;                        % adjust
sb_preStimDir = fullfile(baseDir,indiDir,IDs{iSub},preprocDir, preStimDir);
sb_postStimDir = fullfile(baseDir,indiDir,IDs{iSub},preprocDir, postStimDir);
sb_resultDir = fullfile(baseDir,indiDir,IDs{iSub}, resultDir);
sb_groupDir = fullfile(baseDir, groupDir);

%% 1. Define postTMS period [0.015 0.1] ms
load(fullfile(sb_postStimDir,[num2str(IDs{iSub}) '_postTMS_EMG.mat']));

cfg = [];
cfg.toilim = [0.015 .1]; 
data_postTMS = ft_redefinetrial(cfg, UnRec_datasets_EMG);
% plot(UnRec_datasets_EMG.time{1}, UnRec_datasets_EMG.trial{1})
plot(data_postTMS.time{1}, data_postTMS.trial{6} (1,:))

%% 2. Compute MEPs, mark & exclude outliers
%  MEPs are calculated as peak-to-peak difference for each trial 
%  and stored in vector MEP. The indices of vector MEP correspond to 
%  trial numbers.

trials = 1:length(data_postTMS.trial);          % define trials
MEP = zeros(1,length(trials));                  % prelocate MEP vector
for ltr = trials
    emgval = data_postTMS.trial{ltr};           % emg values of trial ltr
    mep = abs( max(emgval)-min(emgval) );  % respective MEP of trial ltr
    MEP(ltr) = mep;                             % stores all MEPs
end

% find outliers (MEPs >|< 3*MAD (median absolute deviations))
outlier_ind = isoutlier(MEP);
nnz(outlier_ind);                    % # of non-zero elements
trlNum_outlr = find(outlier_ind);    % trial number of nnz elements
MEP(outlier_ind) = 9999;             % outliers are marked as 9999

% plot course of MEPs
MEPlot = 1;
if MEPlot == 1
    figure, plot(MEP, 'o')
    xlabel('Trial number'), ylabel('Magnitude (arbitrary)')
    title('Course of MEPs over trials')
end

% store MEPs for subject i
if exist('MEP_values.mat')
    load(fullfile(sb_groupDir, 'MEP_values.mat'), 'MEP_values')
end

MEP_values{iSub} = struct('subNum', iSub, 'subID', IDs{iSub}, 'MEP', MEP);
save(fullfile(sb_groupDir, 'MEP_values.mat'), 'MEP_values')

%% 3. Divide MEPs by median-split-method
split=median(MEP(MEP~=9999),2); % median over row of MEPs
% MEP_low & MEP_high store MEP and trial index values  as Nx2 vectors
% first column: MEPs, second column: trial number
MEP_low = zeros(length(MEP),2); MEP_high = zeros(length(MEP),2);  % prelocate 
for ii = 1:length(MEP)
    if (MEP(ii) < split) && (MEP(ii) < 9999)
        MEP_low(ii,1) = MEP(ii);      % even or odd number of trials? 
        MEP_low(ii,2) = ii;           % MEP(ii) < split or <= split ?
    elseif (MEP(ii) > split) && (MEP(ii) < 9999)        
        MEP_high(ii,1) = MEP(ii); 
        MEP_high(ii,2) = ii; 
    end
end

% cut vectors: only MEPs and trial indices left
mask = MEP_low(:,1) ~= 0; % indicies of MEP vector ~= 0
MEP_low = MEP_low(mask,:);
mask = MEP_high(:,1) ~= 0;
MEP_high = MEP_high(mask,:);
% save trl indices in 2 vectors for trialdef. in subsequent ft_freqanalysis
trlsHigh = MEP_high(:,2)'; 
trlsLow = MEP_low(:,2)';

fprintf(['outliers at trials: ' num2str(trlNum_outlr) '\n']) 

% save
% ANY MORE VALUES? MEAN MEPs? MEPs for both cond?
save(fullfile(sb_resultDir, [num2str(IDs{iSub}) '_MEP_trlNum.mat']),'trlsHigh', 'trlsLow')


% plot MEPs
check = 1;
if check == 1
    figure
    plot(MEP(trlsHigh), 'o')
    axis tight
    hold on
    plot(MEP(trlsLow), '*')
end


% if trial number is unequal, number of observations will be unequal too
% e.g.: trials = 57
% length(MEP_low) >> 29
% length(MEP_high) >> 28