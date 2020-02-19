% 1. merge EEG datasets
% 2. merge rec EMG datasets
% 3. merge EEG & rec EMG
% 4. merge unrec EMG datasets
%% 1. merge EEG data sets 
load(fullfile(sb_preStimDir, [num2str(IDs{iSub}) '_preproc_EEG_data1.mat']), 'dataEEG_vrejec1');
load(fullfile(sb_preStimDir, [num2str(IDs{iSub}) '_preproc_EEG_data2.mat']), 'dataEEG_vrejec2');

cfg = [];
cfg.keepsampleinfo = 'yes';
datasets_EEG = ft_appenddata(cfg, dataEEG_vrejec1, dataEEG_vrejec2);

wantFigures = 0;
if wantFigures == 1
    figure, subplot(2,1,1)
    plot(dataEEG_vrejec2.time{1}, dataEEG_vrejec2.trial{1} (17,:))
    subplot(2,1,2)
    plot(datasets_EEG.time{1},datasets_EEG.trial{57} (17,:)) 
end

%% 2. merge rectified preTMS EMG data sets
load(fullfile(sb_preStimDir, [num2str(IDs{iSub}) '_preproc_EMG_data1.mat']), 'dataEMG_rejec1');
load(fullfile(sb_preStimDir, [num2str(IDs{iSub}) '_preproc_EMG_data2.mat']), 'dataEMG_rejec2');

datasets_EMG = ft_appenddata(cfg, dataEMG_rejec1, dataEMG_rejec2);

if wantFigures == 0
    figure, subplot(2,1,1)
    plot(dataEMG_rejec2.time{1}, dataEMG_rejec2.trial{1} )
    subplot(2,1,2)
    plot(datasets_EMG.time{1},datasets_EMG.trial{57} ) 
end

%% 3. merge EEG and rectified preTMS EMG datasets
% adjust times
datasets_EMG.time = datasets_EEG.time; 

% append both data sets
data_preproc = ft_appenddata(cfg, datasets_EEG, datasets_EMG);

data_preproc=ft_struct2single(data_preproc); % make it smaller for disk storage
save(fullfile(sb_preStimDir, [num2str(IDs{iSub}) '_preTMS_EEG_EMG.mat']), 'data_preproc');

%% 4. merge unrectified postTMS EMG data sets
load(fullfile(sb_postStimDir, [num2str(IDs{iSub}) '_UnRec_preproc_EMG_data1.mat']), 'UnRec_dataEMG_rejec1');
load(fullfile(sb_postStimDir, [num2str(IDs{iSub}) '_UnRec_preproc_EMG_data2.mat']), 'UnRec_dataEMG_rejec2');

UnRec_datasets_EMG = ft_appenddata(cfg, UnRec_dataEMG_rejec1, UnRec_dataEMG_rejec2);

save(fullfile(sb_postStimDir, [num2str(IDs{iSub}) '_postTMS_EMG.mat']), 'UnRec_datasets_EMG');

if wantFigures == 1
    figure, subplot(2,1,1)
    plot(UnRec_dataEMG_rejec2.time{1}, UnRec_dataEMG_rejec2.trial{1} )
    title('Compate unrectified postTMS EMG')
    subplot(2,1,2)
    plot(UnRec_datasets_EMG.time{61},UnRec_datasets_EMG.trial{61} ) 
end


