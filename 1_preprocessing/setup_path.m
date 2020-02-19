%% setup path
format compact
restoredefaultpath
projectdir = pwd;
resourcePath = fullfile(projectdir, 'Resources');
addpath(resourcePath);
load('DRCMR_EEG_Head.mat');
eeglayout = lay;
%plot the layout if wanted
%cfg = [];
%cfg.layout = eeglayout;   % this is the layout structure that you created before
%ft_layoutplot(cfg);

ftPath = fullfile(resourcePath, 'fieldtrip-20190724');
addpath(ftPath);
ft_defaults;
%sasicaPath = fullfile(resourcePath, 'SASICA-feature-ft_compat');
%addpath(sasicaPath);

addpath(genpath(fullfile(projectdir, 'code'))); % so functions and scripts are callable

%% participant list

IDs = {'X32870' 'X03496' 'X08629' 'X08695' 'X12768' 'X28575' 'X48262' 'X51474' 'X61679' 'X67535' 'X85256'}; % other subjects have to be added here later on

%% folder structure for project

baseDir='data';

indiDir='individual';
groupDir='group';
figDir='figures';

% individual folders
rawDir='raw'; 
preprocDir='preproc'; preStimDir='preStim'; postStimDir='postStim';
resultDir='results'; 

% create folders
% parent folders
if ~exist(fullfile(baseDir, indiDir))
    mkdir(fullfile(baseDir, indiDir));
end

if ~exist(fullfile(baseDir, groupDir))
    mkdir(fullfile(baseDir, groupDir));
end

if ~exist(fullfile(baseDir, figDir))
    mkdir(fullfile(baseDir, figDir));
end

% children folders
for iSub=1:length(IDs)
    
    if ~exist(fullfile(baseDir, indiDir, IDs{iSub}, preprocDir))
        mkdir(fullfile(baseDir, indiDir, IDs{iSub}, preprocDir))
    end
    
    if ~exist(fullfile(baseDir, indiDir, IDs{iSub}, preprocDir, preStimDir))
        mkdir(fullfile(baseDir, indiDir, IDs{iSub}, preprocDir, preStimDir))
    end
    
    if ~exist(fullfile(baseDir, indiDir, IDs{iSub}, preprocDir, postStimDir))
        mkdir(fullfile(baseDir, indiDir, IDs{iSub}, preprocDir, postStimDir))
    end  
    
        
    if ~exist(fullfile(baseDir, indiDir, IDs{iSub}, resultDir))
        mkdir(fullfile(baseDir, indiDir, IDs{iSub}, resultDir))
    end
    
    if ~exist(fullfile(baseDir, indiDir, IDs{iSub}, rawDir))
        mkdir(fullfile(baseDir, indiDir, IDs{iSub}, rawDir))
    end
    
      
       
end
    
    





