clear,clc

%% Folders
bidsFolder     = '/DATAPOOL/VPINHIBITION/BIDS-VP-INHIBITION';
derivFolder    = fullfile(bidsFolder,'derivatives');
backupFolder   = fullfile(derivFolder,'backupRois');
codeFolder     = pwd;

%% Extract Subject List from BIDS
aux = dir(fullfile(bidsFolder,'sub-*'));
subjectList = extractfield(aux,'name');
nSubjects = length(subjectList);
clear aux

%% Init backup folder
if ~exist(backupFolder,'dir')
   mkdir(backupFolder); 
end

%% Iterate on the subjects and runs

for ss = 1:nSubjects
    
    mkdir(fullfile(backupFolder,subjectList{ss}, 'model_task-inhib_run-1_MNI152NLin2009cAsym'))

    copyfile(fullfile(derivFolder,'spm12', subjectList{ss}, 'model_task-inhib_run-1_MNI152NLin2009cAsym','VOI_leftMT_loc_1.mat'),...
             fullfile(backupFolder,subjectList{ss}, 'model_task-inhib_run-1_MNI152NLin2009cAsym','VOI_leftMT_loc_1.mat'));

    copyfile(fullfile(derivFolder,'spm12', subjectList{ss}, 'model_task-inhib_run-1_MNI152NLin2009cAsym','VOI_rightMT_loc_1.mat'),...
             fullfile(backupFolder,subjectList{ss}, 'model_task-inhib_run-1_MNI152NLin2009cAsym','VOI_rightMT_loc_1.mat'));
    
end
