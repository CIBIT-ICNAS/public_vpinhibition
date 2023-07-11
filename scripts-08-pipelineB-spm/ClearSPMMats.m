clear,clc

%% Folders
bidsFolder     = '/DATAPOOL/VPINHIBITION/BIDS-VP-INHIBITION';
derivFolder    = fullfile(bidsFolder,'derivatives');

%% Extract Subject List from BIDS
aux = dir(fullfile(bidsFolder,'sub-*'));
subjectList = extractfield(aux,'name');
clear aux

%% Extract Task List from BIDS
aux = dir(fullfile(bidsFolder,'sub-01','ses-01','func','*_task-*_bold.json'));
taskList = extractfield(aux,'name');

taskList = cellfun(@(x) x(15:end-10), taskList, 'un', 0); % remove trailing and leading info (VERY custom)
nTasks = length(taskList);

%% Iterate on subs and runs
for ss = 1:length(subjectList)
    
    subjectID = subjectList{ss};
    
    fprintf('%s processing started!\n',subjectID)
    
    %% Create spmFolder
    spmFolder = fullfile(derivFolder,'spm12',subjectID);
    
    for tt = 1:nTasks
        
        delete(fullfile(spmFolder,['model_' taskList{tt} '_MNI152NLin2009cAsym'],'SPM.mat'))
            
    end
    
    % delete task-all
    delete(fullfile(spmFolder,'model_task-inhib_run-all_MNI152NLin2009cAsym','SPM.mat'))

    
end
