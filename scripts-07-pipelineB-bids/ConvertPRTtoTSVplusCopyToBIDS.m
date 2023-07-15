
% Requires neuroelf

tsvFolder = 'temp-tsv';

bidsFolder = '/media/alexandresayal/DATA4TB/BIDS-VP-INHIBITION';

%% Subject List
aux = dir(fullfile(bidsFolder,'sub-*'));
subjectList = extractfield(aux,'name');
clear aux

nSubjects = length(subjectList);

%% Task List
taskList = {'task-inhib_run-1','task-inhib_run-2','task-inhib_run-3','task-inhib_run-4','task-inhib_run-5','task-inhib_run-6'};

%% PRT List
% In this case the protocol is the same for all participants
prtFolder = '/media/alexandresayal/HDD4TB/RAW_DATA_VP_INHIBITION/PRTs_CrossInhibition';
prtList = {'run1.prt','run2.prt','run3.prt','run4.prt','run5.prt','run6.prt'};
nRuns = length(prtList);

%% TR
TR = 1;

%% Convert

for rr = 1:nRuns
    
    [ cond_names , intervalsPRT ,~,~,~, blockDur, blockNum ] = readProtocol( fullfile(prtFolder, prtList{rr}) , TR );
    
    trial_type = {};
    onset = [];
    duration = [];
    for cc = 1:length(cond_names)
        trial_type = [trial_type ; repmat({cond_names(cc)},blockNum(cc),1)];
        onset = [onset ; intervalsPRT.(cond_names{cc})(:,1).*TR-TR];
        duration = [duration ; repmat(blockDur(cc).*TR,blockNum(cc),1)];
    end
    [onset,idx] = sort(onset);
    trial_type = trial_type(idx);
    duration = duration(idx);
    
    T = table(onset,duration,trial_type);
    
    export_file = fullfile(tsvFolder,...
        sprintf('task-%s_events.txt',prtList{rr}(1:end-4)));
    
    writetable(T,export_file,'Delimiter','\t');
    movefile(export_file,[export_file(1:end-4) '.tsv']);
    
end

%% Copy to BIDS
% Iterate

for ss = 1:nSubjects
    
    subjectID = subjectList{ss};
    
    for rr = 1:nRuns
        
        subfuncFolder = fullfile(bidsFolder,subjectID,'ses-01','func');
        
        tsvBIDSName = sprintf('%s_ses-01_%s_events.tsv',subjectID,taskList{rr});
        
        % Check if exists to replace
        if exist(fullfile(subfuncFolder,tsvBIDSName),'file')
            
            copyfile(fullfile(tsvFolder,['task-' prtList{rr}(1:end-4) '_events.tsv']),...
                     fullfile(subfuncFolder,tsvBIDSName) )
            
        else
            warning('%s does not exist. Expected?',tsvBIDSName)
        end
        
    end
    
    fprintf('%s done! \n',subjectID)
       
end
