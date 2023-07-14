clear, clc;

addpath(fullfile('..','tools','edf-converter'))

%% Configs

% --- Configuration data
load('Configs_VP_INHIBITION.mat')

% --- Retrieve participants with eyetracker data
subjectsList = datasetConfigs.subjects(logical(datasetConfigs.eyetracker));
nSubjects = length(subjectsList);
nRuns = 6;

% --- I/O Folders
inputFolder = fullfile('..','data','expC-raw-eyetracker');
outputFolder = fullfile('..','data','expC-eyetracker');

%%
f_percond_duration = struct();  
conditionNames = {'Coh_aCoh','Coh_aInCoh','Coh_aNA','InCoh_aCoh','InCoh_aInCoh','InCoh_aNA'};

%% Iterate
for ss = 1:nSubjects
    
    auxDir = dir(fullfile(inputFolder,subjectsList{ss},'*.edf'));
            
    auxArray = zeros(length(auxDir),1);

    % Iterate on the edf files
    for rr = 1:length(auxDir) % six runs
        
        fname = fullfile(auxDir(rr).folder,auxDir(rr).name);
        myDataInfo = Edf2Mat(fname); % Convert edf to mat

        %% syncronization eyetracker <-> MRI
        % since the start of the eyetracker was not syncronized with the start of
        % the MRI, we must cut according to the end of the recording (which was
        % syncronized with the end of the MRI acquisition)
        fs = 1000; % sampling frequency of the eyetracker signal in Hz
        duration_expected = 374 * fs;
        duration_actual = myDataInfo.timeline(end) - myDataInfo.timeline(1);
        
        extra_timepoints = duration_actual-duration_expected;
        %first_timepoint = extra_timepoints(end)+1; % the first timepoint that should be considered the start of the acquisition.
        
        %% Import run protocol
        protocolFolder = fullfile('..','data','input-stimulus');
        protocolFile = load(fullfile(protocolFolder,sprintf('Protocols_RunMRI_D12_R%i.mat',rr)),'intervalsPRT');
        
        
        % Convert to eyetracer timepoints
        for jj = 1:6 % six conditions
           
            protocolFile.timepoints.(conditionNames{jj}) = (protocolFile.intervalsPRT.(conditionNames{jj}) - 1) * fs;
        
        end
        
        %% Read fixation events        
        f_start = myDataInfo.Events.Efix.start - myDataInfo.timeline(1) + extra_timepoints;
        f_end = myDataInfo.Events.Efix.end - myDataInfo.timeline(1) + extra_timepoints;
       
        for jj = 1:6 % six conditions
        
           cond_range = [ protocolFile.timepoints.(conditionNames{jj})(1,1):protocolFile.timepoints.(conditionNames{jj})(1,2) protocolFile.timepoints.(conditionNames{jj})(2,1):protocolFile.timepoints.(conditionNames{jj})(2,2) ];
                
           inter = ismember(f_start, cond_range); % find fixation events that started during this condition
        
           % find events that end after the end of the condition to trim its
           % duration
           inter_end = f_end(inter);
        
           over_ind_1 = (inter_end > protocolFile.timepoints.(conditionNames{jj})(1,2)) & (inter_end < protocolFile.timepoints.(conditionNames{jj})(2,1)); % bigger than the end of this block and smaller than the start of the following
           over_ind_2 = inter_end > protocolFile.timepoints.(conditionNames{jj})(2,2);
        
           inter_end(over_ind_1) = protocolFile.timepoints.(conditionNames{jj})(1,2);
           inter_end(over_ind_2) = protocolFile.timepoints.(conditionNames{jj})(2,2);
        
           % recalculate durations
           inter_durations = inter_end - f_start(inter);
        
           f_percond_duration.(conditionNames{jj})(ss,rr) = sum(inter_durations) / (6*2*fs) * 100;
        
        end
    end
end

%% Stat analysis

[h1,p1] = ttest(sum(f_percond_duration.Coh_aCoh,2),sum(f_percond_duration.InCoh_aCoh,2))

[h2,p2] = ttest(sum(f_percond_duration.Coh_aInCoh,2),sum(f_percond_duration.InCoh_aInCoh,2))


