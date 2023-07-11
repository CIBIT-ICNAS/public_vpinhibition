clear,clc

%% Requires:
% - bidsphysio (https://github.com/cbinyu/bidsphysio) - install with `python3 -m pip install bidsphysio`
% - make sure to start matlab from the terminal where physio2bidsphysio is
% available (in the path)
% - jsonlab

addpath(fullfile('..','tools','jsonlab'))

%% Settings
subjectList = 1:20;

bidsFolder = '/media/alexandresayal/DATA4TB/BIDS-VP-INHIBITION';

nPhysioType = 3; % number of physio .log files per run (e.g. _RESP, _PULS, _Info) - only works for 3 types here

%% Iterate

for ss = 1:length(subjectList)
    
    subID = sprintf('%02i',subjectList(ss));
    rawPhysioFolder = ['/media/alexandresayal/HDD4TB/RAW_DATA_VP_INHIBITION/VPIS' subID '_LOGS'];
    
    %% Fetch runs and acquisition times
    
    D1 = dir(fullfile(bidsFolder,['sub-' subID],'ses-01','func','*_bold.json'));
    
    J = struct();
    
    for ii = 1:length(D1)
        
        aux = loadjson( fullfile(D1(ii).folder,D1(ii).name) );
        
        J(ii).name = D1(ii).name(1:end-10);
        J(ii).time = aux.AcquisitionTime;
        
    end
    
    [~,index] = sortrows({J.time}.'); J = J(index); clear index; % sort by time
    
    nRuns = length(J);
    
    %% Fetch physiofiles
    
    D2 = dir( fullfile(rawPhysioFolder,'*.log') );
    
    % sort list anyway
    [~,index] = sortrows({D2.name}.'); D2 = D2(index); clear index;
    
    nPhysioFiles = length(D2);
    
    %% Check if number of runs match
    
    if nPhysioFiles/nPhysioType ~= nRuns
        error('Mismatch between the number of func and physio files')
    end
    
    %% Convert and Copy files to BIDS directory
    
    for rr = 1:nRuns
        
        idx = 1+(rr-1)*nPhysioType;
        
        cmd1 = sprintf('physio2bidsphysio --infiles %s %s %s --bidsprefix %s -v',...
            fullfile(rawPhysioFolder,D2(idx).name),...
            fullfile(rawPhysioFolder,D2(idx+1).name),...
            fullfile(rawPhysioFolder,D2(idx+2).name),...
            fullfile(bidsFolder,['sub-' subID],'ses-01','func',J(rr).name) );
        
        system(cmd1);
        
    end
    
    fprintf('All physio files created for subject %s\n',subID)
    
end


