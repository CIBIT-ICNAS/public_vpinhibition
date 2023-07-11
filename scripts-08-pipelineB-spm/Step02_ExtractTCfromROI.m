%% -- Step02_ExtractTCfromROI.m ---------------------------------------- %%
% ----------------------------------------------------------------------- %
% Script for executing the second step regarding the SPM analysis of fMRI
% data after preprocessing with fmriPrep v21 - extracting the time course
% of a number of ROIs.
%
% Dataset:
% - Inhibition (Visual Perception)
%
% Warnings:
% - a number of values/steps are custom for this dataset - full code review
% is strongly advised for different datasets
% - this was designed to run on sim01 - a lot of paths must change if run
% at any other computer
% - only compatible with space MNI152NLin2009cAsym
%
% Requirements:
% - SPM12 in path
%
% Author: Alexandre Sayal
% CIBIT, University of Coimbra
% February 2022
% ----------------------------------------------------------------------- %
% ----------------------------------------------------------------------- %
%% Load Packages on sim01

% SPM12
addpath('/SCRATCH/software/toolboxes/spm12')

%% Folders
bidsFolder     = '/DATAPOOL/VPINHIBITION/BIDS-VP-INHIBITION';
derivFolder    = fullfile(bidsFolder,'derivatives');
codeFolder     = pwd;

%% Extract Subject List from BIDS
aux = dir(fullfile(bidsFolder,'sub-*'));
subjectList = extractfield(aux,'name');
nSubjects = length(subjectList);
clear aux

%% Extract Task List from BIDS
aux = dir(fullfile(bidsFolder,'sub-01','ses-01','func','*_task-*_bold.json'));
taskList = extractfield(aux,'name');
clear aux

taskList = cellfun(@(x) x(15:end-10), taskList, 'un', 0); % remove trailing and leading info (VERY custom)
nTasks = length(taskList);

%% ===================================================================== %%
% Manually define the ROIs for the first run. In this case, I defined left
% and right hMT+ based on the peak voxel coordinates. The ROIs are placed
% on /derivatives/spm12/sub-**/model_task-inhib_run-1_MNI152NLin2009cAsym/
% and named VOI_leftMT_loc_1.mat. The peak coordinates are saved in the xY
% matrix.
% ===================================================================== %%

%% ROI List
roiList = {'leftMT','rightMT','leftSPL','rightSPL'};
nROIs = length(roiList);

%% Start coord and TC table
COORD = struct();
TC = struct();
nVolumes = 374;

for rr = 1:nROIs
    
    COORD.(roiList{rr}) = zeros(nSubjects,3);
    TC.(roiList{rr}) = zeros(nSubjects,nTasks,nVolumes);
    
end

%% Iterate on the subjects, ROIs, and runs

for ss = 1:nSubjects % iterate on the subjects
    
    
    for rr = 1:nROIs % iterate on the ROIs
        roiName = roiList{rr};
        
        % Fetch coordinates
        load(fullfile(derivFolder,'spm12', subjectList{ss}, 'model_task-inhib_run-1_MNI152NLin2009cAsym',['VOI_' roiName '_loc_1.mat']))
        COORD.(roiName)(ss,:) = xY.xyz;
        
        clear xY Y
        
        for tt = 1:nTasks % iterate on the runs
            
            clear matlabbatch
            
            %% leftMT
            matlabbatch{1}.spm.util.voi.spmmat = {fullfile(derivFolder,'spm12', subjectList{ss}, ['model_' taskList{tt} '_MNI152NLin2009cAsym'],'SPM.mat')};
            matlabbatch{1}.spm.util.voi.adjust = 5; % index of F-contrast with the effects of interest (all other will be regressed-out from the signal)
            matlabbatch{1}.spm.util.voi.session = 1;
            matlabbatch{1}.spm.util.voi.name = sprintf('%s_run-%i', roiName, tt);
            
            matlabbatch{1}.spm.util.voi.roi{1}.spm.spmmat = {''};
            matlabbatch{1}.spm.util.voi.roi{1}.spm.contrast = 4;
            matlabbatch{1}.spm.util.voi.roi{1}.spm.conjunction = 1;
            matlabbatch{1}.spm.util.voi.roi{1}.spm.threshdesc = 'FWE';
            matlabbatch{1}.spm.util.voi.roi{1}.spm.thresh = 0.05;
            matlabbatch{1}.spm.util.voi.roi{1}.spm.extent = 25;
            matlabbatch{1}.spm.util.voi.roi{1}.spm.mask = struct('contrast', {}, 'thresh', {}, 'mtype', {});
            
            matlabbatch{1}.spm.util.voi.roi{2}.sphere.centre = COORD.(roiName)(ss,:);
            matlabbatch{1}.spm.util.voi.roi{2}.sphere.radius = 5;
            matlabbatch{1}.spm.util.voi.roi{2}.sphere.move.fixed = 1;
            matlabbatch{1}.spm.util.voi.expression = 'i1 & i2';
            
            % RUN
            spm_jobman('run', matlabbatch);
            
            % Check if empty (no voxels found inside the ROI)
            if isempty(Y)
                warning('\n!!!\n\nReducing threshold to find active voxels on %s, run %s, %s!\n\n!!!',subjectList{ss},taskList{tt},roiName)
                
                matlabbatch{1}.spm.util.voi.roi{1}.spm.spmmat = {''};
                matlabbatch{1}.spm.util.voi.roi{1}.spm.contrast = 4;
                matlabbatch{1}.spm.util.voi.roi{1}.spm.conjunction = 1;
                matlabbatch{1}.spm.util.voi.roi{1}.spm.threshdesc = 'none';
                matlabbatch{1}.spm.util.voi.roi{1}.spm.thresh = 0.001;
                matlabbatch{1}.spm.util.voi.roi{1}.spm.extent = 5;
                matlabbatch{1}.spm.util.voi.roi{1}.spm.mask = struct('contrast', {}, 'thresh', {}, 'mtype', {});
                
                % RUN
                spm_jobman('run', matlabbatch);
            end
            
            if isempty(Y)
                warning('\n!!!\n\nFilling with NaNs on %s, run %s, %s!\n\n!!!',subjectList{ss},taskList{tt},roiName)
                Y = nan(nVolumes,1);
            end
            
            % Save
            TC.(roiName)(ss,tt,:) = Y;
            clear xY Y
            
        end % end run iteration
        
    end % end ROI iteration
    
end % end subject iteration

%% Export
cd(codeFolder)
save('TimeSeries_PeakCoordinates.mat')
