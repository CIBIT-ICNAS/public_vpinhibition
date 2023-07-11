%% -- Step04_ExtractTvalues.m ---------------------------------------- %%
% ----------------------------------------------------------------------- %
% Script for executing the forth step regarding the SPM analysis of fMRI
% data after preprocessing with fmriPrep v21 - extracting the T value of
% the various contrasts of interests for the ROIs defined before.
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
% - shadedErrorBar
%
% Author: Alexandre Sayal
% CIBIT, University of Coimbra
% February 2022
% ----------------------------------------------------------------------- %
% ----------------------------------------------------------------------- %

clear,clc

%% Folders
bidsFolder     = '/DATAPOOL/VPINHIBITION/BIDS-VP-INHIBITION';
derivFolder    = fullfile(bidsFolder,'derivatives');
% fmriPrepFolder = fullfile(bidsFolder,'derivatives','fmriprep');
% codeFolder     = pwd;

%% Extract Subject List from BIDS
aux = dir(fullfile(bidsFolder,'sub-*'));
subjectList = extractfield(aux,'name');
nSubjects = length(subjectList);
clear aux

%% ROI List
roiList = {'leftMT','rightMT','leftSPL','rightSPL'};
nROIs = length(roiList);

%% Contrast List
contrastList = [1:4 6:11];
nContrasts = length(contrastList);

%% Start structure to save data
% mind that I am adding 1 to the number of contrasts because contrast
% number 5 is an F contrast - i do not want to mess up with the indexes
ROIdata = zeros(nSubjects,nContrasts+1,nROIs);

%% Iterate on the subjects, ROIs, contrasts

for ss = 1:nSubjects % iterate on the subjects
    
    roiFolder = fullfile(derivFolder,'spm12',subjectList{ss},'model_task-inhib_run-1_MNI152NLin2009cAsym');
    multiRunFolder = fullfile(derivFolder,'spm12',subjectList{ss},'model_task-inhib_run-all_MNI152NLin2009cAsym');
    
    for cc = 1:nContrasts % iterate on the contrasts
        
        contrastFile = fullfile(multiRunFolder,sprintf('con_%04d.nii',contrastList(cc)));
        
        for rr = 1:nROIs % iterate on the ROIs
            
            roiFile = fullfile(roiFolder,sprintf('VOI_%s_run-1_mask.nii',roiList{rr}));
            
            ROIdata(ss,contrastList(cc),rr) = extractROIdata(roiFile,contrastFile);
            
        end % end ROI iteration
          
    end % end contrast iteration
    
end % end subject iteration

%%
figure
notBoxPlot([ROIdata(:,1,1) ROIdata(:,2,1) ROIdata(:,3,1) ROIdata(:,4,1)])
