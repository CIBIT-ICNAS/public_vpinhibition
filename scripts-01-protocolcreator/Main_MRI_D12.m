clear, clc;

outputFolder = fullfile('..','data','prt');

%% PRT Parameters
PRTParameters = struct();

PRTParameters.FileVersion = 2;
PRTParameters.Resolution = 'Volumes';
PRTParameters.ExperimentName = 'VisualPerception_CrossInhibition';
PRTParameters.BackgroundColor = [0 0 0];
PRTParameters.TextColor = [255 255 255];
PRTParameters.TimeCourseColor = [1 1 1];
PRTParameters.TimeCourseThick = 3;
PRTParameters.ReferenceFuncColor = [0 0 80];
PRTParameters.ReferenceFuncThick = 2;

%% PRT Conditions
condNames = {'Static','Coherent','Incoherent','NonAdapt',...
    'Coh_aCoh','Coh_aInCoh','Coh_aNA',...
    'InCoh_aCoh','InCoh_aInCoh','InCoh_aNA',...
    'MAE','Report','Discard'};

blockDuration = [ 6 12 12 12 6 6 6 6 6 6 3 3 4 ]; %in volumes (think TR = 1000ms)

blockColor = [170 170 170 ; 0 115 190 ; 216 83 25 ; 236 177 32 ; ...
                            110 115 190 ; 150 83 25 ; 190 177 32 ; ...
                            150 115 190 ; 100 83 25 ; 120 177 32 ; ...
                            120 160 120; 120 220 120; 50 50 50];

PRTParameters.nCond = length(condNames);

PRTConditions = struct();

for c = 1:PRTParameters.nCond
    
    PRTConditions.(condNames{c}).Color = blockColor(c,:);
    PRTConditions.(condNames{c}).BlockDuration = blockDuration(c);
    PRTConditions.(condNames{c}).Intervals = [];
    PRTConditions.(condNames{c}).NumBlocks = 0;
    
end

%% Trials
trials = [ 1 2 5 11 12 ;
           1 2 8 11 12 ;
           1 3 6 11 12 ;
           1 3 9 11 12 ;
           1 4 7 11 12 ;
           1 4 10 11 12 ];

%% Run MRI D12 R1
SEQ = [13 reshape(trials([1 3 2 6 5 4 1 3 2 6 5 4],:)',1,60) 1 13];

[ PRTConditions_R1 ] = buildIntervals( SEQ , PRTConditions );

generatePRT( PRTParameters , PRTConditions_R1 , 'RunMRI_D12_R1' , outputFolder );

%% Run MRI D12 R2
SEQ = [13 reshape(trials([1 3 5 2 4 6 1 3 5 2 4 6],:)',1,60) 1 13];

[ PRTConditions_R2 ] = buildIntervals( SEQ , PRTConditions );

generatePRT( PRTParameters , PRTConditions_R2 , 'RunMRI_D12_R2' , outputFolder );

%% Run MRI D12 R3
SEQ = [13 reshape(trials([5 3 1 6 4 2 5 3 1 6 4 2],:)',1,60) 1 13];

[ PRTConditions_R3 ] = buildIntervals( SEQ , PRTConditions );

generatePRT( PRTParameters , PRTConditions_R3 , 'RunMRI_D12_R3' , outputFolder );

%% Run MRI D12 R4
SEQ = [13 reshape(trials([4 5 6 2 3 1 4 5 6 2 3 1],:)',1,60) 1 13];

[ PRTConditions_R4 ] = buildIntervals( SEQ , PRTConditions );

generatePRT( PRTParameters , PRTConditions_R4 , 'RunMRI_D12_R4' , outputFolder );

%% Run MRI D12 R5
SEQ = [13 reshape(trials([1 6 2 5 3 4 1 6 2 5 3 4],:)',1,60) 1 13];

[ PRTConditions_R5 ] = buildIntervals( SEQ , PRTConditions );

generatePRT( PRTParameters , PRTConditions_R5 , 'RunMRI_D12_R5' , outputFolder );

%% Run MRI D12 R6
SEQ = [13 reshape(trials([6 1 5 2 4 3 6 1 5 2 4 3],:)',1,60) 1 13];

[ PRTConditions_R6 ] = buildIntervals( SEQ , PRTConditions );

generatePRT( PRTParameters , PRTConditions_R6 , 'RunMRI_D12_R6' , outputFolder );
