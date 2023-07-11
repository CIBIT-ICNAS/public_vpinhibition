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
condNames = {'Static','Coherent','Incoherent','NonAdapt','Ambiguous_aCoh','Ambiguous_aInCoh','Ambiguous_aNA','MAE'};

blockDuration = [ 6 6 6 6 6 6 6 3 ]; %in volumes (think TR = 1000ms)

blockColor = [170 170 170 ; 150 210 240 ; 150 170 200 ; 150 130 160 ; ...
                            200 170 100 ; 200 170 120 ; 200 170 140 ; ...
                            120 120 120];

PRTParameters.nCond = length(condNames);

PRTConditions = struct();

for c = 1:PRTParameters.nCond
    
    PRTConditions.(condNames{c}).Color = blockColor(c,:);
    PRTConditions.(condNames{c}).BlockDuration = blockDuration(c);
    PRTConditions.(condNames{c}).Intervals = [];
    PRTConditions.(condNames{c}).NumBlocks = 0;
    
end

%% Run Behav D6 R1
%SEQ = [ 1 2 5 6 1 3 5 6 1 4 5 6 1 2 5 6 1 3 5 6 1 4 5 6 1 2 5 6 1 3 5 6 1 4 5 6 1 2 5 6 1 3 5 6 1 4 5 6];
SEQ = [ 1 2 5 8 1 3 6 8 1 4 7 8 1 2 5 8 1 3 6 8 1 4 7 8 1 2 5 8 1 3 6 8 1 4 7 8 1 2 5 8 1 3 6 8 1 4 7 8];

[ PRTConditions_R1 ] = buildIntervals( SEQ , PRTConditions );

generatePRT( PRTParameters , PRTConditions_R1 , 'RunBehav_D6_R1' , outputFolder );

%% Run Behav D6 R2
%SEQ = [ 1 2 5 6 1 4 5 6 1 3 5 6 1 2 5 6 1 4 5 6 1 3 5 6 1 2 5 6 1 4 5 6 1 3 5 6 1 2 5 6 1 4 5 6 1 3 5 6 ];
SEQ = [ 1 2 5 8 1 4 7 8 1 3 6 8 1 2 5 8 1 4 7 8 1 3 6 8 1 2 5 8 1 4 7 8 1 3 6 8 1 2 5 8 1 4 7 8 1 3 6 8 ];

[ PRTConditions_R2 ] = buildIntervals( SEQ , PRTConditions );

generatePRT( PRTParameters , PRTConditions_R2 , 'RunBehav_D6_R2' , outputFolder );

%% Run Behav D6 R3
%SEQ = [ 1 3 5 6 1 4 5 6 1 2 5 6 1 3 5 6 1 4 5 6 1 2 5 6 1 3 5 6 1 4 5 6 1 2 5 6 1 3 5 6 1 4 5 6 1 2 5 6 ];
SEQ = [ 1 3 6 8 1 4 7 8 1 2 5 8 1 3 6 8 1 4 7 8 1 2 5 8 1 3 6 8 1 4 7 8 1 2 5 8 1 3 6 8 1 4 7 8 1 2 5 8 ];

[ PRTConditions_R3 ] = buildIntervals( SEQ , PRTConditions );

generatePRT( PRTParameters , PRTConditions_R3 , 'RunBehav_D6_R3' , outputFolder );
