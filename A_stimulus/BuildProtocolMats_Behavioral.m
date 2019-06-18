% clear,clc;

addpath('functions')

TR = 1; % in seconds
fps = 60;

prt_folder = fullfile(pwd,'prt');
output_folder = fullfile(pwd,'input');

%% Run Behav D6 R1
[ framesCond , framesDots , nFrames , condNames , nCond , intervalsPRT ] = extractFramesPRT( prt_folder , 'RunBehav_D6_R1.prt' , TR , fps );
save(fullfile(output_folder,'Protocols_RunBehav_D6_R1.mat'),'framesCond','framesDots','nFrames','condNames','intervalsPRT','nCond');

%% Run Behav D6 R2
[ framesCond , framesDots , nFrames , condNames , nCond , intervalsPRT ] = extractFramesPRT( prt_folder , 'RunBehav_D6_R2.prt' , TR , fps );
save(fullfile(output_folder,'Protocols_RunBehav_D6_R2.mat'),'framesCond','framesDots','nFrames','condNames','intervalsPRT','nCond');

%% Run Behav D6 R3
[ framesCond , framesDots , nFrames , condNames , nCond , intervalsPRT ] = extractFramesPRT( prt_folder , 'RunBehav_D6_R3.prt' , TR , fps );
save(fullfile(output_folder,'Protocols_RunBehav_D6_R3.mat'),'framesCond','framesDots','nFrames','condNames','intervalsPRT','nCond');

%% Run Behav D12 R1
[ framesCond , framesDots , nFrames , condNames , nCond , intervalsPRT ] = extractFramesPRT( prt_folder , 'RunBehav_D12_R1.prt' , TR , fps );
save(fullfile(output_folder,'Protocols_RunBehav_D12_R1.mat'),'framesCond','framesDots','nFrames','condNames','intervalsPRT','nCond');

%% Run Behav D12 R2
[ framesCond , framesDots , nFrames , condNames , nCond , intervalsPRT ] = extractFramesPRT( prt_folder , 'RunBehav_D12_R2.prt' , TR , fps );
save(fullfile(output_folder,'Protocols_RunBehav_D12_R2.mat'),'framesCond','framesDots','nFrames','condNames','intervalsPRT','nCond');

%% Run Behav D12 R3
[ framesCond , framesDots , nFrames , condNames , nCond , intervalsPRT ] = extractFramesPRT( prt_folder , 'RunBehav_D12_R3.prt' , TR , fps );
save(fullfile(output_folder,'Protocols_RunBehav_D12_R3.mat'),'framesCond','framesDots','nFrames','condNames','intervalsPRT','nCond');

%% End
disp('Protocols .mat created.')
