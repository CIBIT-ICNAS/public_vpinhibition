% clear,clc;

addpath('functions')

TR = 1; % in seconds
fps = 60;

prt_folder = fullfile(pwd,'prt');
output_folder = fullfile(pwd,'input');

nRuns = 6;

for rr = 1:nRuns
    
    [ framesCond , framesDots , nFrames , condNames , nCond , intervalsPRT ] = extractFramesPRT_MRI( prt_folder , ['RunMRI_D12_R' num2str(rr) '.prt'] , TR , fps );
    save(fullfile(output_folder,['Protocols_RunMRI_D12_R' num2str(rr) '.mat']),'framesCond','framesDots','nFrames','condNames','intervalsPRT','nCond');
    
end

%% Clear
disp('Protocols .mat created.')
