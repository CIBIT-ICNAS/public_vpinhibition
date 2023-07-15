clear,clc

%% Folder and map files
keypressFolder = fullfile('..','data','expC-raw-keypress-output');

D = dir(fullfile(keypressFolder,'*.mat'));

%% Initialize matrix
% subs, 6 runs
performanceMatrix = zeros(length(D),6);

%% Extract data from all subs and runs
for ss = 1:length(D)

    load(fullfile(D(ss).folder,D(ss).name))

    for rr = 1:6

        performanceMatrix(ss,rr) = Output.(sprintf('RunMRI_D12_R%i',rr)).Grade(1) / 12 * 100;

    end

end

%% Averages
M = mean(performanceMatrix,2);

fprintf('Mean attention task accuracy M = %0.2f SD = %0.2f\n',mean(M),std(M))
% Mean attention task accuracy M = 91.67 SD = 5.89
