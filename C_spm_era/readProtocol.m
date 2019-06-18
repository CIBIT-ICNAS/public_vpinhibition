function [ cond_names , intervalsPRT , intervals , baseCondIndex , colors , blockDurations , blockNumber ] = readProtocol( path , name , TR )
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here

Temp_prtFile = BVQXfile(fullfile(path , name));

type = Temp_prtFile.ResolutionOfTime;

switch type
    case 'msec'
        prtFile = Temp_prtFile.ConvertToVol(TR*1000);
%         prtFile.SaveAs(fullfile(path , [name(1:end-4) '_vol.prt']));
    case 'Volumes'
        prtFile = Temp_prtFile;
end

cond_names = prtFile.ConditionNames;

baseCondIndex = 0;
colors = zeros(length(cond_names),3);
blockDurations = zeros(length(cond_names),1);
blockNumber = zeros(length(cond_names),1);
v = 1;

for cond = 1:length(cond_names)
    
    cond_names{cond} = strtrim(cond_names{cond});
    
    if strcmpi(cond_names{cond},'baseline') || strcmpi(cond_names{cond},'Neutral')
        baseCondIndex = cond;
        intervals = ones(prtFile.Cond(baseCondIndex).OnOffsets(end),1);
    end
    colors(cond,:) = prtFile.Cond(cond).Color;

end

for cond = 1:length(cond_names)
    
    try
        intervalsPRT.(cond_names{cond}) = prtFile.Cond(cond).OnOffsets;
    catch %Invalid condition name for struct
        cond_names{cond} = ['C' cond_names{cond}];
        fprintf('Renaming Condition name to %s\n',cond_names{cond});
        intervalsPRT.(cond_names{cond}) = prtFile.Cond(cond).OnOffsets;
    end
    
    for int = 1:size(intervalsPRT.(cond_names{cond}),1)
        
        intervals(intervalsPRT.(cond_names{cond})(int,1) : intervalsPRT.(cond_names{cond})(int,2) ) = v;
        
    end
    
    blockDurations(cond) = intervalsPRT.(cond_names{cond})(end,2) - intervalsPRT.(cond_names{cond})(end,1) + 1;
    blockNumber(cond) = size(intervalsPRT.(cond_names{cond}),1);
    v = v + 1;
    
end

end


%%
%         switch type
%             case 'msec'
%                 mod1 = mod(intervalsPRT.(cond_names{cond})(int,1),TR*1000);
%                 if mod1 > 0
%                     start = ceil(intervalsPRT.(cond_names{cond})(int,1) / (TR*1000));
%                 else
%                     start = floor(intervalsPRT.(cond_names{cond})(int,1) / (TR*1000)) + 1;
%                 end
%
%                 mod2 = mod(intervalsPRT.(cond_names{cond})(int,2),TR*1000);
%                 if mod2 > 0
%                     finish = ceil(intervalsPRT.(cond_names{cond})(int,2) / (TR*1000));
%                 else
%                     finish = floor(intervalsPRT.(cond_names{cond})(int,2) / (TR*1000));
%                 end
%
%                 intervals( start : finish ) = v;
%
%             case 'Volumes'
%                 intervals(intervalsPRT.(cond_names{cond})(int,1) : intervalsPRT.(cond_names{cond})(int,2) ) = v;
%         end

