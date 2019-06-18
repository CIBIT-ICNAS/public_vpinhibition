function [ framesCond , framesDots , nFrames , condNames , nCond , intervalsPRT ] = extractFramesPRT( path , name , TR , fps )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

[ condNames , intervalsPRT , intervals ] = readProtocol( path , name , TR );

nVols = length(intervals);
nFrames = nVols*fps*TR;
nCond = length(condNames);

framesCond = zeros(nFrames,1);
framesDots = zeros(nFrames,1);

% Create intervals for dots where NonAdapt changes at every volume
intervalsDots = intervals;
idx = find(intervalsDots == 4); % doublecheck if NonAdapt is 4th condition

intervalsDots(idx(1:6:end)) = 2;       % Pattern Down
intervalsDots(idx(2:6:end)) = 333;     % Component Out Vertical
intervalsDots(idx(3:6:end)) = 22;      % Pattern Up
intervalsDots(idx(4:6:end)) = 3333;    % Component In Vertical
intervalsDots(idx(5:6:end)) = 3;       % Component Out
intervalsDots(idx(6:6:end)) = 33;      % Component In
% intervalsDots(idx(7:8:end)) = 222;     % Pattern Right
% intervalsDots(idx(8:8:end)) = 2222;    % Pattern Left

for t = 0:nVols-1
    
    framesCond(t*fps*TR+1:(t+1)*fps*TR) = intervals(t+1);
    
    if any(intervalsDots(t+1) == 2)                    % Pattern Down
        framesDots(t*fps*TR+1:(t+1)*fps*TR) = 22;
        
    elseif any(intervalsDots(t+1) == 3)                % Component Out
        framesDots(t*fps*TR+1:(t+1)*fps*TR) = 33;
        
    elseif any(intervalsDots(t+1) == 22)               % Pattern Up
        framesDots(t*fps*TR+1:(t+1)*fps*TR) = 222;
        
    elseif any(intervalsDots(t+1) == 33)               % Component In
        framesDots(t*fps*TR+1:(t+1)*fps*TR) = 333;
        
    elseif any(intervalsDots(t+1) == 222)              % Pattern Right
        framesDots(t*fps*TR+1:(t+1)*fps*TR) = 2222;
        
    elseif any(intervalsDots(t+1) == 333)              % Component Out Vertical
        framesDots(t*fps*TR+1:(t+1)*fps*TR) = 3333;
        
    elseif any(intervalsDots(t+1) == 2222)             % Pattern Left
        framesDots(t*fps*TR+1:(t+1)*fps*TR) = 22222;
        
    elseif any(intervalsDots(t+1) == 3333)             % Component In Vertical
        framesDots(t*fps*TR+1:(t+1)*fps*TR) = 33333;
        
    elseif any(intervalsDots(t+1) == [5,6,7])             % Ambiguous (no dots)
        framesDots(t*fps*TR+1:(t+1)*fps*TR) = 55;
        
    else                                               % Static Texture
        framesDots(t*fps*TR+1:(t+1)*fps*TR) = 11;
        
    end
    
end

end