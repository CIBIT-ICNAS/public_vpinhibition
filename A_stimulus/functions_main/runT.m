function [ OutputRunT ] = runT( run_name , S , T )
%RUNAB123 Stimulus function for run B
% Usage: [ OutputRunB ] = runAB( run_name , S , T )
%
% Inputs:
%   : run_name - string with run name and number
%   : S - struct containing very important info (screen, subject, colors,
%   trigger, response box, eyetracker, etc)
%   : T - struct containg the textures
%
% Outputs:
%   : OutputRunB - struct containing time of start, end and trigger of the
%   stimulation, key presses and codes, user responses.
%

OutputRunT = struct(); % Output struct
D = struct(); % Dots struct

% Screen('Preference', 'SkipSyncTests', 1);

%% Clean the COMs
% Never trust what happened before...
if S.RSPBOX
    IOPort('Purge',S.response_box_handle);
end
if S.TRIGGER
    IOPort('Purge', S.syncbox_handle);
end

%% Read PRT
load(fullfile(S.input_path,['Protocols_' run_name '.mat']));

%% Search for Storage and Black condition
% During these conditions, only the fixation cross is drawn.
% strCondIdx = find(ismember(condNames, 'Storage')) ;
% blackCondIdx = find(ismember(condNames, 'Black')) ;

%% Find Report Condition
% drawFixationCross needs it to turn orange
% S.reportCondIdx = find(ismember(condNames, 'Report')) ;

%% Define Attention Task
% During the trial period (adapt+test) the fixation cross can
% increase its size from 1 to 4 times. In this section are defined the
% frames in which this happens.
% 
% eventsAtt = []; % volumes in which the cross increases.
% nEvents = zeros(size(intervalsPRT.Report,1),1); % Events (number of times the cross increased)
% nEventsUser = zeros(size(intervalsPRT.Report,1),2); % User responses. This will be compared to nEvents to determine the feedback.
% 
% ReportBlockDur = intervalsPRT.Report(1,2)-intervalsPRT.Report(1,1)+1;
% TestBlockDur = intervalsPRT.Test_aPatt(1,2)-intervalsPRT.Test_aPatt(1,1)+1;
% 
% framesQ1= zeros(size(intervalsPRT.Report,1)-1,ReportBlockDur*S.fps*S.TR); % Frames during which the user responses will be interpreted.
% 
% framesEvents = zeros(nFrames,1); % Frames in which the cross increases.
% supersizeTime = 0.25; % In seconds
% 
% for j = 1:size(intervalsPRT.Report,1) % Iterate on the Report blocks 
%     
%     int = intervalsPRT.Static(j,2)+2 : intervalsPRT.Report(j,1)-2; % volumes during which the cross can increase.
%     
%     testStartVol = intervalsPRT.Report(j,1)-2 - TestBlockDur; % test block start volume.
% 
%     nEvents(j) = ceil(4*rand(1,1));
%     auxVols = int(randperm(length(int), nEvents(j)));
%     
%     % In the case of 4 events, check if all of them are outside the test
%     % block. If they are, move the last of them to the test block. This
%     % avoids loss of attention during the trial (the participant knows a
%     % priori that the cross can only increase a maximum of 4 times).
%     if max(auxVols) < testStartVol && nEvents(j) == 4 
%         auxVols(auxVols==max(auxVols)) =  testStartVol + 1;
%         fprintf('[runB] Displaced one event to Test block - volume %i.\n',testStartVol + 1);
%     end
%     
%     eventsAtt = [eventsAtt auxVols] ;
%     
%     framesQ1(j,:) = (intervalsPRT.Report(j,1)-1)*S.fps*S.TR + 1 : intervalsPRT.Report(j,2)*S.fps*S.TR;
%     
% end
% 
% % -- Create framesEvents based on eventsAtt and supersizeTime
% for j = 1:length(eventsAtt)
%     
%     framesEvents((eventsAtt(j)-1)*S.fps*S.TR+1 : (eventsAtt(j)-1)*S.fps*S.TR+1 + ceil(supersizeTime*S.fps)) = 1;
%     
% end

%% Initialise key-related stuff
% nFrames = size(nFrames,1);
keysPressed = zeros(nFrames , 2);
nButtons = 2;
key_codes = zeros(nButtons,1);
KbName('UnifyKeyNames');

%% Stim

try
    % -- Open Window
    [ windowID , winRect ] = Screen('OpenWindow', S.screenNumber , S.screenBackground );
    
    % -- Rendering Options
    Screen('BlendFunction', windowID, GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    
    % -- Flip Interval
    ifi = Screen('GetFlipInterval', windowID);
    
    % -- Select specific text font, style and size
    Screen('TextFont', windowID , 'Arial');
    Screen('TextSize', windowID , 20);
    
    % -- Center Variables
    [S.xCenter, S.yCenter] = RectCenter(winRect);
    
    % -- Build Dots
    [ D ] = buildDots( D , T );
    
    % -- Hide Cursor
    HideCursor;
    
    % -- Initialise EyeTracker
    if S.EYETRACKER
        eyeTrackerInit( windowID , S , run_name(end-1:end));
    end
    
    % -- Instruction
    text = 'Press any key to continue!';
    DrawFormattedText(windowID, text, 'center', 'center', S.white);
    Screen('Flip',windowID);
    KbStrokeWait;

    % -- Map Keys / Identify Key Codes
    btUnique = false;
    
    while ~btUnique
        
        % Button 1
        DrawFormattedText(windowID, 'Press Button 1 (Coherent)', 'center', 'center', S.white);
        Screen('Flip',windowID);
        if S.RSPBOX
            pr = 1;
            while pr
                key = IOPort('Read',S.response_box_handle);
                if ~isempty(key) && (length(key) == 1)
                    key_codes(1) = key;
                    pr = 0;
                end
                IOPort('Flush',S.response_box_handle);
            end
        else
            [~,code] = KbPressWait;
            key_codes(1) = find(code==1);
        end
        
        % Button 2
        DrawFormattedText(windowID, 'Press Button 2 (Incoherent)', 'center', 'center', S.white);
        Screen('Flip',windowID);
        if S.RSPBOX
            pr = 1;
            while pr
                key = IOPort('Read',S.response_box_handle);
                if ~isempty(key) && (length(key) == 1)
                    key_codes(2) = key;
                    pr = 0;
                end
                IOPort('Flush',S.response_box_handle);
            end
        else
            [~,code] = KbPressWait;
            key_codes(2) = find(code==1);
        end
        
        if length(unique(key_codes)) == nButtons
            btUnique = true;
        else
            key_codes = zeros(nButtons,1);
        end
        
        clear code key pr
        
    end
    
    % -- Start Eyetracker
    if S.EYETRACKER
        eyeTrackerStart();
    end
    
    % -- Wait for Key Press or Trigger
    if S.TRIGGER
        DrawFormattedText(windowID, 'Waiting to start...', 'center', 'center', S.white);
        Screen('Flip',windowID);
        disp('[runT] Waiting for trigger...')
        
        [gotTrigger, timeStamp] = waitForTrigger(S.syncbox_handle,1,300); % timeOut = 5 min (300s)
        if gotTrigger
            disp('[runT] Trigger Received.')
            IOPort('Flush', S.syncbox_handle);
            IOPort('Purge', S.syncbox_handle);
        else
            disp('[runT] Trigger Not Received. Aborting!')
            return
        end
    else
        DrawFormattedText(windowID, 'Press Enter to Start', 'center', 'center', S.white);
        Screen('Flip',windowID);
        disp('[runT] Waiting to start...')
        KbPressWait;
    end
    
    % -- Start parameters
    t = 0;
    textIndexX = 1;
    textIndexY = 1;
%     idxStop = 1;
    
    % -- Frame Iteration
    vbl = Screen('Flip', windowID);
    disp('[runT] Starting iteration...')
    
    init = GetSecs;
    
    while t < nFrames % Iteration on the frames
        
        % Make Texture
        windowtext = Screen('MakeTexture', windowID, T.Textures{textIndexY,textIndexX});

%         if framesCond(t+1) ~= strCondIdx && framesCond(t+1) ~= blackCondIdx % Not Storage or Black
            
            % Draw Lines
            Screen('DrawTextures', windowID, windowtext)
            
            % Draw Dots if not ambigous condition
            if any(framesCond(t+1) == [5,6,7])
                no_dots_flag = true;
            else
                no_dots_flag = false;
            end
            [ D ] = drawDots( windowID , D , T , S , textIndexX , textIndexY , no_dots_flag);
            
            % Update Texture Index
            if any(framesCond(t+1) == [2,3,5,6,7]) || any(framesDots(t+1) == [22,33])
                textIndexY = textIndexY + 1;
            elseif any(framesDots(t+1) == [222,333])
                textIndexY = textIndexY - 1;
            elseif any(framesDots(t+1) == [2222,3333])
                textIndexX = textIndexX - 1;
            elseif any(framesDots(t+1) == [22222,33333])
                textIndexX = textIndexX + 1;
            end
            
            % Restrict Texture Index
            if textIndexX > T.nTextX
                textIndexX = 1;
            elseif textIndexX <= 0
                textIndexX = T.nTextX;
            end
            if textIndexY > T.nTextY
                textIndexY = 1;
            elseif textIndexY <= 0
                textIndexY = T.nTextY;
            end

%         end % End If Not Storage
        
%         if any(t+1 == framesQ1(idxStop,:)) % Ask for Number of times
%             
%             % Fixation Cross
%             drawFixationCross( windowID , S , framesCond(t+1) , framesEvents(t+1) );
%             
%             % -------- KEYS --------
%             [keyPress,~,keyCode] = KbCheck();
%             if keyCode(KbName('escape')) == 1 %Quit if "Esc" is pressed
%                 throw(MException('user:escape','Aborted by escape key.'))
%             end
%             
%             if S.RSPBOX
%                 [key,timestamp] = IOPort('Read',S.response_box_handle);
%                 
%                 if ~isempty(key) && (length(key) == 1)
%                     IOPort('Flush',S.response_box_handle);
%                     
%                     keysPressed(t+1,1) = key;
%                     nEventsUser(idxStop,1) = find(key_codes==key);
%                     nEventsUser(idxStop,2) = timestamp;
%                     
%                 end
%                 
%                 IOPort('Flush',S.response_box_handle);
%             else
%                 if keyPress
%                     a = find(keyCode==1);
%                     if length(a) == 1
%                         keysPressed(t+1,1) = a;
%                         if any(a==key_codes) % User is numb and can press any key...
%                             nEventsUser(idxStop,1) = find(key_codes==a);
%                         end
%                         nEventsUser(idxStop,2) = idxStop;
%                     end
%                 end
%             end
%             
%             keysPressed(t+1,2) = GetSecs;
%             
%             % If the user has answered
%             if nEventsUser(idxStop,1) ~= 0
%                 if nEventsUser(idxStop,1) == nEvents(idxStop) % Correct answer
%                     drawFixationCross( windowID , S , 100 , framesEvents(t+1) );
%                 else % Incorrect answer
%                     drawFixationCross( windowID , S , 101 , framesEvents(t+1) );
%                 end
%             end
%             
%             % Do it
%             vbl = Screen('Flip', windowID, vbl + 0.5*ifi);
%             
%             % Reaching the end of this response trial, increase idxStop
%             if ( t+1 == framesQ1(idxStop,end) ) && ( idxStop < size(framesQ1,1) )
%                 idxStop = idxStop+1;
%             end
%             
%         else % Trial itself (do not record key responses)
%             
            % Fixation Cross
            drawFixationCross( windowID , S );
            
            % Do it
            vbl = Screen('Flip', windowID, vbl + 0.5*ifi);
%             
            % -------- KEYS --------
            [keyPress,~,keyCode] = KbCheck();
            if keyPress
                if keyCode(KbName('escape')) == 1 %Quit if "Esc" is pressed
                    throw(MException('user:escape','Aborted by escape key.')) 
                else
                    a = find(keyCode==1);
                    if length(a) == 1
                        keysPressed(t+1,1) = a;
                    end
                end
            end
            
            % Clean the Response Box COM at every iteration. This allows
            % the participant to freely play with the buttons when there
            % should be no responses. But just don't try it.
            if S.RSPBOX
                IOPort('Purge',S.response_box_handle);
            end
            
            keysPressed(t+1,2) = GetSecs;
%             
%         end % End of If Ask for number of times
        
%         % Record frames
%         rect = round([S.xCenter-S.height/2 ; S.yCenter-S.height/2 ; S.xCenter+S.height/2 ; S.yCenter+S.height/2]);
%         frameImage=Screen('GetImage', windowID, rect, [], [], []);
%         imwrite(frameImage,[pwd '\output_frames\RunB_' num2str(t+1000) '.png'],'png');
        
        % Clear Texture
        Screen('Close', windowtext);
        
        % Move dots
        [ D ] = moveDots( framesDots(t+1) , D , T );
        
        % -- Iterate frame
        t = t+1;
        
    end % End of frame interation
    
    finit = GetSecs;
    
    % -- Stop EyeTracker
    if S.EYETRACKER
        eyeTrackerStop( 0 );
    end
    
    % -- Close window
    Screen('CloseAll');
    ShowCursor;
    commandwindow;
    
    % -- Attention Task Grade
%     correct_incorrect = [ sum(nEvents==nEventsUser(:,1))  sum(nEvents~=nEventsUser(:,1))  ];
    
    % -- Export Log
    OutputRunT.Subject = S.SUBJECT;
    OutputRunT.Start = init;
    OutputRunT.End = finit;
    if S.TRIGGER
        OutputRunT.TriggerTime = timeStamp;
    end
    OutputRunT.Keys = keysPressed;
%     OutputRunT.FramesQ1 = framesQ1;
    OutputRunT.KeyCodes = key_codes;
%     OutputRunT.EventsVol = eventsAtt;
%     OutputRunT.Events = nEvents;
%     OutputRunT.EventsUser = nEventsUser;
%     OutputRunT.Grade = correct_incorrect;
    
    output_filename = [S.SUBJECT '_' run_name '_' datestr(now,'HHMM_ddmmmmyyyy')];
    save(fullfile(S.output_path,output_filename),'OutputRunT')
    
    disp('[runT] Done.')
    
catch ME
    
    finit = GetSecs;
    
    % -- Stop EyeTracker
    if S.EYETRACKER
        eyeTrackerStop( 1 );
    end
    
    % -- Close window
    Screen('CloseAll');
    ShowCursor;
    commandwindow;
    
    % -- Deal with it
    switch ME.identifier
        case 'user:escape'
            OutputRunT.Subject = S.SUBJECT;
            OutputRunT.Start = init;
            OutputRunT.End = finit;
            OutputRunT.Keys = keysPressed;
%             OutputRunT.FramesQ1 = framesQ1;
            OutputRunT.KeyCodes = key_codes;
%             OutputRunT.EventsVol = eventsAtt;
%             OutputRunT.Events = nEvents;
%             OutputRunT.EventsUser = nEventsUser;
            
            output_filename = [S.SUBJECT '_' run_name '_' datestr(now,'HHMM_ddmmmmyyyy') '_Aborted'];
            save(fullfile(S.output_path,output_filename),'OutputRunT')
            
            disp('[runT] Aborted by escape key.')
        otherwise
            rethrow(ME);
            % psychrethrow(psychlasterror);
    end
    
end % End try/catch

end % End function