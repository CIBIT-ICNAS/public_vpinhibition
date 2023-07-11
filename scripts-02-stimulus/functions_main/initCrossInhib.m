function [ S ] = initCrossInhib( S , sNumber )
%INITADAPTCHECK Initialise parameters for the cross ihnibition protocol.
% Usage: [ S ] = initCrossInhib( S , sNumber );
%
% Inputs:
%   : S - struct containing very important info (screen, subject, colors,
%   trigger, response box, eyetracker, etc)
%   : sNumber - screen number (if set to 50, the screen with the highest
%   index is selected)
% Outputs:
%   : S - same as input, a lot fatter :)
%

% --- Define screen number based on input
screens=Screen('Screens');
if sNumber == 50
    S.screenNumber=max(screens);
else
  S.screenNumber=sNumber;  
end

% --- Settings
Screen('Preference', 'SkipSyncTests', 0);
KbName('UnifyKeyNames');

% --- Determine screen size and set colors
[S.screenX, S.screenY] = Screen('WindowSize', S.screenNumber);

S.white = WhiteIndex(S.screenNumber);
S.black = BlackIndex(S.screenNumber);
S.grey = S.white / 2;

S.screenBackground = 0;
S.textBackground = 50;
S.lines = 130;

% --- Stimulus/Experiment Settings
S.TR = 1; % Repetition time in seconds
% S.fps = Screen('FrameRate',S.screenNumber); % Screen Frame Rate
S.fps = 60; % force 60Hz
S.dist = 156;  % Distance from eye to screen in cm
S.width = 70; % Width of the screen in cm

% S.reportCondIdx = 99;

% --- Select Texture Size in pixels
S.factor = S.screenX / 1920; % relative to fullHD screen
% S.height = 749; % Corresponds to 10 deg in visual angle
% S.height = 674; % Corresponds to 9 deg in visual angle
% S.height = 598; % Corresponds to 8 deg in visual angles
S.height = round(674 * S.factor); % Corresponds to 9 deg in visual angle 

% --- Folders
S.input_path = fullfile('..','data','input-stimulus');
S.output_path = fullfile('..','data','output-stimulus');
%S.outputImages_path = fullfile(pwd,'outputImages');

% --- Open COM Ports for Response box and Trigger
if S.RSPBOX
    S.response_box_handle = IOPort('OpenSerialPort','COM3');
    IOPort('Purge',S.response_box_handle);
end

if S.TRIGGER
    S.syncbox_handle = IOPort('OpenSerialPort', 'COM2', 'BaudRate=57600 DataBits=8 Parity=None StopBits=1 FlowControl=None');
    IOPort('Purge',S.syncbox_handle);
end

end % End function