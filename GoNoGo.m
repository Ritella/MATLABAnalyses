% Codes for a Go NoGo Response Inhibition Task
% Includes infrequent Go sounds so that Go and NoGo sounds
% are matched on novelty.
% Uses non-verbalizable sounds that are equally disimilar to 
% each other in order to ensure that language is not employed
% and to ensure that both the novelty consideration is 
% enforced and that the NoGo sound is not too similar to a Go sound
% RL
% 06/16

function [] = GoNoGo (subNum, isBlind, testOrPracticeLetter, runNum, scanNum) 

%% CLEAN UP
clearvars -except subNum isBlind runNum scanNum testOrPracticeLetter
PsychPortAudio('Close');% Close all audio handles
CedrusResponseBox('CloseAll');
C.possibleCedrusPorts = {'/dev/cu.usbserial-000012FD', '/dev/cu.usbserial-000013FA', '/dev/cu.usbserial-142',  '/dev/cu.usbserial-141'}; % Check in terminal, 'ls /dev' 
C.versionNum = 2;

%% DECLARE CONSTANTS THAT MAY VARY
% There are 4 Practice Runs and unlimited Task Runs, 
% but typically 6 are run
if strcmp(testOrPracticeLetter, 'P') && runNum==1
        C.percentFreqGo = 0.50;  % Check that this adds to 0 and gives trials that are a natural number
        C.percentInfreqGo = 0.25;
        C.percentNoGo = 0.25;
        C.giveFeedback = 1;
        C.soundDuration = 0.45; % in seconds
        C.soundITI = 0.35; 
        C.feedbackDuration = 0.5; % this should be the max of the duration for all the feedback sounds 
        C.feedbackITI = 0.2;
        C.trialsPerRun = 80;
        C.restPeriods = 0;
        C.trigger = 1;
        C.soundWait = 0.1; % Has to be less than C.soundITI or C.feedbackITI or sound won't close
elseif strcmp(testOrPracticeLetter, 'P') && runNum==2
        C.percentFreqGo = 0.50;  % Check that this adds to 0 and gives trials that are a natural number
        C.percentInfreqGo = 0.25;
        C.percentNoGo = 0.25;
        C.giveFeedback = 1;
        C.soundDuration = 0.45; % in seconds
        C.soundITI = 0.25; 
        C.feedbackDuration = 0.5;  % this should be the max of the duration for all the feedback sounds 
        C.feedbackITI = 0.2;
        C.trialsPerRun = 80;
        C.restPeriods = 0;
        C.trigger = 1;
        C.soundWait = 0.1; % Has to be less than C.soundITI or C.feedbackITI or sound won't close
elseif strcmp(testOrPracticeLetter, 'P') && runNum == 3
        C.percentFreqGo = 0.50;  % Check that this adds to 0 and gives trials that are a natural number
        C.percentInfreqGo = 0.25;
        C.percentNoGo = 0.25;
        C.giveFeedback = 1;
        C.soundDuration = 0.45; % in seconds
        C.soundITI = 0.15; 
        C.feedbackDuration = 0.5; % this should be the max of the duration for all the feedback sounds 
        C.feedbackITI = 0.2;
        C.trialsPerRun = 80;
        C.restPeriods = 0;
        C.trigger = 1;
        C.soundWait = 0.1; % Has to be less than C.soundITI or C.feedbackITI or sound won't close
elseif strcmp(testOrPracticeLetter, 'P') && runNum == 4
        C.percentFreqGo = 0.50;  % Check that this adds to 0 and gives trials that are a natural number
        C.percentInfreqGo = 0.25;
        C.percentNoGo = 0.25;
        C.giveFeedback = 1;
        C.soundDuration = 0.45; % in seconds
        C.soundITI = 0.15; 
        C.feedbackDuration = 0.5; % this should be the max of the duration for all the feedback sounds 
        C.feedbackITI = 0.2;
        C.trialsPerRun = 80;
        C.restPeriods = 0;
        C.trigger = 0;
        C.soundWait = 0.1; % Has to be less than C.soundITI or C.feedbackITI or sound won't close
elseif strcmp(testOrPracticeLetter, 'T')
        C.percentFreqGo = 0.50;  % Check that this adds to 0 and gives trials that are a natural number
        C.percentInfreqGo = 0.25;
        C.percentNoGo = 0.25;
        C.giveFeedback = 0;
        C.soundDuration = 0.45; % in seconds
        C.soundITI = 0.5; 
        C.feedbackDuration = 0;
        C.feedbackITI = 0;
        C.trialsPerRun = 400;
        C.restPeriods = 1;
        C.trigger = 1;
        C.soundWait = 0.1; % Has to be less than C.soundITI or C.feedbackITI or sound won't close
else
        error('Run Number not valid for TestOrPractice Flag');
end

%% DECLARE CONSTANTS THAT DONT VARY ACCORDING TO CONTEXT
C.restWait = 20;
C.falseStartTime = 0.15;
C.freqGoNum = 1; C.infreqGoNum = -1; C.noGoNum = 0;
C.maxReps = 3; % No more than 3 InfreqGo or NoGo in a row
C.typesToLimit = [C.infreqGoNum; C.noGoNum];
C.restBeforeTrials = [1,round(C.trialsPerRun/3),round(C.trialsPerRun*2/3),C.trialsPerRun+1];

%% SETUP STRUCTS
out = struct('runNum',runNum,'subNum',subNum,'runStart', [], 'trialNumber', zeros(1,C.trialsPerRun), 'trialType', zeros(1,C.trialsPerRun), ...
    'trialSoundStart',zeros(1,C.trialsPerRun), 'keyPress', zeros(1,C.trialsPerRun),  ...
    'RT', zeros(1,C.trialsPerRun), 'trialEnd', zeros(1,C.trialsPerRun), 'feedbackSoundStart',zeros(1,C.trialsPerRun), ...
    'keyPressInterim', [], 'keyPressFalseStarts', [], 'keyPressOnset', zeros(1,C.trialsPerRun), 'keyPressTrialExtras', [],  'keyPressFeedback', [], 'tempTimeCheck', [], 'pressCheck', []);  % cell arrays can be added later but cannot be initialized
out.currentSound = {}; % cell arrays can be added later but cannot be initialized
out.currentTrialType = {};
out.trialSoundPlayed = {};
out.feedbackSoundPlayed = {};
stats = struct;
C.exptStr = 'GNGC';
if isBlind == 0
    C.subPre = 'S';
else
    C.subPre = 'CB';
end
if strcmp(testOrPracticeLetter, 'T')
    C.runType = 'R'; % R stands for Run
elseif strcmp(testOrPracticeLetter, 'P')
    C.runType = 'Practice';
else
    error('Please enter T or P for Task or Practice, respectively.');
end
out.subName = [C.exptStr '_' C.subPre '_' sprintf('%02d',subNum)];
out.taskName = 'gonogo';
out.runName = [C.runType, sprintf('%02d', runNum)];
out.matName = strcat('Data/', out.subName, '/', out.subName, '-', out.taskName, '_', out.runName, '-D', datestr(now, 'mmdd'), '.mat');
if exist(out.matName, 'file')
    error('File with these input parameters already exists!');
end
if ~exist('scanNum','var')
    out.scanNum = 0;
else
    out.scanNum = scanNum;
    clear scanNum;
end
out.testOrPractice = testOrPracticeLetter;
out.audioMatName = strcat('Data/', out.subName, '/', out.subName, '-', out.taskName, '_', 'sounds', '-D', datestr(now, 'mmdd'), '.mat');
out.audioFileTypes =load(out.audioMatName, 'AUDIO_FILE_TYPES');
out.audioFileTypes = out.audioFileTypes.AUDIO_FILE_TYPES;

%% PRELIMINARY SETUP
warning off;
rand('twister',sum(100*clock));
AssertOpenGL;

%% CONDITION LOGIC
C.numberFreqGo = C.trialsPerRun*C.percentFreqGo;
C.numberInfreqGo = C.trialsPerRun*C.percentInfreqGo;
C.numberNoGo = C.trialsPerRun*C.percentNoGo;
C.isGo = vertcat(C.freqGoNum*ones(C.numberFreqGo,1), C.infreqGoNum*ones(C.numberInfreqGo,1), C.noGoNum*ones(C.numberNoGo,1));
C.isGoPermute = permuteWithMaxReps(C.isGo, C.maxReps, C.typesToLimit);
[~, C.isGoPermuteTranspose] = sort(C.isGoPermute); %necessary to get the transpose or the indexing will be the opposite of how you want it.
out.trialType = C.isGo(C.isGoPermuteTranspose,:)'; %Add other to trialArray or rename?
Priority(9);

%% LOAD AND INITIALIZE AUDIO FILES   
for i = 1:length(out.audioFileTypes) %for each filetype
    filetype = out.audioFileTypes(i).type;
    filename = out.audioFileTypes(i).name;
    audioFiles.(filetype) = dir(filename);  %create a structure containing all associated files        
    for j = 1:length(audioFiles.(filetype)) %for each file, load and get info
        [audioFiles.(filetype)(j).wavedata, audioFiles.(filetype)(j).freq] = audioread(fullfile(fileparts(filename), audioFiles.(filetype)(j).name)); % load sound file
        audioFiles.(filetype)(j).duration = length(audioFiles.(filetype)(j).wavedata) ./ audioFiles.(filetype)(j).freq;
        audioFiles.(filetype)(j).channels = size(audioFiles.(filetype)(j).wavedata,2);
        audioFiles.(filetype)(j).pahandle = PsychPortAudio('Open', [], [], 2, audioFiles.(filetype)(j).freq, audioFiles.(filetype)(j).channels, 0); % opens sound buffer at a different frequency
        PsychPortAudio('FillBuffer', audioFiles.(filetype)(j).pahandle, audioFiles.(filetype)(j).wavedata'); % loads data into buffer
    end
    clear filetype filename j;
end
clear i;
InitializePsychSound(1); %initializes sound driver...the 1 pushes for low latency   

%%Setup Trials to Save Time 
for i = 1:C.trialsPerRun
    switch out.trialType(i)
        case C.freqGoNum           
            out.currentTrialType{i} = 'Go1';
        case C.infreqGoNum 
            out.currentTrialType{i} = 'Go2';
        case C.noGoNum
            out.currentTrialType{i} = 'NoGo';
    end
    out.currentSound{i} = audioFiles.(out.currentTrialType{i}).name;
end
clear i;

%% CONNECT TO CEDRUS BUTTON BOX
display('Connecting to Response Box');
try
    C.serial = instrhwinfo ('serial');
    portIdx = find(ismember(C.possibleCedrusPorts, C.serial.AvailableSerialPorts));
    C.cedrusPort = C.possibleCedrusPorts{portIdx};
    C.cedrusHandle = CedrusResponseBox('Open', C.cedrusPort);
    CedrusResponseBox('SetConnectorMode', C.cedrusHandle, 'ReflectiveSinglePulse');
    clear portIdx
catch e
    sca;
    rethrow(e) % Try restarting Matlab after connecting and/or restarting box
end

%% WAIT FOR TRIGGER
display('Start the fMRI Scan Now');
switch C.trigger
    case 0
        display('No Trigger Needed.');
        WaitSecs(5);
        CedrusResponseBox('ClearQueues', C.cedrusHandle);
        CedrusResponseBox('FlushEvents', C.cedrusHandle);
        out.runStart = GetSecs;
    case 1
        display('Waiting for Trigger to Start Experiment');
        CedrusResponseBox('ClearQueues', C.cedrusHandle);
        CedrusResponseBox('FlushEvents', C.cedrusHandle);
        while 1
            evt = CedrusResponseBox('WaitButtons', C.cedrusHandle);
            out.runStart = GetSecs;
            if ~isempty(evt) && evt.button == 7 %sync pulse happens
                CedrusResponseBox('FlushEvents', C.cedrusHandle);
                clear evt;
                break;
            end
        end 
end

%% RUN TRIALS
disp('Experiment Starting!');
for currentTrial = 1:C.trialsPerRun+1 
    if C.restPeriods==1 && ismember(currentTrial,C.restBeforeTrials)
        tempTime = CedrusResponseBox('ResetRTTimer', C.cedrusHandle) - out.runStart; 
        WaitSecs(C.restWait); 
        % Get additional button presses (and flush queue) before the next trial  
        while 1
            evt = CedrusResponseBox('GetButtons', C.cedrusHandle);
            if isempty(evt)
                clear tempTime;
                break;
            else
                if evt.action == 1
                    out.keyPressInterim(end+1) = evt.rawtime + tempTime;  % can't use this method all the time b/c there is a button-box time drift.
                    %disp('key press during interim');
                end
                clear evt;
            end  
        end
    end
    if currentTrial==C.trialsPerRun+1 % Need the loop to run through trials + 1
        break;
    end
    stopSound = 0;
    stopFeedbackSound = 0;
    out.trialNumber(currentTrial) = currentTrial;  % This needs to be kept inside this loop in order to count trials actually implemented 
    [stimHandle, out.trialSoundStart(currentTrial), out.trialSoundPlayed{currentTrial}] = playSound(audioFiles, out.currentTrialType{currentTrial}, 0, out.runStart);
    while GetSecs - out.runStart - out.trialSoundStart(currentTrial) < C.soundDuration + C.soundITI % Warning: This setup means that the trial can extend slightly beyond its intended duration
        if(GetSecs - out.runStart - out.trialSoundStart(currentTrial) >= C.soundDuration + C.soundWait && stopSound==0); PsychPortAudio('Stop', stimHandle); stopSound=1; clear ans; end;
        evt=CedrusResponseBox('GetButtons', C.cedrusHandle);         
        if ~isempty(evt) && evt.action == 1    
            evt.time = evt.ptbfetchtime - out.runStart;
            out.pressCheck(end+1) = GetSecs - out.runStart;
            if evt.time - out.trialSoundStart(currentTrial) < C.falseStartTime
                out.keyPressFalseStarts(end+1) = evt.time;
                %disp('false start!');
            else
                if out.keyPress(currentTrial) == 0 % If this is the first response per trial
                    out.keyPressOnset(currentTrial) = evt.time; 
                    out.keyPress(currentTrial) = 1;
                    out.RT(currentTrial) = out.keyPressOnset(currentTrial) - out.trialSoundStart(currentTrial);
                    %disp('key press');
                else
                    out.keyPressTrialExtras(end+1) = evt.time;
                    %disp('extra key press');
                end
           end
           clear evt;
        end
    end
    switch C.giveFeedback
        case 1
            switch out.keyPress(currentTrial)
                case 1       
                    switch out.trialType(currentTrial)
                        case C.noGoNum
                            [feedbackHandle, out.feedbackSoundStart(currentTrial), out.feedbackSoundPlayed{currentTrial}] = playSound(audioFiles, 'Wrong', 0, out.runStart);
                        otherwise
                            [feedbackHandle, out.feedbackSoundStart(currentTrial), out.feedbackSoundPlayed{currentTrial}] = playSound(audioFiles, 'Click', 0, out.runStart); % Only play a click if it's within the allocated time
                    end
                case 0
                    switch out.trialType(currentTrial)
                        case C.noGoNum
                            [feedbackHandle, out.feedbackSoundStart(currentTrial), out.feedbackSoundPlayed{currentTrial}] = playSound(audioFiles, 'Correct',0,out.runStart);
                        otherwise
                            [feedbackHandle, out.feedbackSoundStart(currentTrial), out.feedbackSoundPlayed{currentTrial}] = playSound(audioFiles, 'Timeout',0,out.runStart); 
                    end
            end
    end
    while GetSecs - out.runStart - out.trialSoundStart(currentTrial) < C.soundDuration + C.soundITI + C.feedbackDuration + C.feedbackITI
        if(GetSecs - out.runStart - out.feedbackSoundStart(currentTrial) >= C.feedbackDuration + C.soundWait && stopFeedbackSound==0); PsychPortAudio('Stop', feedbackHandle); stopFeedbackSound=1; clear ans; end;   
        evt=CedrusResponseBox('GetButtons', C.cedrusHandle);         
        if ~isempty(evt) && evt.action == 1
            evt.time = evt.ptbfetchtime - out.runStart;
            out.keyPressFeedback(end+1) = evt.time;
            clear evt;
            %disp('key press during feedback');
        end
    end
    clear stopSound stopFeedbackSound stimHandle feedbackHandle;
    out.trialEnd(currentTrial) = GetSecs - out.runStart;  
end
out.finalTrial = currentTrial;
clear runNum subNum currentTrial evt isBlind testOrPracticeLetter ans;
out.runDuration = GetSecs - out.runStart;
save(out.matName);  % Don't save every trial to save time

%% GET STATISTICS
[stats.FreqGo, stats.InfreqGo, stats.NoGo, stats.FalseStarts] = GetStats(out.matName);
display('Done!');
disp(['Go 1 = ' mat2str(round(stats.FreqGo.CorrectRate*100)) '%']);
disp(['Go 2 = ' mat2str(round(stats.InfreqGo.CorrectRate*100)) '%']);
disp(['NoGo = ' mat2str(round(stats.NoGo.CorrectRate*100)) '%']);
disp(['False Starts = ' mat2str(round(stats.FalseStarts.Rate*100)) '%']);
save(out.matName, 'stats', '-append'); 

%% CLEAN UP
clear all
PsychPortAudio('Close');% Close the audio device. Without a second argument it closes all handles
CedrusResponseBox('CloseAll')

end


%% SUBFUNCTIONS

function [handle, startSoundTime, filePlayed] = playSound(audioFiles, neededFileType, waitToFinish, runStart) % Can put an additional argument that requests a specific file number of each filetype or could just pass the specific struct of interest rather than entire struct
 if ~isfield(audioFiles, neededFileType); error('Missing a needed audio file. Check Code.'); end; % make sure the file type we need has been added!
 %i = randi(length(audioFiles.(neededFileType))); % pick one of the exemplar files if a particular one has not been specified
 i=1;
 startSoundTime = GetSecs - runStart;
 PsychPortAudio('Start', audioFiles.(neededFileType)(i).pahandle, 1,0); % starts sound immediately
 switch waitToFinish
     case 1
        WaitSecs(audioFiles.(neededFileType)(i).duration); % waits for the whole duration of sound for it to play,if this wait is too short then sounds will be cutoff
        PsychPortAudio('Stop', audioFiles.(neededFileType)(i).pahandle); % Stop sound playback
     case 0
 end
 handle = audioFiles.(neededFileType)(i).pahandle;
 filePlayed = neededFileType;
end

function isGoPermuteArray = permuteWithMaxReps(anIsGoArray, numMaxReps, typesToLimit)
    %check the number of types and issue a warning
    tempPerm = randperm(length(anIsGoArray));
    types = struct;
    types(1).logical = ~ismember(anIsGoArray, typesToLimit);  % Make the first index the types that are not limited
    types(1).tempPerm = tempPerm(types(1).logical);
    % For each of the trial types, get their associated perms
    for type = 1:length(typesToLimit)
        types(type+1).logical = (anIsGoArray==typesToLimit(type));
        types(type+1).tempPerm = tempPerm(types(type+1).logical);
        types(type+1).sortedPerm = sort(types(type+1).tempPerm);
        
        flag = 1;
        while flag == 1
            %Find places where there are repetitions
            k = [true;diff(types(type+1).sortedPerm(:))~=1 ];
            s = cumsum(k);
            x =  histc(s,1:s(end));
            idx = find(k);
            out = types(type+1).sortedPerm(idx(x>numMaxReps));
            switch length(out) 
                case 0
                    flag = 0;
                otherwise
                    %Swap out the permutation indices
                    repeatToReplace = out+numMaxReps;
                    findReplacements = randperm(length(types(1).tempPerm));
                    nonRepeatToReplace = types(1).tempPerm(findReplacements(1:length(repeatToReplace)));
                    types(1).tempPerm(findReplacements(1:length(repeatToReplace))) = repeatToReplace;
                    types(type+1).tempPerm(ismember(types(type+1).tempPerm,repeatToReplace)) = nonRepeatToReplace;
                    types(type+1).sortedPerm = sort(types(type+1).tempPerm);
            end
        end     
    end
    isGoPermuteArray = [];
    for i = 1:length(types)
     isGoPermuteArray = [isGoPermuteArray, types(i).tempPerm];
    end
end

%% UNUSED FUNCTIONS
function responseMade = getResponse(responseMade)
[keyIsDown, ~, ~] = KbCheck(C.thisKeyboardIndex); % Can use the RT from here too
         if(keyIsDown && out.keyPress(currentTrial)==0) % If this is the first key press
            
         elseif(keyIsDown && out.keyPress(currentTrial)==1)
             
         end  
end

function out = getKeyResponse     
    [keyIsDown,seconds,keyCode] = KbCheck(-1); 
    if keyIsDown 
        out = 1; 
        response = find(keyCode); 
        if response == pauseKey % Pause
            cell2csv(['vwfaLoc_' SUBJ_ID '_run' num2str(RUN) '_quitAtTrial' iI],subjData); 
            save(['vwfaLoc_' SUBJ_ID '_run' num2str(RUN) '_quitAtTrial' iI],'subjData','RT_raw','responses');  
            while 1   
                [kId,sec,kC] = KbCheck(-1); 
                if kId 
                    resp = find(kC); 
                    if resp == continueKey % Continue 
                        out = 0;
                        return
                    elseif resp == quitKey 
                        cell2csv(['vwfaLoc_' SUBJ_ID '_run' num2str(RUN) '_quitAtTrial' iI],subjData); 
                        save(['vwfaLoc_' SUBJ_ID '_run' num2str(RUN) '_quitAtTrial' iI],'subjData','RT_raw','responses'); 
                        Screen('CloseAll');
                        fprintf('Experiment quit by pressing ESCAPE\n');
                        ShowCursor
                        out = 1;
                        break
                    end 
                end 
            end 
        elseif response == quitKey % Quit 
            cell2csv(['vwfaLoc_' SUBJ_ID '_run' num2str(RUN) '_quitAtTrial' iI],subjData); 
            save(['vwfaLoc_' SUBJ_ID '_run' num2str(RUN) '_quitAtTrial' iI],'subjData','RT_raw','responses'); 
            Screen('CloseAll');
            fprintf('Experiment quit by pressing ESCAPE\n');
            ShowCursor
            return
        end 
    else 
        out = 0; 
    end 
end 
                        
function [] = playTone (aFrequency)
samplingRate = 20500;
duration = 4; % seconds
values = 0:(1/samplingRate):duration;
amplitude = 10; 
y = amplitude*sin(2*pi*aFrequency*values);
sound(y, samplingRate); % Have to provide a sampling rate (i.e. X value) or it will use default and play the wrong tone
end