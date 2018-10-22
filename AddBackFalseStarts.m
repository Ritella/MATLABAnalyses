function [] = AddBackFalseStarts ()

mats = dir('Data/GNGC_*/GNGC_*-gonogo_R0*');

for m = 1:length(mats)  
  % only get runs 1-3 and don't get mats that are already cleaned   
  if ~isempty(regexp(mats(m).name, '_R0[123]', 'match')) && isempty(regexp(mats(m).name, '-clean', 'match'))
      
    load([mats(m).folder '/' mats(m).name]);
    
    % Fix outdated naming conventions for first 2 subjects
    if ~isempty(regexp(mats(m).name, 'GNGC_S_0[12]', 'match'))        
       out = output;
       out.runNum = out.runNo-4; % I initially made the first run = 4, since 3 practice runs. Renamed the files, but needs changing in the mat.
       C.runType = 'R';
       out.runName =[C.runType, sprintf('%02d', out.runNum)];
       out.matName = regexprep(out.matName, 'R[1-8]', ['R0' num2str(out.runNum)]);
       out.matName = regexprep(out.matName, '_S0', '_S_0');
       out.matName = regexprep(out.matName, 'GoNoGo_', '');
       out.matName = regexprep(out.matName, '', '_S_0');
       out.matName = regexprep(out.matName, '_R', '-gonogo_R');
       out.subName = [C.EXPT_STR '_' C.SUB_PRE '_' sprintf('%02d',out.sub)];
       C.restBeforeTrials = C.REST_BEFORE_TRIALS;
       C.restWait = C.REST_WAIT;
       C.soundITI = C.SOUND_ITI;
       C.feedbackDuration = C.FEEDBACK_DURATION;
       C.feedbackITI = C.FEEDBACK_ITI;
       C.giveFeedback = C.GIVE_FEEDBACK;
       C.soundDuration = C.SOUND_DURATION;
       out.taskName = 'gonogo';
       out.testOrPractice = 'T';
       clear output 
       out = rmfield(out, 'runNo');  % can't just clear a field
       C = rmfield(C, {'REST_BEFORE_TRIALS', 'REST_WAIT', 'SOUND_ITI', 'FEEDBACK_DURATION', 'FEEDBACK_ITI', 'GIVE_FEEDBACK', 'SOUND_DURATION'});
    end
    
    % Fix accidental naming of task as gonono for 3rd and 4th participants
    if strcmp(out.taskName, 'gonono')
        out.taskName = 'gonogo';
        out.matName = regexprep(out.matName, 'gonono', 'gonogo');
    end
    % Remove second (end of run) sync; not an actual key press from participant
    out.keyPressInterim = out.keyPressInterim(1:end-1);
    out.trialFalseStart = zeros(1, size(out.trialNumber,2));

    % There should be no false starts during Test Runs because no Feedback
    trialNumPrior = 0;  % In some cases, a participant can press more than once during a FS
    toRemove=[];
    for i = 1:length(out.keyPressFalseStarts)
       trialNum = out.trialNumber((out.trialSoundStart < out.keyPressFalseStarts(i))&(out.trialEnd > out.keyPressFalseStarts(i)));
       out.trialFalseStart(trialNum) = 1;
       if trialNum == trialNumPrior
           out.keyPressTrialExtras = out.keyPressFalseStarts(i);
           toRemove(end+1)=i;
       end
       trialNumPrior=trialNum;
    end
    out.keyPressFalseStarts(toRemove) = [];
    out.trialFalseStart = logical(out.trialFalseStart);
    % if participant pressed again for the same trial, count as
    % a double press. 
    % Note that you cannot have the following code: 
    % FSExtraPress = out.keyPressOnset(out.keyPressOnset(out.trialFalseStart)~=0);
    % The reason is because that code indexes into out.keyPressOnset (using
    % the logical out.trialFalseStart). Once you index into 
    FSExtraPress = out.keyPressOnset(((out.keyPressOnset .* out.trialFalseStart)~=0));
    out.keyPressTrialExtras = sort([FSExtraPress out.keyPressTrialExtras]);
    clear i trialNum FSExtraPress toRemove trialNumPrior
    % now add back the false starts
    out.keyPressOriginalNoFS=out.keyPress;
    %fprintf(out.matName);
    out.keyPress(out.trialFalseStart) = 1;
    out.keyPressOnset(out.trialFalseStart) = out.keyPressFalseStarts; % this assumes that out.keyPressFalseStarts is in order
    out.RT(out.trialFalseStart) = out.keyPressOnset(out.trialFalseStart) - out.trialSoundStart(out.trialFalseStart);
    

    %% GET STATISTICS
    out.newMatName = regexprep(out.matName, '.mat', '-clean.mat');
    clear ans m 
    save(out.newMatName,'-regexp', '\<(?!mats\>)\w*');
    [cleanStats.FreqGo, cleanStats.InfreqGo, cleanStats.NoGo, cleanStats.FalseStarts] = GetStatsClean(out.newMatName);
    save(out.newMatName, 'cleanStats', '-append'); 
    
    clearvars -except m mats;
  end
 end
 
end