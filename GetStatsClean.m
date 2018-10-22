function [FreqGo, InfreqGo, NoGo, FalseStarts] = GetStatsClean (matfile)

DATA_DIR = 'Data/';
OVERWRITE = 0;

csv.Filename = strcat(DATA_DIR, 'AllSubjectsStats', '.csv');

if ~exist(csv.Filename, 'file')
   csv.out = {'Subject', 'TestOrPractice', 'RunNo', '%FreqGo', '%InFreqGo', '%NoGo', 'TotalTrials', 'hadFeedback', 'SoundDuration', 'SoundITI', 'FeedbackDuration', 'FeedbackITI', 'FreqGoCorrectRate', 'InfreqGoCorrectRate', 'NoGoCorrectRate', 'FreqGoRTsMean', 'InfreqGoRTsMean', 'NoGoRTsMean', 'FreqGoRTsSD', ...
            'InfreqGoRTsSD', 'NoGoRTsSD', 'FreqGoSound', 'InfreqGoSound', 'NoGoSound', 'AverageTrialDuration', 'AverageOnsetDuration', 'FalseStartRate', 'FSFreqGoRate', 'FSInfreqGoRate', 'FSNoGoRate', 'FSPriorFreqGoRate', 'FSInfreqGoRate', 'FSNoGoRate', 'FSRespondedTo'};
   cell2csv(csv.Filename,csv.out,',',0,0,'','');
end                      

matName = strcat(matfile);
load(matName);
%% Calculate Stats for Each Subject
FreqGo.isTrial = (out.trialNumber~=0)&(out.trialType==1);
FreqGo.Total = sum(FreqGo.isTrial);
FreqGo.isCorrect = (out.trialNumber~=0)&(out.trialType==1)&(out.keyPress==1);
FreqGo.CorrectNo = sum(FreqGo.isCorrect);
FreqGo.CorrectRate = FreqGo.CorrectNo/FreqGo.Total;
FreqGo.isIncorrect = (out.trialNumber~=0)&(out.trialType==1)&(out.keyPress==0);
FreqGo.IncorrectNo = sum(FreqGo.isIncorrect);
FreqGo.IncorrectRate = FreqGo.IncorrectNo/FreqGo.Total;
FreqGo.RTs = out.RT((out.trialNumber~=0)&(out.trialType==1)&(out.keyPress==1));
FreqGo.RTsMean = mean(FreqGo.RTs);
FreqGo.RTsSD = std(FreqGo.RTs);

InfreqGo.isTrial = (out.trialNumber~=0)&(out.trialType==-1);
InfreqGo.Total = sum(InfreqGo.isTrial);
InfreqGo.isCorrect = (out.trialNumber~=0)&(out.trialType==-1)&(out.keyPress==1);
InfreqGo.CorrectNo = sum(InfreqGo.isCorrect);
InfreqGo.CorrectRate = InfreqGo.CorrectNo/InfreqGo.Total;
InfreqGo.isIncorrect = (out.trialNumber~=0)&(out.trialType==-1)&(out.keyPress==0);
InfreqGo.IncorrectNo = sum(InfreqGo.isIncorrect);
InfreqGo.IncorrectRate = InfreqGo.IncorrectNo/InfreqGo.Total;
InfreqGo.RTs = out.RT((out.trialNumber~=0)&(out.trialType==-1)&(out.keyPress==1));
InfreqGo.RTsMean = mean(InfreqGo.RTs);
InfreqGo.RTsSD = std(InfreqGo.RTs);
%(InfreqGo.CorrectNo+InfreqGo.IncorrectNo)~=InfreqGo.Total && warning('Incorrect & Correct Trials Do Not Sum to Total Trials for Trial Type');
% Can also check that number of trials makes sense

NoGo.isTrial = (out.trialNumber~=0)&(out.trialType==0);
NoGo.Total = sum(NoGo.isTrial);
NoGo.isCorrect = (out.trialNumber~=0)&(out.trialType==0)&(out.keyPress==0);
NoGo.CorrectNo = sum(NoGo.isCorrect);
NoGo.CorrectRate = NoGo.CorrectNo/NoGo.Total;
NoGo.isIncorrect = (out.trialNumber~=0)&(out.trialType==0)&(out.keyPress==1);
NoGo.IncorrectNo = sum(NoGo.isIncorrect);
NoGo.IncorrectRate = NoGo.IncorrectNo/NoGo.Total;
NoGo.RTs = out.RT((out.trialNumber~=0)&(out.trialType==0)&(out.keyPress==1));
NoGo.RTsMean = mean(NoGo.RTs);
NoGo.RTsSD = std(NoGo.RTs);
%(NoGo.CorrectNo+NoGo.IncorrectNo)~=NoGo.Total && warning('Incorrect & Correct Trials Do Not Sum to Total Trials for Trial Type');(NoGo.CorrectNo+NoGo.Incorrect~=NoGo.Total) && warning('Incorrect & Correct Trials Do Not Sum to Total Trials for Trial Type');
% Can also check that number of trials makes sense
Trials.Total = FreqGo.Total + InfreqGo.Total + NoGo.Total; % this is according to attempted
FreqGo.Percentage = FreqGo.Total/Trials.Total;
InfreqGo.Percentage = InfreqGo.Total/Trials.Total;
NoGo.Percentage = NoGo.Total/Trials.Total;

AverageTrialDuration = mean(out.trialEnd(out.trialEnd~=0)-out.trialSoundStart(out.trialSoundStart~=0));

% Have to account for rest periods
startIdx = 0;
diffArray = [];
for j = 1:length(C.restBeforeTrials)
    if C.restBeforeTrials(j) ~=1
        diffArray = [diffArray diff(out.trialSoundStart(startIdx+1:C.restBeforeTrials(j)-1))];
        startIdx = C.restBeforeTrials(j);
    end
end
diffArray = [diffArray diff(out.trialSoundStart(C.restBeforeTrials(j):end))];
AverageOnsetDuration = mean(diffArray);
clear j startIdx diffArray

FalseStarts.Number = length(out.keyPressFalseStarts) + length(out.keyPressFeedback);  % When feedback is given, you can false-start before the trial during the feedback wait period
FalseStarts.Rate = FalseStarts.Number/sum(out.trialNumber~=0);  % can't use Trials.Total b/c you want real total rather than attempted total
FalseStarts.isFS = out.trialFalseStart;
% These are the trials on which the false starts occur
FalseStarts.FreqGoRate = sum((out.trialNumber~=0)&(out.trialType==1)&(out.trialFalseStart==1))/sum((out.trialNumber~=0)&(out.trialType==1));
FalseStarts.InfreqGoRate = sum((out.trialNumber~=0)&(out.trialType==-1)&(out.trialFalseStart==1))/sum((out.trialNumber~=0)&(out.trialType==-1));
FalseStarts.NoGoRate = sum((out.trialNumber~=0)&(out.trialType==0)&(out.trialFalseStart==1))/sum((out.trialNumber~=0)&(out.trialType==0));
% These are the trials preceding the false starts
priorToFS = find(out.trialFalseStart==1) - 1;
priorToFS(priorToFS==0) = [];
priorToFSArray = zeros(1, size(out.trialFalseStart,2));
priorToFSArray(priorToFS) = 1;
FalseStarts.PriorFreqGoRate = sum((out.trialNumber~=0)&(out.trialType==1)&(priorToFSArray==1))/sum((out.trialNumber~=0)&(out.trialType==1));
FalseStarts.PriorInfreqGoRate = sum((out.trialNumber~=0)&(out.trialType==-1)&(priorToFSArray==1))/sum((out.trialNumber~=0)&(out.trialType==-1));
FalseStarts.PriorNoGoRate = sum((out.trialNumber~=0)&(out.trialType==0)&(priorToFSArray==1))/sum((out.trialNumber~=0)&(out.trialType==0));
% This gives you a sense if the participant responded subsequent to the
% false start (for the same trial); out.keyPress refers to if key was
% pressed after the false start period
FalseStarts.Responses = sum((out.keyPressOriginalNoFS==1)&(out.trialFalseStart==1))/sum(out.trialFalseStart==1);

csv.out = {[out.subName], out.testOrPractice, out.runNum, FreqGo.Percentage, InfreqGo.Percentage, NoGo.Percentage, Trials.Total, C.giveFeedback, C.soundDuration, C.soundITI, C.feedbackDuration, C.feedbackITI, FreqGo.CorrectRate, InfreqGo.CorrectRate, NoGo.CorrectRate, FreqGo.RTsMean, InfreqGo.RTsMean, NoGo.RTsMean, FreqGo.RTsSD, ...
    InfreqGo.RTsSD, NoGo.RTsSD, audioFiles.Go1.name, audioFiles.Go2.name, audioFiles.NoGo.name, AverageTrialDuration, AverageOnsetDuration, FalseStarts.Rate, FalseStarts.FreqGoRate, FalseStarts.InfreqGoRate, FalseStarts.NoGoRate, FalseStarts.PriorFreqGoRate, FalseStarts.PriorInfreqGoRate, FalseStarts.PriorNoGoRate,  FalseStarts.Responses};
clearvars -except i listing DATA_DIR csv OVERWRITE FreqGo InfreqGo NoGo FalseStarts;

cell2csv(csv.Filename,csv.out,',',0,1,'','');
end