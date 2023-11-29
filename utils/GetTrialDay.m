function [trialDay, sdate] = GetTrialDay(sessiondir)
% compute the trial day since the implant day of participants
% 
% works for session directory or full path directory
%
% return:
% 
%     trialday - trial day since implant day for the participant
%     sdate    - session date
%
% written by Ewina Pun, 2022.02.13

% initiate all participants and their implant days
participantList = {'t11','t5'};
implantdayList = {{'2019.09.24'}
                {'2016.08.17'}};
% identify participant
participantID = [];
pi = 1;
while pi <= length(participantList)
    if any(regexp(lower(sessiondir), participantList{pi}))
        participantID = pi;
        break
    end
    pi = pi + 1;
end

% error if still unable to identify
if isempty(participantID)
    error('Unable to identify the participant')
end

% find session date
sdate = regexp(sessiondir,'\d\d\d\d.\d\d\.\d\d','match');
if isempty(sdate)
    error('Unable to find a session date')
end

% compute trial day
trialDay = datenum(sdate, 'yyyy.mm.dd') ...
         - datenum(implantdayList{participantID}, 'yyyy.mm.dd');
if trialDay < 0
    error('Trial day cannot be prior to implant day')
end
sdate = sdate{:};
end