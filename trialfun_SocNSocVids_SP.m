function [trl,event] = trialfun_SocNSocVids_SP(data, TEinfo)
% RH oct 21
% for Safe Passage

% defines the trial events for the social and non-social videos
% valid trials using task engine info

% output: 
% trl file with trl = [begsample endsample offset task];
% event variable containing all the events in the raw data

%% read the header (needed for the samping rate) and the events
hdr        = data.hdr;
event      = data.events;

% for the events of interest, find the sample numbers (these are integers)
% for the events of interest, find the trigger values (these are strings in the case of BrainVision)
EVsample   = [event.sample]';
EVvalue    = cell2mat({event.value}');
numEvents = length(EVvalue);

%% Selecting event specific for braintools

% Begin samples
    StimCodes = [10 11]; % social/ non-social EEG markers
    XOnset = zeros(numEvents,length(StimCodes));
    for ii = 1:length(StimCodes)
        XOnset((StimCodes(1,ii) == EVvalue),ii) = 1;
    end
    XOnset(end,1) = 0; %set last XOnset to 0 (impossible to begin at last event)
    XOnset = sum(XOnset,2); IndTStim = find(XOnset==1);

% End samples
    StimOffset = 19; % social/ non-social
    XOffset = zeros(numEvents,1);
    XOffset((StimOffset == EVvalue),1) = 1; 
    XOffset(1,1) = 0; %set first XOffset to 0 (impossible to end at first event)
    IndTOff = find(XOffset==1); % select all events with matching stim codes

% define the begin, and endsample
    begsample = EVsample(IndTStim);
    endsample = EVsample(IndTStim+1) - 1;

    if isequal(size(IndTStim,1),size(IndTOff,1))
        % create offset and task variables
        offset = zeros(size(begsample,1),1);
        task = EVvalue(IndTStim,1);
        % create trl variable
        trl = [begsample endsample offset task];
        % test validity using task engine info
        if isequal(size(TEinfo,1),size(trl,1))
            validStim = ones(size(trl,1),1);
            for e = 1:size(validStim,1)
                    % sanity check EEG - TE duration
                    deltaEEG = (endsample(e) - begsample(e))/hdr.fs;
                    deltaTE = TEinfo(e,2);
                    if abs(deltaTE - deltaEEG) > 2
                        validStim(e,1) = 0;
                    end
                    % duration must be more than 10 seconds
                    if deltaEEG < 10 
                        validStim(e,1) = 0;
                    end
                    % looking at least 20%
                    if TEinfo(e,1) < .20
                        validStim(e,1) = 0;
                    end
            end
            % select valid trials
            IndValid = validStim==1;
            trl = trl(IndValid,:);
        else
            error('Uneven number of events in EEG and TE')
        end

    else 
        warning ('Unable to find matching begin and end samples')
    end

end % function