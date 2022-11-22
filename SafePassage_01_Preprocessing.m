%% Safe Passage: 1) preprocessing EEG data

% Preproc
% 1) Segment data for videos: from 10 soc / 11 nsoc to 19, \
%   a) check if duration is at least 10 sec but not more than 70sec
%   b) check if looking time is at least 20%
%   if both a and b, then count video segment as valid, otherwise not
% 2) Re-segment video data into 2 second epochs with 50% overlap
% 3) Filter band pas .1-48Hz, and detrend epochs
% 4) Automatic artefact identification: eogStat
% 5) exclude eogStat epochs
% 6) Identify scalp: flat channels, threshold [-150, 150] µV, jumps
% 7) Interpolate whole channel if >= 80% trials are bad (other 20% may be subthreshold but still bad)
% 8) Interpolate on a trl x ch basis
% 9) Identify residual artefacts
% 10) Calculate power for each channel for each trial
% 11) Reject trl x ch based on residual ARs: 
% - set ch x trl containing AR to NaN
% - set ch to NaN if fewer than 20 trials are clean

% 12) Calculate further derivatives: 
% 	- absolute power
% 	- log power
% 	- relative power
% 	- pow across frequencies of interest
% > frontal, & occipital 
% > for social, non-social, all trials

% Created by Rianne Haartsen, October 2021


%% Set up local paths to scripts
clear variables

%add eeglab path
addpath(genpath('XXX/eeglab14_1_2b')); 
%add fieldtrip path and set to defaults
addpath('XXX/MATLAB/fieldtrip-20180925'); 
ft_defaults
%add LM code TaskEngine2
addpath(genpath('XXX/TaskEngine2'))
addpath(genpath('XXX/lm_tools'));
addpath(genpath('XXX/eegtools'));
%add RH code
addpath('XXX/Safe_Passage');
% add SP path
addpath('XXX DATAPATH XXX/SafePassage')

%% Loop through all participants

cd('XXX DATAPATH XXX/SafePassage')
load XXXXX %ReliabilityTable: table with IDs and paths to folders with data
load SP_preproc_Power.mat

for ss = 1:height(ReliabilityTable)
    
    fprintf('Currently nr %i out of %i\n',ss,height(ReliabilityTable))

    Subj = ReliabilityTable.ID{ss}; %ppt code
    fprintf('Subject %s\n',Subj)
    
    % check if the subject has already been preprocessed
    Found = 0;
    for rr = 1:height(SP_preproc_Power_tab)
        if strcmp(SP_preproc_Power_tab.ID{rr},Subj)
            Found = 1;
        end
    end
        clear rr

if Found == 1 % data has already been saved for this subject continue to the next one
    disp('Data already read in and saved for subject:')
    disp(Subj)
    
else % preprocess the data 

    % 1) Read in data
    root_folder = 'XXX/SafePassage'; % folder with all ppt folders /Volumes/data_passage
    path_Fieldtrip = strcat(ReliabilityTable.Path_ses{ss},'/fieldtrip/');
    Session = ReliabilityTable.Session{ss};
    pathSession = ReliabilityTable.Path_ses{ss};

    % 2) Segment the data 
    
    if ~strcmp(Session, 'Test_Combined') % check if doesn data not come from combined enobio data
    
        % read in info from the session
        cd(path_Fieldtrip)
        load Session_FT.mat
        
        % check the channel layout
        if strcmp(WS_FT_data.label{1},'Ch1')  
            load('XXX/SP_20ch_layout_labels.mat')
            WS_FT_data.label = SP_20ch_layout_labels;
        end
        
        % get ET info about looking times and duration
        ses = teSession(pathSession);
        tab = teLogFilter(ses.Log.LogArray, 'task', 'restingvideos2', 'topic', 'trial_log_data');
        if isempty(tab)
          fprintf('Task not found in dataset.\n');
          TEinfo = NaN;
        else
          TEinfo = cell2mat([tab.et_looking_prop tab.movie_duration]);
        end
        
        % Data during soc and n-soc videos %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            [trl] = trialfun_SocNSocVids_SP(WS_FT_data, TEinfo);
            cfg = [];
            cfg.trl = trl;
            Video_seg = ft_redefinetrial(cfg, WS_FT_data);   
            clear cfg
        % Segment video data into 2-sec epochs
            cfg = [];
            cfg.length          = 2; % length is 2 second
            cfg.overlap         = .5; % overlap is 1 second (50% overlap)
            Epochs_2s = ft_redefinetrial(cfg, Video_seg);
            clear cfg
        % Numbers for presented videos and trials in total
            DataRetention.Npres_vid = size(trl,1);  
            DataRetention.NpresS_vid = size(find(trl(:,4)==10),1);  
            DataRetention.NpresNS_vid = size(find(trl(:,4)==11),1);  
            DataRetention.Npres_Eps = size(Epochs_2s.trialinfo,1); 
            DataRetention.NpresS_Eps = sum(Epochs_2s.trialinfo == 10,1); 
            DataRetention.NpresNS_Eps = sum(Epochs_2s.trialinfo == 11,1); 
        disp('Data have been segmented into 2-sec epochs with 50% overlap')
    else  % combined data
        
        pathEpochsdata = strcat(root_folder,'/',Subj,'/Test_Combined');
        load(strcat(pathEpochsdata,'/',Subj,'_SP_Epochs_2s.mat'))
        % Numbers for presented videos and trials in total
            DataRetention.Npres_vid = Npres_vid;  
            DataRetention.NpresS_vid = NpresS_vid;  
            DataRetention.NpresNS_vid = NpresNS_vid;  
            DataRetention.Npres_Eps = Npres_Eps; 
            DataRetention.NpresS_Eps = NpresS_Eps; 
            DataRetention.NpresNS_Eps = NpresNS_Eps; 
        clear Npres_vid NpresS_vid NpresNS_vid Npres_Eps NpresS_Eps NpresNS_Eps 

    end
    
    %% 3) Filter data
    % Filter the data and detrend
        cfg = [];
        cfg.bpfilter        = 'yes'; % bandpass filter
        cfg.bpfreq          = [.1, 48]; % bandpass freqs: .1 to 48 Hz 
        cfg.dftfilter       = 'yes'; % dft filter for line noise (shallow roll-off of bp filter)
        cfg.padding         = 3; % add 3 seconds of padding for filtering
        cfg.padtype         = 'mirror'; %mirror data for padding
        cfg.detrend         = 'yes'; % detrend the segments
        Epochs_2s_filt = ft_preprocessing(cfg, Epochs_2s);
    disp('Data have been filtered and detrended')    
    clear cfg trl
    
    %% 4) Exclude the O2/LS channel
        cfg = [];
        cfg.channel = {'all','-O2'};
        Epochs_2s_19ch = ft_selectdata(cfg, Epochs_2s_filt);
        clear cfg    
   
    %% 5) Prepare for artefact identification   
        
    % Prepare layout
    % layout
        cfg = [];
        cfg.layout = 'eeg1010.lay';
        cfg.channel = Epochs_2s_19ch.label;
        SP_layout1010Enobio = ft_prepare_layout(cfg);
        clear cfg
    % re-order channels in the order of the layout    
        Ch_ord_new = zeros(length(Epochs_2s_19ch.label),1);
        Layout_ord = SP_layout1010Enobio.label(1:19,1);
        for ch_lo = 1:length(Layout_ord)
            Ch_ord_new(ch_lo,1) = find(strcmp(Layout_ord{ch_lo},Epochs_2s_19ch.label)==1);
        end
        data1 = Epochs_2s_19ch;
        for tr = 1:length(data1.trial)
            data1.trial{tr} = Epochs_2s_19ch.trial{tr}(Ch_ord_new,:);
        end
        data1.label = Epochs_2s_19ch.label(Ch_ord_new);
        clear tr cfg ch_lo Ch_ord_new Layout_ord
        
        Data1_Epochs_2s_19ch = data1;
        if exist('tab','var')
            Data1_Epochs_2s_19ch.TEtable = table2struct(tab);
        end
        clear TEinfo WS_FT_data Video_seg tab Epochs_2s Epochs_2s_filt Epochs_2s_19ch
        clear data1
                  
    
    %% 6) Identify EOG artefacts and exclude trials
    
        % a) Identify eog and exclude trials
        % eogmethod LM; with 3-10Hz bp, blinklen = 50, zcrit = 2.5
            NotFrontChs = true(length(Data1_Epochs_2s_19ch.label),1);
            NotFrontChs(ismember(Data1_Epochs_2s_19ch.label,{'Fpz','AF7','AF8'}),1) = false;
            [Data1_AR1] = eegAR_Detect(Data1_Epochs_2s_19ch,'method','eogstat','maxsd',2.5,'excluded_channels',NotFrontChs);
        % select GOOD trials (no eog)
            cfg = [];
            cfg.trials = ~any(Data1_AR1.art,1);
            Data2_noeog = ft_selectdata(cfg, Data1_AR1);
            clear data dataAR1 NotFrontChs cfg

        DataRetention.Ntrls_1_noeog = size(Data2_noeog.trial,2);   
        DataRetention.NtrlsS_1_noeog = sum(Data2_noeog.trialinfo == 10,1);   
        DataRetention.NtrlsNS_1_noeog = sum(Data2_noeog.trialinfo == 11,1);   

    %% 7) Identify scalp artefacts and bad channels
        % b) Identify scalp artefacts; flat, threshold, jumps
        data_stripped = rmfieldIfPresent(Data2_noeog,...
                        {'art', 'art_type'});
        % i) flat
            Data2_AR2 = eegAR_Detect(data_stripped, 'method', 'flat'); 
        % ii) threshold
            Data2_AR2 = eegAR_Detect(Data2_AR2, 'method', 'minmax', 'threshold', [-150, 150]); 
        % iii) jumps: peak to peak exceeding, difference within 4ms
            Data2_AR2 = eegAR_Detect_RH(Data2_AR2, 'method', 'jumps', 'jumps_threshold', [200, 100]); 

        % c) find BAD channels with flat and/or exceeding thresholds/ jumps for 80% or
        % more of trials
            arts = any(Data2_AR2.art,3);
            PropBADtrls_x_ch = sum(arts~=0,2) / length(arts);
            BAD_ch = PropBADtrls_x_ch >= .8;
            if any(BAD_ch)
                artnew = zeros(size(arts));
                artnew(BAD_ch,:) = 1;
            else
                artnew = zeros(size(arts));
            end
            % add in BADch into art
            Data2_AR2.art = cat(3,Data2_AR2.art, artnew);
            Data2_AR2.art_type = cat(2,Data2_AR2.art_type, 'BADch');

            dataAR_cur = Data2_AR2;
            clear data_stripped arts artnew
            clear PropBADtrls_x_ch BAD_ch
         
    %% 8) Interpolate on a trl x channel basis 
    
    % neighbours template
            cfg = [];
            cfg.channel = dataAR_cur.label;
            cfg.method = 'distance';
            cfg.neighbourdist = .25;
            cfg.layout = SP_layout1010Enobio;
            nb_curr = ft_prepare_neighbours(cfg, dataAR_cur);   
            
            % flags to store which trials/channels were interpolated/excluded
            numTrials = length(dataAR_cur.trial);
            numChannels = length(dataAR_cur.label);
            interp = false(length(dataAR_cur.label), length(dataAR_cur.trial));
            cantInterp = false(length(dataAR_cur.label), length(dataAR_cur.trial));
            excl = false(length(dataAR_cur.label), length(dataAR_cur.trial));
            interpNeigh = cell(numChannels, numTrials);
            art = any(dataAR_cur.art, 3);

            % loop through trials
            tmp_trial = dataAR_cur.trial;
            tmp_canInterp = cell(numTrials, 1);
            tmp_bad = cell(numTrials, 1);
            parfor tr = 1:numTrials

                % check that there are some channels with artefacts on this current
                % trial
                if ~any(art(:, tr)), continue, end

                % select data from current trial        
                cfg = [];
                cfg.trials = false(numTrials, 1);
                cfg.trials(tr) = true;
                data_stripped = rmfieldIfPresent(dataAR_cur,...
                    {'interp', 'interpNeigh', 'art', 'chanExcl', 'art_type'});
                tmp = ft_selectdata(cfg, data_stripped);

                % extract channels with artefacts on this trial
                bad = art(:, tr);

                % find non-bad neighbours
                [canInterp, canInterpLabs, canInterpNb, canInterpSmry] =...
                    eegAR_FindInterpChans(dataAR_cur, bad, false, nb_curr);

                % store indices of channels that can't be interpolated
                cantInterp(:, tr) = bad & ~canInterp;

                if any(canInterp)
                    % interpolate
                    cfg = [];
                    cfg.method = 'spline';
                    cfg.badchannel = canInterpLabs;
                    cfg.neighbours = canInterpNb;
                    cfg.layout = SP_layout1010Enobio;
                    tmpi = ft_channelrepair(cfg, tmp);

                    % store interpolated data in temp structure
                    tmp_trial{tr} = tmpi.trial{1};            
                    tmp_bad{tr} = bad;
                    tmp_canInterp{tr} = canInterp;
                end

            end
            clear tr

            % update flags
            for tr = 1:numTrials
                dataAR_cur.trial{tr} = tmp_trial{tr};
                interpNeigh(tmp_canInterp{tr}, tr) = tmp_canInterp(tr);    
                interp(tmp_canInterp{tr}, tr) = true;
                excl(tmp_bad{tr} & ~tmp_canInterp{tr}, tr) = true;   
            end
            clear tr

            % summarise interpolation
            dataAR_cur.interp = interp;
            dataAR_cur.cantInterp = cantInterp;
            dataAR_cur.excl = excl;
            dataAR_cur.interp_summary.trialsIntPerChan = sum(interp, 2);                    % num channels with any trials interpolated
            dataAR_cur.interp_summary.chansIntPerTrial = sum(interp, 1);                   % num trials with any channels interpolated
            dataAR_cur.interp_summary.totalNumIntSegs = sum(interp(:));                     % total num of chan x trial interpolations
            dataAR_cur.interp_summary.propSegsInt = dataAR_cur.interp_summary.totalNumIntSegs / length(interp(:));     % prop of chan x trial interpolations
            dataAR_cur.interp_summary.intNeighbours = interpNeigh;
            
            Data3_interp = dataAR_cur;
            clear dataAR_cur interpNeigh interp excl cantInterp
            clear numTrials numChannels tmp_bad tmp_canInterp tmp_trial 
            clear art nb_curr cfg SP_layout1010Enobio
            
    %% 9) Identify residual artefacts/ bad channels that were not interpolated
    
            % run another AR detection to check how many ARs are left
            data_stripped3 = rmfieldIfPresent(Data3_interp,...
                    {'interp', 'interpNeigh', 'art', 'chanExcl', 'art_type','interp_summary'});
            % i) flat
                dataARres = eegAR_Detect(data_stripped3, 'method', 'flat'); 
            % ii) threshold
                dataARres = eegAR_Detect(dataARres, 'method', 'minmax', 'threshold', [-150, 150]); 
            % iii) jumps: peak to peak exceeding, difference within 4ms
                dataARres = eegAR_Detect_RH(dataARres, 'method', 'jumps', 'jumps_threshold', [200, 100]); 
            
            % find BAD channels:
            % i) with flat and/or exceeding thresholds/ jumps for 80% or more of trials
                arts2 = any(dataARres.art,3);
                artnew2 = zeros(size(arts2));
                PropBADtrls_x_ch = sum(arts2~=0,2) / length(arts2);
                BADi_ch2 = PropBADtrls_x_ch >= .8;
                if any(BADi_ch2) 
                    artnew2(BADi_ch2,:) = 1;
                end
             % ii) with fewer than 20 clean trials 
                NumBADtrls_x_ch = sum(arts2==0,2);
                BADii_ch2 = NumBADtrls_x_ch <= 20;
                if any(BADii_ch2) 
                    artnew2(BADii_ch2,:) = 1;
                end
                dataARres.art = cat(3,dataARres.art, artnew2);
                dataARres.art_type = cat(2,dataARres.art_type, 'BADch');
                
                Data3_ARres = dataARres;
                clear dataARres data_stripped3 arts2 artnew2 
                clear PropBADtrls_x_ch BADii_ch2 NumBADtrls_x_ch BADi_ch2

     %% 9b) re-ref on a trial basis to all clean channels

        numTrials = size(Data3_ARres.trial,2);
        art = any(Data3_ARres.art,3);
        data = Data3_ARres;
    
        % create temporary trials
                tmp_trial = cell(numTrials, 1);
                Rerefm_trial = cell(numTrials,1);

            % check in the number of trials is larger than 1
                if numTrials > 1
                    % loop through trials
                    for tr = 1:numTrials

%                         if ~any(art(:, tr,:)), continue, end

                        artCurrTrl = sum(squeeze(art(:, tr,:)),2);
                        % check if there more than 10 good channels
                        NgoodCh = size(find(artCurrTrl == 0),1);
                        if NgoodCh >= 15

                            % select data from current trial 
                            CurrTrlind = false(numTrials, 1);
                            CurrTrlind(tr) = true;

                            cfg = [];
                            cfg.trials = CurrTrlind;
                            chans = data.label;
                            cfg.channel = chans;
                            data_stripped = rmfieldIfPresent(data,...
                                {'art_type', 'art', 'summary'});
                            tmp = ft_selectdata(cfg, data_stripped);

                            cfg = [];
                            cfg.reref         = 'yes';
                            cfg.refchannel    = data.label(artCurrTrl == 0); % select clean channels only
                            cfg.refmethod     = 'avg';
                            tmpi = ft_preprocessing(cfg, tmp);

                            % store re-ref data in temp structure
                            tmp_trial{tr} = tmpi.trial{1};    

                            % save number of channels reref across
                            Rerefm_trial{tr,1} = strcat('all: ',num2str(NgoodCh));

                        else
                            disp('Not enough good channels for avg re-ref')
                            tmp_trial{tr} = data.trial{tr};
                            Rerefm_trial{tr,1} = 'unsuccessful';
                        end

                    end

                elseif numTrials == 1 % if there is only 1 trial

                    tr = 1;
                    artCurrTrl = sum(squeeze(art(:, tr,:)),2);
                        if ~any(artCurrTrl) % check if there any good channels

                        % check if there more than 10 good channels
                            NgoodCh = size(find(artCurrTrl == 0),1);
                            if NgoodCh >= 15

                                cfg = [];
                                chans = data.label;
                                cfg.channel = chans;
                                data_stripped = rmfieldIfPresent(data,...
                                    {'art_type', 'art', 'summary'});
                                tmp = ft_selectdata(cfg, data_stripped);

                                cfg = [];
                                cfg.reref         = 'yes';
                                cfg.refchannel    = data.label(artCurrTrl == 0);
                                cfg.refmethod     = 'avg';
                                tmpi = ft_preprocessing(cfg, tmp);

                                % store re-ref data in temp structure
                                tmp_trial{1} = tmpi.trial{1};   
                                Rerefm_trial{tr,1} = strcat('all: ',num2str(NgoodCh));
                            else
                                disp('Not enough good channels for avg re-ref')
                                tmp_trial{tr} = data.trial{tr};
                                Rerefm_trial{tr,1} = 'unsuccessful';
                            end
                        end      
                else
                    error('Unrecognised number of trials in dataset to re-ref to average')

                end

                % update flags and replace old trial with new re-referenced ones if
                % present
                for tr = 1:numTrials
                    if ~isempty(tmp_trial{tr}) % put in new trials
                        data.trial{tr} = tmp_trial{tr};
                    end
                end
            data.Reref_meth_trl = Rerefm_trial;
            Data3b_ReREF = data;
            clear tr tmp_trial tmpi tmp 
 
    %% 10) Power analysis across all trials and channels (to avoid error if substitute bad data with NaN)
        % strip data from LM fields
        data_stripped4 = rmfieldIfPresent(Data3b_ReREF,...
                    {'interp', 'interpNeigh', 'art', 'chanExcl', 'art_type','interp_summary'});
                
        % run power analysis        
        cfg                 = [];
        cfg.output          = 'pow';
        cfg.method          = 'mtmfft';
        cfg.taper           = 'hanning';
        cfg.foi             = 1:.5:48;  
        cfg.keeptrials      = 'yes'; 
        Power_trls = ft_freqanalysis(cfg, data_stripped4);
        Data4b_pow = Power_trls;
        clear cfg
        
    %% 11) Reject BAD data on a channel by trial basis 
        % set channel timeseries with artefacts to NaN 
        BAD_chxtrl = any(Data3_ARres.art, 3);
        % check if re-referencing was successful, exclude trial if
        % unsuccessful
        for tr = 1:numTrials
            if strcmp(Rerefm_trial{tr,1},'unsuccessful')
                BAD_chxtrl(:,tr) = 1;
            end
        end
        clear tr   

        % assign bad chxtrl NaN values
        for tr = 1:size(BAD_chxtrl,2)
            bad_chcurr = BAD_chxtrl(:,tr);
            Power_trls.powspctrm(tr, bad_chcurr,:) = NaN;
        end
        clear tr bad_chcurr   
        
        DataPowb_clean = Power_trls;
        clear Power_trls data_stripped4 BAD_chxtrl Rerefm_trial
        
        
    %% 12) Metrics of interest; power spectra for 4 RoIs and across all channels
    
    if ~isempty(DataPowb_clean)
        % find frequencies
            Freqs = DataPowb_clean.freq;
        % find indices for ROIs
            Chind_front = find(ismember(DataPowb_clean.label,{'Fz','F3','F4'}));
            Chind_cent = find(ismember(DataPowb_clean.label,{'Cz','C3','C4'}));
            Chind_pari = find(ismember(DataPowb_clean.label,{'Pz','P3','P4'}));
            Chind_occip = find(ismember(DataPowb_clean.label,{'PO7','Oz','PO8'}));
        
        % Across all trials
            Powerb_alltrls.abs = squeeze(nanmean(DataPowb_clean.powspctrm,1));
            Powerb_alltrls.log = log(Powerb_alltrls.abs);
                        
            Powerb_alltrls.Abspowspec.Front = mean(Powerb_alltrls.abs(Chind_front,:),1,'omitnan');
            Powerb_alltrls.Abspowspec.Cent = mean(Powerb_alltrls.abs(Chind_cent,:),1,'omitnan');
            Powerb_alltrls.Abspowspec.Pari = mean(Powerb_alltrls.abs(Chind_pari,:),1,'omitnan');
            Powerb_alltrls.Abspowspec.Occip = mean(Powerb_alltrls.abs(Chind_occip,:),1,'omitnan');
            Powerb_alltrls.Abspowspec.WholeScalp = mean(Powerb_alltrls.abs,1,'omitnan');
            
            Powerb_alltrls.Logpowspec.Front = mean(Powerb_alltrls.log(Chind_front,:),1,'omitnan');
            Powerb_alltrls.Logpowspec.Cent = mean(Powerb_alltrls.log(Chind_cent,:),1,'omitnan');
            Powerb_alltrls.Logpowspec.Pari = mean(Powerb_alltrls.log(Chind_pari,:),1,'omitnan');
            Powerb_alltrls.Logpowspec.Occip = mean(Powerb_alltrls.log(Chind_occip,:),1,'omitnan');
            Powerb_alltrls.Logpowspec.WholeScalp = mean(Powerb_alltrls.log,1,'omitnan');
         
            % check data retention
            ARfree = ~isnan(DataPowb_clean.powspctrm(:,:,1));
            GOODtrls = sum(ARfree,2) ~= 0;
            GOODtrls_ARfree = ARfree(GOODtrls,:);

            Powerb_alltrls.Ntrls = size(GOODtrls_ARfree,1);
            Powerb_alltrls.Ntrlxch_clean = sum(GOODtrls_ARfree,'all');
            Powerb_alltrls.Proptrlxch_clean = Powerb_alltrls.Ntrlxch_clean./numel(DataPowb_clean.powspctrm(:,:,1));
            
            DataRetention.Alltrlsb.Ntrls = Powerb_alltrls.Ntrls;
            DataRetention.Alltrlsb.Ntrlxch_clean = Powerb_alltrls.Ntrlxch_clean;
            DataRetention.Alltrlsb.Proptrlxch_clean = Powerb_alltrls.Proptrlxch_clean;
            clear ARfree Goodtrls GOODtrls_ARfree
            
            
            Powerb_alltrls.Pow_ft_struc = DataPowb_clean;

        % For social trials only
            % select social trials only
            cfg = [];
            cfg.trials = DataPowb_clean.trialinfo == 10;
            Powerb_soctrials = ft_selectdata(cfg,DataPowb_clean);
        
            Powerb_soctrls.abs = squeeze(nanmean(Powerb_soctrials.powspctrm,1));
            Powerb_soctrls.log = log(Powerb_soctrls.abs);
            
            Powerb_soctrls.Abspowspec.Front = mean(Powerb_soctrls.abs(Chind_front,:),1,'omitnan');
            Powerb_soctrls.Abspowspec.Cent = mean(Powerb_soctrls.abs(Chind_cent,:),1,'omitnan');
            Powerb_soctrls.Abspowspec.Pari = mean(Powerb_soctrls.abs(Chind_pari,:),1,'omitnan');
            Powerb_soctrls.Abspowspec.Occip = mean(Powerb_soctrls.abs(Chind_occip,:),1,'omitnan');
            Powerb_soctrls.Abspowspec.WholeScalp = mean(Powerb_soctrls.abs,1,'omitnan');
            
            Powerb_soctrls.Logpowspec.Front = mean(Powerb_soctrls.log(Chind_front,:),1,'omitnan');
            Powerb_soctrls.Logpowspec.Cent = mean(Powerb_soctrls.log(Chind_cent,:),1,'omitnan');
            Powerb_soctrls.Logpowspec.Pari = mean(Powerb_soctrls.log(Chind_pari,:),1,'omitnan');
            Powerb_soctrls.Logpowspec.Occip = mean(Powerb_soctrls.log(Chind_occip,:),1,'omitnan');
            Powerb_soctrls.Logpowspec.WholeScalp = mean(Powerb_soctrls.log,1,'omitnan');
            
            % check data retention
            ARfree = ~isnan(Powerb_soctrials.powspctrm(:,:,1));
            GOODtrls = sum(ARfree,2) ~= 0;
            GOODtrls_ARfree = ARfree(GOODtrls,:);

            Powerb_soctrls.Ntrls = size(GOODtrls_ARfree,1);
            Powerb_soctrls.Ntrlxch_clean = sum(GOODtrls_ARfree,'all');
            Powerb_soctrls.Proptrlxch_clean = Powerb_soctrls.Ntrlxch_clean./numel(Powerb_soctrials.powspctrm(:,:,1));
            
            DataRetention.soctrlsb.Ntrls = Powerb_soctrls.Ntrls;
            DataRetention.soctrlsb.Ntrlxch_clean = Powerb_soctrls.Ntrlxch_clean;
            DataRetention.soctrlsb.Proptrlxch_clean = Powerb_soctrls.Proptrlxch_clean;
            clear ARfree Goodtrls GOODtrls_ARfree
            
            Powerb_soctrls.Pow_ft_struc = Powerb_soctrials;
            clear Power_soctrials cfg

        % For non-social trials only
            % select social trials only
            cfg = [];
            cfg.trials = DataPowb_clean.trialinfo == 11;
            Powerb_nsoctrials = ft_selectdata(cfg,DataPowb_clean);
        
            Powerb_nsoctrls.abs = squeeze(nanmean(Powerb_nsoctrials.powspctrm,1));
            Powerb_nsoctrls.log = log(Powerb_nsoctrls.abs);
            
            Powerb_nsoctrls.Abspowspec.Front = mean(Powerb_nsoctrls.abs(Chind_front,:),1,'omitnan');
            Powerb_nsoctrls.Abspowspec.Cent = mean(Powerb_nsoctrls.abs(Chind_cent,:),1,'omitnan');
            Powerb_nsoctrls.Abspowspec.Pari = mean(Powerb_nsoctrls.abs(Chind_pari,:),1,'omitnan');
            Powerb_nsoctrls.Abspowspec.Occip = mean(Powerb_nsoctrls.abs(Chind_occip,:),1,'omitnan');
            Powerb_nsoctrls.Abspowspec.WholeScalp = mean(Powerb_nsoctrls.abs,1,'omitnan');
            
            Powerb_nsoctrls.Logpowspec.Front = mean(Powerb_nsoctrls.log(Chind_front,:),1,'omitnan');
            Powerb_nsoctrls.Logpowspec.Cent = mean(Powerb_nsoctrls.log(Chind_cent,:),1,'omitnan');
            Powerb_nsoctrls.Logpowspec.Pari = mean(Powerb_nsoctrls.log(Chind_pari,:),1,'omitnan');
            Powerb_nsoctrls.Logpowspec.Occip = mean(Powerb_nsoctrls.log(Chind_occip,:),1,'omitnan');
            Powerb_nsoctrls.Logpowspec.WholeScalp = mean(Powerb_nsoctrls.log,1,'omitnan');

            % check data retention
            ARfree = ~isnan(Powerb_nsoctrials.powspctrm(:,:,1));
            GOODtrls = sum(ARfree,2) ~= 0;
            GOODtrls_ARfree = ARfree(GOODtrls,:);

            Powerb_nsoctrls.Ntrls = size(GOODtrls_ARfree,1);
            Powerb_nsoctrls.Ntrlxch_clean = sum(GOODtrls_ARfree,'all');
            Powerb_nsoctrls.Proptrlxch_clean = Powerb_nsoctrls.Ntrlxch_clean./numel(Powerb_nsoctrials.powspctrm(:,:,1));
            
            DataRetention.nsoctrlsb.Ntrls = Powerb_nsoctrls.Ntrls;
            DataRetention.nsoctrlsb.Ntrlxch_clean = Powerb_nsoctrls.Ntrlxch_clean;
            DataRetention.nsoctrlsb.Proptrlxch_clean = Powerb_nsoctrls.Proptrlxch_clean;
            clear ARfree Goodtrls GOODtrls_ARfree

            Powerb_nsoctrls.Pow_ft_struc = Powerb_nsoctrials;
            clear Power_nsoctrials

            clear Chind_front Chind_cent Chind_pari Chind_occip cfg
        
    end
    
      

%% 13)  Save fieldtrip data in relevant folder  
        % save the data
        NameData = SP_preproc_Power_tab.PreprocEEG{ss};
        Name = extractBefore(NameData,'.mat');
        NewNameData = strcat(Name,'B_ReREF.mat');
        
        save(NewNameData,'Subj','DataRetention','Data3b_ReREF' ,'Data4b_pow',...
            'DataPowb_clean', 'Powerb_alltrls', 'Powerb_soctrls','Powerb_nsoctrls')
        
        disp('Fieldtrip EEG data saved')

    %% 14)   Save relevant data in tracker

       % Keep track of data
            NewRow.ID = {Subj};
            NewRow.PreprocEEG = {NewNameData};
            NewRow.Npr_Vid = DataRetention.Npres_vid;
            NewRow.Npr_SVid = DataRetention.NpresS_vid;
            NewRow.Npr_NSVid = DataRetention.NpresNS_vid;
            NewRow.Npr_Eps = DataRetention.Npres_Eps;
            NewRow.Npr_SEps = DataRetention.NpresS_Eps;
            NewRow.Npr_NSEps = DataRetention.NpresNS_Eps;
            NewRow.Ntr1_noEOG = DataRetention.Ntrls_1_noeog;
            NewRow.Ntr1_SnoEOG = DataRetention.NtrlsS_1_noeog;
            NewRow.Ntr1_NSnoEOG = DataRetention.NtrlsNS_1_noeog;
            NewRow.Ntr2b_clean = Powerb_alltrls.Ntrls;
            NewRow.Ntr2b_Sclean = Powerb_soctrls.Ntrls;
            NewRow.Ntr2b_NSclean = Powerb_nsoctrls.Ntrls;
            NewRow.Ntrch2b_clean = Powerb_alltrls.Ntrlxch_clean;
            NewRow.Ntrch2b_Sclean = Powerb_soctrls.Ntrlxch_clean;
            NewRow.Ntrch2b_NSclean = Powerb_nsoctrls.Ntrlxch_clean;
            NewRow.Proptrch2b_clean = Powerb_alltrls.Proptrlxch_clean;
            NewRow.Proptrch2b_Sclean = Powerb_soctrls.Proptrlxch_clean;
            NewRow.Proptrch2b_NSclean = Powerb_nsoctrls.Proptrlxch_clean;
            % power spectra
            % all trls
            NewRow.PowSp_AbsAll_Front = {Powerb_alltrls.Abspowspec.Front};
            NewRow.PowSp_AbsAll_Cent = {Powerb_alltrls.Abspowspec.Cent};
            NewRow.PowSp_AbsAll_Pari = {Powerb_alltrls.Abspowspec.Pari};
            NewRow.PowSp_AbsAll_Occip = {Powerb_alltrls.Abspowspec.Occip};
            NewRow.PowSp_AbsAll_WholeScalp = {Powerb_alltrls.Abspowspec.WholeScalp};
            NewRow.PowSp_LogAll_Front = {Powerb_alltrls.Logpowspec.Front};
            NewRow.PowSp_LogAll_Cent = {Powerb_alltrls.Logpowspec.Cent};
            NewRow.PowSp_LogAll_Pari = {Powerb_alltrls.Logpowspec.Pari};
            NewRow.PowSp_LogAll_Occip = {Powerb_alltrls.Logpowspec.Occip};
            NewRow.PowSp_LogAll_WholeScalp = {Powerb_alltrls.Logpowspec.WholeScalp};
            % soc trls
            NewRow.PowSp_AbsSoc_Front = {Powerb_soctrls.Abspowspec.Front};
            NewRow.PowSp_AbsSoc_Cent = {Powerb_soctrls.Abspowspec.Cent};
            NewRow.PowSp_AbsSoc_Pari = {Powerb_soctrls.Abspowspec.Pari};
            NewRow.PowSp_AbsSoc_Occip = {Powerb_soctrls.Abspowspec.Occip};
            NewRow.PowSp_AbsSoc_WholeScalp = {Powerb_soctrls.Abspowspec.WholeScalp};
            NewRow.PowSp_LogSoc_Front = {Powerb_soctrls.Logpowspec.Front};
            NewRow.PowSp_LogSoc_Cent = {Powerb_soctrls.Logpowspec.Cent};
            NewRow.PowSp_LogSoc_Pari = {Powerb_soctrls.Logpowspec.Pari};
            NewRow.PowSp_LogSoc_Occip = {Powerb_soctrls.Logpowspec.Occip};
            NewRow.PowSp_LogSoc_WholeScalp = {Powerb_soctrls.Logpowspec.WholeScalp};
            % nsoc trls
            NewRow.PowSp_AbsNSoc_Front = {Powerb_nsoctrls.Abspowspec.Front};
            NewRow.PowSp_AbsNSoc_Cent = {Powerb_nsoctrls.Abspowspec.Cent};
            NewRow.PowSp_AbsNSoc_Pari = {Powerb_nsoctrls.Abspowspec.Pari};
            NewRow.PowSp_AbsNSoc_Occip = {Powerb_nsoctrls.Abspowspec.Occip};
            NewRow.PowSp_AbsNSoc_WholeScalp = {Powerb_nsoctrls.Abspowspec.WholeScalp};
            NewRow.PowSp_LogNSoc_Front = {Powerb_nsoctrls.Logpowspec.Front};
            NewRow.PowSp_LogNSoc_Cent = {Powerb_nsoctrls.Logpowspec.Cent};
            NewRow.PowSp_LogNSoc_Pari = {Powerb_nsoctrls.Logpowspec.Pari};
            NewRow.PowSp_LogNSoc_Occip = {Powerb_nsoctrls.Logpowspec.Occip};
            NewRow.PowSp_LogNSoc_WholeScalp = {Powerb_nsoctrls.Logpowspec.WholeScalp};
            
            clear Powerb_alltrls Powerb_soctrls Powerb_nsoctrls DataRetention
            clear NewNameData Data3b_ReREF Data4b_pow 

        % add new data to tracker table and save 
            cd('XXX/SafePassage')
            TrackerNew = struct2table(NewRow);
                      
            SP_preprocReREFavg_Power(ss,:) = TrackerNew;
            save('SP_preprocReREFavg_Power.mat','SP_preprocReREFavg_Power')

        disp(Subj)
        disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
        disp('%%                Saving done                 %%')

        clear NewRow TrackerNew
        clear cfg d ext full_Sespath NameData 
        clear Subj Session pathSession ses 
        
    

end % subject already preprocessed

    clear Subj Found


    disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')    
    disp('%%          Ready for next subject            %%')    
    disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')

    
end



