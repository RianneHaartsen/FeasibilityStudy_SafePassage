%% Safe Passage; 4) Extraction of EEG metrics for statistical analyses in R

% This script extracts the EEG metrics from the fieldtrip data and the 
% FOOOF results for statistical analyses in R. 

% SP list of EEG Measures 
% For each condition
% •	N trials
% For each region (frontal, central, parietal, occipital) * condition (social, non-social):
% •	Data amount ch x trl time series
% •	Canonical theta log power (averaged across 4-7Hz)
% •	Canonical alpha log power (averaged across 8-12Hz)
% •	FOOOF model fit: R2 
% •	FOOOF model fit: fit of the error
% •	FOOOF model: 1/f intercept
% •	FOOOF model: 1/f slope
% •	FOOOF model: Peak 1 frequency (Hz)
% •	FOOOF model: Peak 1 amplitude (log power)
% •	FOOOF model: Peak 2 frequency (Hz)
% •	FOOOF model: Peak 2 amplitude (log power)
% •	FOOOF model: Peak 3 frequency (Hz)
% •	FOOOF model: Peak 3 amplitude (log power)

% Output:
% 1 table with variables for social condition: SP_Soc_EEGmeasures_perRegion
% 1 table with variables for non-social condition: SP_NSoc_EEGmeasures_perRegion

% Created by Rianne Haartsen, February 2022


%% Get data and put into tables

cd XXX
% load data
load XXX/SP_preprocReREFavg_Power.mat
FOOOFdata = readtable('XXX/SP_FOOOF_wholesample/Fit_1_30Hz_data.csv');

Ind_data = zeros(height(SP_preprocReREFavg_Power),1);
for ss = 1:height(SP_preprocReREFavg_Power)
    if ~isempty(SP_preprocReREFavg_Power.PreprocEEG{ss})
        Ind_data(ss) = 1;
    end
end
SP_tab = SP_preprocReREFavg_Power(Ind_data == 1,:);



    
for ss = 1:height(SP_tab)
        
    Subj = SP_tab.ID{ss};
    % load social trials
        path_part = extractAfter(SP_tab.PreprocEEG{ss},'SafePassage/SP');
        load(strcat('XXX/SP',path_part),'Powerb_soctrls')
    % find frequencies
        Freqs = Powerb_soctrls.Pow_ft_struc.freq;
        Chind_FreqTheta = find(Freqs == 4) : find(Freqs == 7);
        Chind_FreqAlpha = find(Freqs == 8) : find(Freqs == 12);
    % find indices for ROIs
        Chind_front = find(ismember(Powerb_soctrls.Pow_ft_struc.label,{'Fz','F3','F4'}));
        Chind_cent = find(ismember(Powerb_soctrls.Pow_ft_struc.label,{'Cz','C3','C4'}));
        Chind_pari = find(ismember(Powerb_soctrls.Pow_ft_struc.label,{'Pz','P3','P4'}));
        Chind_occip = find(ismember(Powerb_soctrls.Pow_ft_struc.label,{'PO7','Oz','PO8'}));
    
    % Social trials %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % data availability
        NewRow.ID = {Subj};
        NewRow.Npr_SVid = SP_tab.Npr_SVid(ss);
        NewRow.Npr_SEps = SP_tab.Npr_SEps(ss);
        NewRow.Neps_Sclean = SP_tab.Ntr2b_Sclean(ss);
        
         % add new data to tracker table and save 
            cd('XXX')
            if exist('SP_Soc_DataAvailability.mat','file') == 2
                load SP_Soc_DataAvailability.mat
            end
            TrackerNew = struct2table(NewRow);
            if ss == 1
                SP_Soc_DataAvailability= TrackerNew;
            else
                SP_Soc_DataAvailability(ss,:) = TrackerNew;
            end
            save('SP_Soc_DataAvailability.mat','SP_Soc_DataAvailability')
          % clear up
            clear SP_Soc_DataAvailability NewRow TrackerNew
            
    % For separate regions %%%%%%%%%%%%%%%%%%%%%%%%%%
     NewRow.ID = {Subj};    
            
    % Frontal data %%%%%%%%%%%%%%%%%%%%%%
    % a) data availability
        Data_1Hz = Powerb_soctrls.Pow_ft_struc.powspctrm(:,Chind_front,Powerb_soctrls.Pow_ft_struc.freq == 1);
        NewRow.Fr_PercCleanEpxCh = (numel(find(~isnan(Data_1Hz))))/(numel(Data_1Hz))*100;
        clear Data_1Hz
    % b) canonical power
        NewRow.Fr_CanLogPow_theta = mean(Powerb_soctrls.Logpowspec.Front(1,Chind_FreqTheta),2);
        NewRow.Fr_CanLogPow_alpha = mean(Powerb_soctrls.Logpowspec.Front(1,Chind_FreqAlpha),2);
    % c) FOOOF parameters and AUC
        % find the row for subj front Soc
        Ind_FOOOF = zeros(height(FOOOFdata),1);
        SpectrumName = strcat('/',Subj,'_Ab_S_Fr');
        for ii = 1:height(FOOOFdata)
            if strcmp(FOOOFdata.Var1{ii},SpectrumName)
                Ind_FOOOF(ii) = 1;
            end
        end
        clear ii
        Ind_currSubj = find(Ind_FOOOF == 1);
        if ~isempty(Ind_currSubj)
            % FOOOF output
            NewRow.Fr_rSq = FOOOFdata.rSquared(Ind_currSubj);
            NewRow.Fr_fitError = FOOOFdata.fitError(Ind_currSubj);
            NewRow.Fr_intercept = FOOOFdata.intercept(Ind_currSubj);
            NewRow.Fr_slope = FOOOFdata.slope(Ind_currSubj);
            NewRow.Fr_peak1Freq = FOOOFdata.peak1Freq(Ind_currSubj);
            NewRow.Fr_peak1Ampl = FOOOFdata.peak1Amplitude(Ind_currSubj);
            NewRow.Fr_peak2Freq = FOOOFdata.peak2Freq(Ind_currSubj);
            NewRow.Fr_peak2Ampl = FOOOFdata.peak2Amplitude(Ind_currSubj);
            NewRow.Fr_peak3Freq = FOOOFdata.peak3Freq(Ind_currSubj);
            NewRow.Fr_peak3Ampl = FOOOFdata.peak3Amplitude(Ind_currSubj);
            % AUC
            % Predicted power
                offset = FOOOFdata.intercept(Ind_currSubj);
                exponent = FOOOFdata.slope(Ind_currSubj);
                freqs_cur = 1:.5:30;
                Predicted = zeros(length(freqs_cur),1);
                for ii = 1:length(Predicted)
                    freq_cur = freqs_cur(1,ii);
                    Predicted(ii,1) = offset - log10(0+freq_cur^(exponent));
                end
                clear ii
            % Observed power
                Freq_lb = find(Powerb_soctrls.Pow_ft_struc.freq == freqs_cur(1,1));
                Freq_ub = find(Powerb_soctrls.Pow_ft_struc.freq == freqs_cur(1,end));
                PowLog10 = log10(Powerb_soctrls.Abspowspec.Front(1,Freq_lb:Freq_ub));
            % For theta
                Freq_th_lb = find(freqs_cur == 4);
                Freq_th_ub = find(freqs_cur == 7);
                Predicted = Predicted';
                NewRow.Fr_AUCth_predFOOF = trapz(Predicted(1,Freq_th_lb:Freq_th_ub));
                NewRow.Fr_AUCth_obs = trapz(PowLog10(1,Freq_th_lb:Freq_th_ub));
                NewRow.Fr_AUCth_diff = NewRow.Fr_AUCth_obs - NewRow.Fr_AUCth_predFOOF;
            % For alpha
                Freq_al_lb = find(freqs_cur == 8);
                Freq_al_ub = find(freqs_cur == 12);
                NewRow.Fr_AUCal_predFOOF = trapz(Predicted(1,Freq_al_lb:Freq_al_ub));
                NewRow.Fr_AUCal_obs = trapz(PowLog10(1,Freq_al_lb:Freq_al_ub));
                NewRow.Fr_AUCal_diff = NewRow.Fr_AUCal_obs - NewRow.Fr_AUCal_predFOOF;
           clear offset exponent Predicted PowLog10 Freq_lb Freq_ub Freq_th_lb Freq_th_ub 
           clear Freq_al_lb Freq_al_ub freq_cur
        else
            warning('Unable to find spectrum for current subject')
            NewRow.Fr_rSq = NaN;
            NewRow.Fr_fitError = NaN;
            NewRow.Fr_intercept = NaN;
            NewRow.Fr_slope = NaN;
            NewRow.Fr_peak1Freq = NaN;
            NewRow.Fr_peak1Ampl = NaN;
            NewRow.Fr_peak2Freq = NaN;
            NewRow.Fr_peak2Ampl = NaN;
            NewRow.Fr_peak3Freq = NaN;
            NewRow.Fr_peak3Ampl = NaN;
            NewRow.Fr_AUCth_predFOOF = NaN;
            NewRow.Fr_AUCth_obs = NaN;
            NewRow.Fr_AUCth_diff = NaN;
            NewRow.Fr_AUCal_predFOOF = NaN;
            NewRow.Fr_AUCal_obs = NaN;
            NewRow.Fr_AUCal_diff = NaN;    
        end   
        clear Ind_FOOOF Ind_currSubj SpectrumName

        
    % Central data %%%%%%%%%%%%%%%%%%%%%%
    % a) data availability
        Data_1Hz = Powerb_soctrls.Pow_ft_struc.powspctrm(:,Chind_cent,Powerb_soctrls.Pow_ft_struc.freq == 1);
        NewRow.Ce_PercCleanEpxCh = (numel(find(~isnan(Data_1Hz))))/(numel(Data_1Hz))*100;
        clear Data_1Hz
    % b) canonical power
        NewRow.Ce_CanLogPow_theta = mean(Powerb_soctrls.Logpowspec.Cent(1,Chind_FreqTheta),2);
        NewRow.Ce_CanLogPow_alpha = mean(Powerb_soctrls.Logpowspec.Cent(1,Chind_FreqAlpha),2);
    % c) FOOOF parameters and AUC
        % find the row for subj front Soc
        Ind_FOOOF = zeros(height(FOOOFdata),1);
        SpectrumName = strcat('/',Subj,'_Ab_S_Ce');
        for ii = 1:height(FOOOFdata)
            if strcmp(FOOOFdata.Var1{ii},SpectrumName)
                Ind_FOOOF(ii) = 1;
            end
        end
        Ind_currSubj = find(Ind_FOOOF == 1);
        if ~isempty(Ind_currSubj)
            % FOOOF output
            NewRow.Ce_rSq = FOOOFdata.rSquared(Ind_currSubj);
            NewRow.Ce_fitError = FOOOFdata.fitError(Ind_currSubj);
            NewRow.Ce_intercept = FOOOFdata.intercept(Ind_currSubj);
            NewRow.Ce_slope = FOOOFdata.slope(Ind_currSubj);
            NewRow.Ce_peak1Freq = FOOOFdata.peak1Freq(Ind_currSubj);
            NewRow.Ce_peak1Ampl = FOOOFdata.peak1Amplitude(Ind_currSubj);
            NewRow.Ce_peak2Freq = FOOOFdata.peak2Freq(Ind_currSubj);
            NewRow.Ce_peak2Ampl = FOOOFdata.peak2Amplitude(Ind_currSubj);
            NewRow.Ce_peak3Freq = FOOOFdata.peak3Freq(Ind_currSubj);
            NewRow.Ce_peak3Ampl = FOOOFdata.peak3Amplitude(Ind_currSubj);
            % AUC
            % Predicted power
                offset = FOOOFdata.intercept(Ind_currSubj);
                exponent = FOOOFdata.slope(Ind_currSubj);
                freqs_cur = 1:.5:30;
                Predicted = zeros(length(freqs_cur),1);
                for ii = 1:length(Predicted)
                    freq_cur = freqs_cur(1,ii);
                    Predicted(ii,1) = offset - log10(0+freq_cur^(exponent));
                end
                clear ii
            % Observed power
                Freq_lb = find(Powerb_soctrls.Pow_ft_struc.freq == freqs_cur(1,1));
                Freq_ub = find(Powerb_soctrls.Pow_ft_struc.freq == freqs_cur(1,end));
                PowLog10 = log10(Powerb_soctrls.Abspowspec.Cent(1,Freq_lb:Freq_ub));
            % For theta
                Freq_th_lb = find(freqs_cur == 4);
                Freq_th_ub = find(freqs_cur == 7);
                Predicted = Predicted';
                NewRow.Ce_AUCth_predFOOF = trapz(Predicted(1,Freq_th_lb:Freq_th_ub));
                NewRow.Ce_AUCth_obs = trapz(PowLog10(1,Freq_th_lb:Freq_th_ub));
                NewRow.Ce_AUCth_diff = NewRow.Ce_AUCth_obs - NewRow.Ce_AUCth_predFOOF;
            % For alpha
                Freq_al_lb = find(freqs_cur == 8);
                Freq_al_ub = find(freqs_cur == 12);
                NewRow.Ce_AUCal_predFOOF = trapz(Predicted(1,Freq_al_lb:Freq_al_ub));
                NewRow.Ce_AUCal_obs = trapz(PowLog10(1,Freq_al_lb:Freq_al_ub));
                NewRow.Ce_AUCal_diff = NewRow.Ce_AUCal_obs - NewRow.Ce_AUCal_predFOOF;
           clear offset exponent Predicted PowLog10 Freq_lb Freq_ub Freq_th_lb Freq_th_ub 
           clear Freq_al_lb Freq_al_ub freq_cur
        else
            warning('Unable to find spectrum for current subject')
            NewRow.Ce_rSq = NaN;
            NewRow.Ce_fitError = NaN;
            NewRow.Ce_intercept = NaN;
            NewRow.Ce_slope = NaN;
            NewRow.Ce_peak1Freq = NaN;
            NewRow.Ce_peak1Ampl = NaN;
            NewRow.Ce_peak2Freq = NaN;
            NewRow.Ce_peak2Ampl = NaN;
            NewRow.Ce_peak3Freq = NaN;
            NewRow.Ce_peak3Ampl = NaN;
            NewRow.Ce_AUCth_predFOOF = NaN;
            NewRow.Ce_AUCth_obs = NaN;
            NewRow.Ce_AUCth_diff = NaN;
            NewRow.Ce_AUCal_predFOOF = NaN;
            NewRow.Ce_AUCal_obs = NaN;
            NewRow.Ce_AUCal_diff = NaN;    
        end
        clear Ind_FOOOF Ind_currSubj SpectrumName
    
    
    % Parietal data %%%%%%%%%%%%%%%%%%%%%%
    % a) data availability
        Data_1Hz = Powerb_soctrls.Pow_ft_struc.powspctrm(:,Chind_pari,Powerb_soctrls.Pow_ft_struc.freq == 1);
        NewRow.Pa_PercCleanEpxCh = (numel(find(~isnan(Data_1Hz))))/(numel(Data_1Hz))*100;
        clear Data_1Hz
    % b) canonical power
        NewRow.Pa_CanLogPow_theta = mean(Powerb_soctrls.Logpowspec.Pari(1,Chind_FreqTheta),2);
        NewRow.Pa_CanLogPow_alpha = mean(Powerb_soctrls.Logpowspec.Pari(1,Chind_FreqAlpha),2);
    % c) FOOOF parameters and AUC
        % find the row for subj front Soc
        Ind_FOOOF = zeros(height(FOOOFdata),1);
        SpectrumName = strcat('/',Subj,'_Ab_S_Pa');
        for ii = 1:height(FOOOFdata)
            if strcmp(FOOOFdata.Var1{ii},SpectrumName)
                Ind_FOOOF(ii) = 1;
            end
        end
        Ind_currSubj = find(Ind_FOOOF == 1);
        if ~isempty(Ind_currSubj)
            % FOOOF output
            NewRow.Pa_rSq = FOOOFdata.rSquared(Ind_currSubj);
            NewRow.Pa_fitError = FOOOFdata.fitError(Ind_currSubj);
            NewRow.Pa_intercept = FOOOFdata.intercept(Ind_currSubj);
            NewRow.Pa_slope = FOOOFdata.slope(Ind_currSubj);
            NewRow.Pa_peak1Freq = FOOOFdata.peak1Freq(Ind_currSubj);
            NewRow.Pa_peak1Ampl = FOOOFdata.peak1Amplitude(Ind_currSubj);
            NewRow.Pa_peak2Freq = FOOOFdata.peak2Freq(Ind_currSubj);
            NewRow.Pa_peak2Ampl = FOOOFdata.peak2Amplitude(Ind_currSubj);
            NewRow.Pa_peak3Freq = FOOOFdata.peak3Freq(Ind_currSubj);
            NewRow.Pa_peak3Ampl = FOOOFdata.peak3Amplitude(Ind_currSubj);
            % AUC
            % Predicted power
                offset = FOOOFdata.intercept(Ind_currSubj);
                exponent = FOOOFdata.slope(Ind_currSubj);
                freqs_cur = 1:.5:30;
                Predicted = zeros(length(freqs_cur),1);
                for ii = 1:length(Predicted)
                    freq_cur = freqs_cur(1,ii);
                    Predicted(ii,1) = offset - log10(0+freq_cur^(exponent));
                end
                clear ii
            % Observed power
                Freq_lb = find(Powerb_soctrls.Pow_ft_struc.freq == freqs_cur(1,1));
                Freq_ub = find(Powerb_soctrls.Pow_ft_struc.freq == freqs_cur(1,end));
                PowLog10 = log10(Powerb_soctrls.Abspowspec.Pari(1,Freq_lb:Freq_ub));
            % For theta
                Freq_th_lb = find(freqs_cur == 4);
                Freq_th_ub = find(freqs_cur == 7);
                Predicted = Predicted';
                NewRow.Pa_AUCth_predFOOF = trapz(Predicted(1,Freq_th_lb:Freq_th_ub));
                NewRow.Pa_AUCth_obs = trapz(PowLog10(1,Freq_th_lb:Freq_th_ub));
                NewRow.Pa_AUCth_diff = NewRow.Pa_AUCth_obs - NewRow.Pa_AUCth_predFOOF;
            % For alpha
                Freq_al_lb = find(freqs_cur == 8);
                Freq_al_ub = find(freqs_cur == 12);
                NewRow.Pa_AUCal_predFOOF = trapz(Predicted(1,Freq_al_lb:Freq_al_ub));
                NewRow.Pa_AUCal_obs = trapz(PowLog10(1,Freq_al_lb:Freq_al_ub));
                NewRow.Pa_AUCal_diff = NewRow.Pa_AUCal_obs - NewRow.Pa_AUCal_predFOOF;
           clear offset exponent Predicted PowLog10 Freq_lb Freq_ub Freq_th_lb Freq_th_ub 
           clear Freq_al_lb Freq_al_ub freq_cur
        else
            warning('Unable to find spectrum for current subject')
            NewRow.Pa_rSq = NaN;
            NewRow.Pa_fitError = NaN;
            NewRow.Pa_intercept = NaN;
            NewRow.Pa_slope = NaN;
            NewRow.Pa_peak1Freq = NaN;
            NewRow.Pa_peak1Ampl = NaN;
            NewRow.Pa_peak2Freq = NaN;
            NewRow.Pa_peak2Ampl = NaN;
            NewRow.Pa_peak3Freq = NaN;
            NewRow.Pa_peak3Ampl = NaN;
            NewRow.Pa_AUCth_predFOOF = NaN;
            NewRow.Pa_AUCth_obs = NaN;
            NewRow.Pa_AUCth_diff = NaN;
            NewRow.Pa_AUCal_predFOOF = NaN;
            NewRow.Pa_AUCal_obs = NaN;
            NewRow.Pa_AUCal_diff = NaN;    
        end
        clear Ind_FOOOF Ind_currSubj SpectrumName


    % Occipital data %%%%%%%%%%%%%%%%%%%%%%
    % a) data availability
        Data_1Hz = Powerb_soctrls.Pow_ft_struc.powspctrm(:,Chind_occip,Powerb_soctrls.Pow_ft_struc.freq == 1);
        NewRow.Oc_PercCleanEpxCh = (numel(find(~isnan(Data_1Hz))))/(numel(Data_1Hz))*100;
        clear Data_1Hz
    % b) canonical power
        NewRow.Oc_CanLogPow_theta = mean(Powerb_soctrls.Logpowspec.Occip(1,Chind_FreqTheta),2);
        NewRow.Oc_CanLogPow_alpha = mean(Powerb_soctrls.Logpowspec.Occip(1,Chind_FreqAlpha),2);
    % c) FOOOF parameters and AUC
        % find the row for subj front Soc
        Ind_FOOOF = zeros(height(FOOOFdata),1);
        SpectrumName = strcat('/',Subj,'_Ab_S_Oc');
        for ii = 1:height(FOOOFdata)
            if strcmp(FOOOFdata.Var1{ii},SpectrumName)
                Ind_FOOOF(ii) = 1;
            end
        end
        Ind_currSubj = find(Ind_FOOOF == 1);
        if ~isempty(Ind_currSubj)
            % FOOOF output
            NewRow.Oc_rSq = FOOOFdata.rSquared(Ind_currSubj);
            NewRow.Oc_fitError = FOOOFdata.fitError(Ind_currSubj);
            NewRow.Oc_intercept = FOOOFdata.intercept(Ind_currSubj);
            NewRow.Oc_slope = FOOOFdata.slope(Ind_currSubj);
            NewRow.Oc_peak1Freq = FOOOFdata.peak1Freq(Ind_currSubj);
            NewRow.Oc_peak1Ampl = FOOOFdata.peak1Amplitude(Ind_currSubj);
            NewRow.Oc_peak2Freq = FOOOFdata.peak2Freq(Ind_currSubj);
            NewRow.Oc_peak2Ampl = FOOOFdata.peak2Amplitude(Ind_currSubj);
            NewRow.Oc_peak3Freq = FOOOFdata.peak3Freq(Ind_currSubj);
            NewRow.Oc_peak3Ampl = FOOOFdata.peak3Amplitude(Ind_currSubj);
            % AUC
            % Predicted power
                offset = FOOOFdata.intercept(Ind_currSubj);
                exponent = FOOOFdata.slope(Ind_currSubj);
                freqs_cur = 1:.5:30;
                Predicted = zeros(length(freqs_cur),1);
                for ii = 1:length(Predicted)
                    freq_cur = freqs_cur(1,ii);
                    Predicted(ii,1) = offset - log10(0+freq_cur^(exponent));
                end
                clear ii
            % Observed power
                Freq_lb = find(Powerb_soctrls.Pow_ft_struc.freq == freqs_cur(1,1));
                Freq_ub = find(Powerb_soctrls.Pow_ft_struc.freq == freqs_cur(1,end));
                PowLog10 = log10(Powerb_soctrls.Abspowspec.Occip(1,Freq_lb:Freq_ub));
            % For theta
                Freq_th_lb = find(freqs_cur == 4);
                Freq_th_ub = find(freqs_cur == 7);
                Predicted = Predicted';
                NewRow.Oc_AUCth_predFOOF = trapz(Predicted(1,Freq_th_lb:Freq_th_ub));
                NewRow.Oc_AUCth_obs = trapz(PowLog10(1,Freq_th_lb:Freq_th_ub));
                NewRow.Oc_AUCth_diff = NewRow.Oc_AUCth_obs - NewRow.Oc_AUCth_predFOOF;
            % For alpha
                Freq_al_lb = find(freqs_cur == 8);
                Freq_al_ub = find(freqs_cur == 12);
                NewRow.Oc_AUCal_predFOOF = trapz(Predicted(1,Freq_al_lb:Freq_al_ub));
                NewRow.Oc_AUCal_obs = trapz(PowLog10(1,Freq_al_lb:Freq_al_ub));
                NewRow.Oc_AUCal_diff = NewRow.Oc_AUCal_obs - NewRow.Oc_AUCal_predFOOF;
           clear offset exponent Predicted PowLog10 Freq_lb Freq_ub Freq_th_lb Freq_th_ub 
           clear Freq_al_lb Freq_al_ub freq_cur
        else
            warning('Unable to find spectrum for current subject')
            NewRow.Oc_rSq = NaN;
            NewRow.Oc_fitError = NaN;
            NewRow.Oc_intercept = NaN;
            NewRow.Oc_slope = NaN;
            NewRow.Oc_peak1Freq = NaN;
            NewRow.Oc_peak1Ampl = NaN;
            NewRow.Oc_peak2Freq = NaN;
            NewRow.Oc_peak2Ampl = NaN;
            NewRow.Oc_peak3Freq = NaN;
            NewRow.Oc_peak3Ampl = NaN;
            NewRow.Oc_AUCth_predFOOF = NaN;
            NewRow.Oc_AUCth_obs = NaN;
            NewRow.Oc_AUCth_diff = NaN;
            NewRow.Oc_AUCal_predFOOF = NaN;
            NewRow.Oc_AUCal_obs = NaN;
            NewRow.Oc_AUCal_diff = NaN;    
        end
        clear Ind_FOOOF Ind_currSubj SpectrumName


        
         % add new data to tracker table and save 
            cd('XXX')
            if exist('SP_Soc_EEGmeasures_perRegion.mat','file') == 2
                load SP_Soc_EEGmeasures_perRegion.mat
            end
            TrackerNew = struct2table(NewRow);
            if ss == 1
                SP_Soc_EEGmeasures_perRegion= TrackerNew;
            else
                SP_Soc_EEGmeasures_perRegion(ss,:) = TrackerNew;
            end
            save('SP_Soc_EEGmeasures_perRegion.mat','SP_Soc_EEGmeasures_perRegion')
          % clear up
            clear SP_Soc_EEGmeasures_perRegion NewRow TrackerNew
end





%% Non-social trials



cd XXX
% load data
load XXX/SP_preprocReREFavg_Power.mat
FOOOFdata = readtable('XXX/SP_FOOOF_wholesample/Fit_1_30Hz_data.csv');

Ind_data = zeros(height(SP_preprocReREFavg_Power),1);
for ss = 1:height(SP_preprocReREFavg_Power)
    if ~isempty(SP_preprocReREFavg_Power.PreprocEEG{ss})
        Ind_data(ss) = 1;
    end
end
SP_tab = SP_preprocReREFavg_Power(Ind_data == 1,:);



    
for ss = 1:height(SP_tab)
        
    Subj = SP_tab.ID{ss};
    % load social trials
        path_part = extractAfter(SP_tab.PreprocEEG{ss},'SafePassage/SP');
        load(strcat('XXX/SP',path_part),'Powerb_nsoctrls')
    % find frequencies
        Freqs = Powerb_nsoctrls.Pow_ft_struc.freq;
        Chind_FreqTheta = find(Freqs == 4) : find(Freqs == 7);
        Chind_FreqAlpha = find(Freqs == 8) : find(Freqs == 12);
    % find indices for ROIs
        Chind_front = find(ismember(Powerb_nsoctrls.Pow_ft_struc.label,{'Fz','F3','F4'}));
        Chind_cent = find(ismember(Powerb_nsoctrls.Pow_ft_struc.label,{'Cz','C3','C4'}));
        Chind_pari = find(ismember(Powerb_nsoctrls.Pow_ft_struc.label,{'Pz','P3','P4'}));
        Chind_occip = find(ismember(Powerb_nsoctrls.Pow_ft_struc.label,{'PO7','Oz','PO8'}));
    
    % Non-social trials %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % data availability
        NewRow.ID = {Subj};
        NewRow.Npr_NSVid = SP_tab.Npr_NSVid(ss);
        NewRow.Npr_NSEps = SP_tab.Npr_NSEps(ss);
        NewRow.Neps_NSclean = SP_tab.Ntr2b_NSclean(ss);
        
         % add new data to tracker table and save 
            cd('XXX')
            if exist('SP_NSoc_DataAvailability.mat','file') == 2
                load SP_NSoc_DataAvailability.mat
            end
            TrackerNew = struct2table(NewRow);
            if ss == 1
                SP_NSoc_DataAvailability= TrackerNew;
            else
                SP_NSoc_DataAvailability(ss,:) = TrackerNew;
            end
            save('SP_NSoc_DataAvailability.mat','SP_NSoc_DataAvailability')
          % clear up
            clear SP_NSoc_DataAvailability NewRow TrackerNew
            
    % For separate regions %%%%%%%%%%%%%%%%%%%%%%%%%%
     NewRow.ID = {Subj};    
            
    % Frontal data %%%%%%%%%%%%%%%%%%%%%%
    % a) data availability
        Data_1Hz = Powerb_nsoctrls.Pow_ft_struc.powspctrm(:,Chind_front,Powerb_nsoctrls.Pow_ft_struc.freq == 1);
        NewRow.Fr_PercCleanEpxCh = (numel(find(~isnan(Data_1Hz))))/(numel(Data_1Hz))*100;
        clear Data_1Hz
    % b) canonical power
        NewRow.Fr_CanLogPow_theta = mean(Powerb_nsoctrls.Logpowspec.Front(1,Chind_FreqTheta),2);
        NewRow.Fr_CanLogPow_alpha = mean(Powerb_nsoctrls.Logpowspec.Front(1,Chind_FreqAlpha),2);
    % c) FOOOF parameters and AUC
        % find the row for subj front Soc
        Ind_FOOOF = zeros(height(FOOOFdata),1);
        SpectrumName = strcat('/',Subj,'_Ab_NS_Fr');
        for ii = 1:height(FOOOFdata)
            if strcmp(FOOOFdata.Var1{ii},SpectrumName)
                Ind_FOOOF(ii) = 1;
            end
        end
        clear ii
        Ind_currSubj = find(Ind_FOOOF == 1);
        if ~isempty(Ind_currSubj)
            % FOOOF output
            NewRow.Fr_rSq = FOOOFdata.rSquared(Ind_currSubj);
            NewRow.Fr_fitError = FOOOFdata.fitError(Ind_currSubj);
            NewRow.Fr_intercept = FOOOFdata.intercept(Ind_currSubj);
            NewRow.Fr_slope = FOOOFdata.slope(Ind_currSubj);
            NewRow.Fr_peak1Freq = FOOOFdata.peak1Freq(Ind_currSubj);
            NewRow.Fr_peak1Ampl = FOOOFdata.peak1Amplitude(Ind_currSubj);
            NewRow.Fr_peak2Freq = FOOOFdata.peak2Freq(Ind_currSubj);
            NewRow.Fr_peak2Ampl = FOOOFdata.peak2Amplitude(Ind_currSubj);
            NewRow.Fr_peak3Freq = FOOOFdata.peak3Freq(Ind_currSubj);
            NewRow.Fr_peak3Ampl = FOOOFdata.peak3Amplitude(Ind_currSubj);
            % AUC
            % Predicted power
                offset = FOOOFdata.intercept(Ind_currSubj);
                exponent = FOOOFdata.slope(Ind_currSubj);
                freqs_cur = 1:.5:30;
                Predicted = zeros(length(freqs_cur),1);
                for ii = 1:length(Predicted)
                    freq_cur = freqs_cur(1,ii);
                    Predicted(ii,1) = offset - log10(0+freq_cur^(exponent));
                end
                clear ii
            % Observed power
                Freq_lb = find(Powerb_nsoctrls.Pow_ft_struc.freq == freqs_cur(1,1));
                Freq_ub = find(Powerb_nsoctrls.Pow_ft_struc.freq == freqs_cur(1,end));
                PowLog10 = log10(Powerb_nsoctrls.Abspowspec.Front(1,Freq_lb:Freq_ub));
            % For theta
                Freq_th_lb = find(freqs_cur == 4);
                Freq_th_ub = find(freqs_cur == 7);
                Predicted = Predicted';
                NewRow.Fr_AUCth_predFOOF = trapz(Predicted(1,Freq_th_lb:Freq_th_ub));
                NewRow.Fr_AUCth_obs = trapz(PowLog10(1,Freq_th_lb:Freq_th_ub));
                NewRow.Fr_AUCth_diff = NewRow.Fr_AUCth_obs - NewRow.Fr_AUCth_predFOOF;
            % For alpha
                Freq_al_lb = find(freqs_cur == 8);
                Freq_al_ub = find(freqs_cur == 12);
                NewRow.Fr_AUCal_predFOOF = trapz(Predicted(1,Freq_al_lb:Freq_al_ub));
                NewRow.Fr_AUCal_obs = trapz(PowLog10(1,Freq_al_lb:Freq_al_ub));
                NewRow.Fr_AUCal_diff = NewRow.Fr_AUCal_obs - NewRow.Fr_AUCal_predFOOF;
           clear offset exponent Predicted PowLog10 Freq_lb Freq_ub Freq_th_lb Freq_th_ub 
           clear Freq_al_lb Freq_al_ub freq_cur
        else
            warning('Unable to find spectrum for current subject')
            NewRow.Fr_rSq = NaN;
            NewRow.Fr_fitError = NaN;
            NewRow.Fr_intercept = NaN;
            NewRow.Fr_slope = NaN;
            NewRow.Fr_peak1Freq = NaN;
            NewRow.Fr_peak1Ampl = NaN;
            NewRow.Fr_peak2Freq = NaN;
            NewRow.Fr_peak2Ampl = NaN;
            NewRow.Fr_peak3Freq = NaN;
            NewRow.Fr_peak3Ampl = NaN;
            NewRow.Fr_AUCth_predFOOF = NaN;
            NewRow.Fr_AUCth_obs = NaN;
            NewRow.Fr_AUCth_diff = NaN;
            NewRow.Fr_AUCal_predFOOF = NaN;
            NewRow.Fr_AUCal_obs = NaN;
            NewRow.Fr_AUCal_diff = NaN;    
        end   
        clear Ind_FOOOF Ind_currSubj SpectrumName

        
    % Central data %%%%%%%%%%%%%%%%%%%%%%
    % a) data availability
        Data_1Hz = Powerb_nsoctrls.Pow_ft_struc.powspctrm(:,Chind_cent,Powerb_nsoctrls.Pow_ft_struc.freq == 1);
        NewRow.Ce_PercCleanEpxCh = (numel(find(~isnan(Data_1Hz))))/(numel(Data_1Hz))*100;
        clear Data_1Hz
    % b) canonical power
        NewRow.Ce_CanLogPow_theta = mean(Powerb_nsoctrls.Logpowspec.Cent(1,Chind_FreqTheta),2);
        NewRow.Ce_CanLogPow_alpha = mean(Powerb_nsoctrls.Logpowspec.Cent(1,Chind_FreqAlpha),2);
    % c) FOOOF parameters and AUC
        % find the row for subj front Soc
        Ind_FOOOF = zeros(height(FOOOFdata),1);
        SpectrumName = strcat('/',Subj,'_Ab_NS_Ce');
        for ii = 1:height(FOOOFdata)
            if strcmp(FOOOFdata.Var1{ii},SpectrumName)
                Ind_FOOOF(ii) = 1;
            end
        end
        Ind_currSubj = find(Ind_FOOOF == 1);
        if ~isempty(Ind_currSubj)
            % FOOOF output
            NewRow.Ce_rSq = FOOOFdata.rSquared(Ind_currSubj);
            NewRow.Ce_fitError = FOOOFdata.fitError(Ind_currSubj);
            NewRow.Ce_intercept = FOOOFdata.intercept(Ind_currSubj);
            NewRow.Ce_slope = FOOOFdata.slope(Ind_currSubj);
            NewRow.Ce_peak1Freq = FOOOFdata.peak1Freq(Ind_currSubj);
            NewRow.Ce_peak1Ampl = FOOOFdata.peak1Amplitude(Ind_currSubj);
            NewRow.Ce_peak2Freq = FOOOFdata.peak2Freq(Ind_currSubj);
            NewRow.Ce_peak2Ampl = FOOOFdata.peak2Amplitude(Ind_currSubj);
            NewRow.Ce_peak3Freq = FOOOFdata.peak3Freq(Ind_currSubj);
            NewRow.Ce_peak3Ampl = FOOOFdata.peak3Amplitude(Ind_currSubj);
            % AUC
            % Predicted power
                offset = FOOOFdata.intercept(Ind_currSubj);
                exponent = FOOOFdata.slope(Ind_currSubj);
                freqs_cur = 1:.5:30;
                Predicted = zeros(length(freqs_cur),1);
                for ii = 1:length(Predicted)
                    freq_cur = freqs_cur(1,ii);
                    Predicted(ii,1) = offset - log10(0+freq_cur^(exponent));
                end
                clear ii
            % Observed power
                Freq_lb = find(Powerb_nsoctrls.Pow_ft_struc.freq == freqs_cur(1,1));
                Freq_ub = find(Powerb_nsoctrls.Pow_ft_struc.freq == freqs_cur(1,end));
                PowLog10 = log10(Powerb_nsoctrls.Abspowspec.Cent(1,Freq_lb:Freq_ub));
            % For theta
                Freq_th_lb = find(freqs_cur == 4);
                Freq_th_ub = find(freqs_cur == 7);
                Predicted = Predicted';
                NewRow.Ce_AUCth_predFOOF = trapz(Predicted(1,Freq_th_lb:Freq_th_ub));
                NewRow.Ce_AUCth_obs = trapz(PowLog10(1,Freq_th_lb:Freq_th_ub));
                NewRow.Ce_AUCth_diff = NewRow.Ce_AUCth_obs - NewRow.Ce_AUCth_predFOOF;
            % For alpha
                Freq_al_lb = find(freqs_cur == 8);
                Freq_al_ub = find(freqs_cur == 12);
                NewRow.Ce_AUCal_predFOOF = trapz(Predicted(1,Freq_al_lb:Freq_al_ub));
                NewRow.Ce_AUCal_obs = trapz(PowLog10(1,Freq_al_lb:Freq_al_ub));
                NewRow.Ce_AUCal_diff = NewRow.Ce_AUCal_obs - NewRow.Ce_AUCal_predFOOF;
           clear offset exponent Predicted PowLog10 Freq_lb Freq_ub Freq_th_lb Freq_th_ub 
           clear Freq_al_lb Freq_al_ub freq_cur
        else
            warning('Unable to find spectrum for current subject')
            NewRow.Ce_rSq = NaN;
            NewRow.Ce_fitError = NaN;
            NewRow.Ce_intercept = NaN;
            NewRow.Ce_slope = NaN;
            NewRow.Ce_peak1Freq = NaN;
            NewRow.Ce_peak1Ampl = NaN;
            NewRow.Ce_peak2Freq = NaN;
            NewRow.Ce_peak2Ampl = NaN;
            NewRow.Ce_peak3Freq = NaN;
            NewRow.Ce_peak3Ampl = NaN;
            NewRow.Ce_AUCth_predFOOF = NaN;
            NewRow.Ce_AUCth_obs = NaN;
            NewRow.Ce_AUCth_diff = NaN;
            NewRow.Ce_AUCal_predFOOF = NaN;
            NewRow.Ce_AUCal_obs = NaN;
            NewRow.Ce_AUCal_diff = NaN;    
        end
        clear Ind_FOOOF Ind_currSubj SpectrumName
    
    
    % Parietal data %%%%%%%%%%%%%%%%%%%%%%
    % a) data availability
        Data_1Hz = Powerb_nsoctrls.Pow_ft_struc.powspctrm(:,Chind_pari,Powerb_nsoctrls.Pow_ft_struc.freq == 1);
        NewRow.Pa_PercCleanEpxCh = (numel(find(~isnan(Data_1Hz))))/(numel(Data_1Hz))*100;
        clear Data_1Hz
    % b) canonical power
        NewRow.Pa_CanLogPow_theta = mean(Powerb_nsoctrls.Logpowspec.Pari(1,Chind_FreqTheta),2);
        NewRow.Pa_CanLogPow_alpha = mean(Powerb_nsoctrls.Logpowspec.Pari(1,Chind_FreqAlpha),2);
    % c) FOOOF parameters and AUC
        % find the row for subj front Soc
        Ind_FOOOF = zeros(height(FOOOFdata),1);
        SpectrumName = strcat('/',Subj,'_Ab_NS_Pa');
        for ii = 1:height(FOOOFdata)
            if strcmp(FOOOFdata.Var1{ii},SpectrumName)
                Ind_FOOOF(ii) = 1;
            end
        end
        Ind_currSubj = find(Ind_FOOOF == 1);
        if ~isempty(Ind_currSubj)
            % FOOOF output
            NewRow.Pa_rSq = FOOOFdata.rSquared(Ind_currSubj);
            NewRow.Pa_fitError = FOOOFdata.fitError(Ind_currSubj);
            NewRow.Pa_intercept = FOOOFdata.intercept(Ind_currSubj);
            NewRow.Pa_slope = FOOOFdata.slope(Ind_currSubj);
            NewRow.Pa_peak1Freq = FOOOFdata.peak1Freq(Ind_currSubj);
            NewRow.Pa_peak1Ampl = FOOOFdata.peak1Amplitude(Ind_currSubj);
            NewRow.Pa_peak2Freq = FOOOFdata.peak2Freq(Ind_currSubj);
            NewRow.Pa_peak2Ampl = FOOOFdata.peak2Amplitude(Ind_currSubj);
            NewRow.Pa_peak3Freq = FOOOFdata.peak3Freq(Ind_currSubj);
            NewRow.Pa_peak3Ampl = FOOOFdata.peak3Amplitude(Ind_currSubj);
            % AUC
            % Predicted power
                offset = FOOOFdata.intercept(Ind_currSubj);
                exponent = FOOOFdata.slope(Ind_currSubj);
                freqs_cur = 1:.5:30;
                Predicted = zeros(length(freqs_cur),1);
                for ii = 1:length(Predicted)
                    freq_cur = freqs_cur(1,ii);
                    Predicted(ii,1) = offset - log10(0+freq_cur^(exponent));
                end
                clear ii
            % Observed power
                Freq_lb = find(Powerb_nsoctrls.Pow_ft_struc.freq == freqs_cur(1,1));
                Freq_ub = find(Powerb_nsoctrls.Pow_ft_struc.freq == freqs_cur(1,end));
                PowLog10 = log10(Powerb_nsoctrls.Abspowspec.Pari(1,Freq_lb:Freq_ub));
            % For theta
                Freq_th_lb = find(freqs_cur == 4);
                Freq_th_ub = find(freqs_cur == 7);
                Predicted = Predicted';
                NewRow.Pa_AUCth_predFOOF = trapz(Predicted(1,Freq_th_lb:Freq_th_ub));
                NewRow.Pa_AUCth_obs = trapz(PowLog10(1,Freq_th_lb:Freq_th_ub));
                NewRow.Pa_AUCth_diff = NewRow.Pa_AUCth_obs - NewRow.Pa_AUCth_predFOOF;
            % For alpha
                Freq_al_lb = find(freqs_cur == 8);
                Freq_al_ub = find(freqs_cur == 12);
                NewRow.Pa_AUCal_predFOOF = trapz(Predicted(1,Freq_al_lb:Freq_al_ub));
                NewRow.Pa_AUCal_obs = trapz(PowLog10(1,Freq_al_lb:Freq_al_ub));
                NewRow.Pa_AUCal_diff = NewRow.Pa_AUCal_obs - NewRow.Pa_AUCal_predFOOF;
           clear offset exponent Predicted PowLog10 Freq_lb Freq_ub Freq_th_lb Freq_th_ub 
           clear Freq_al_lb Freq_al_ub freq_cur
        else
            warning('Unable to find spectrum for current subject')
            NewRow.Pa_rSq = NaN;
            NewRow.Pa_fitError = NaN;
            NewRow.Pa_intercept = NaN;
            NewRow.Pa_slope = NaN;
            NewRow.Pa_peak1Freq = NaN;
            NewRow.Pa_peak1Ampl = NaN;
            NewRow.Pa_peak2Freq = NaN;
            NewRow.Pa_peak2Ampl = NaN;
            NewRow.Pa_peak3Freq = NaN;
            NewRow.Pa_peak3Ampl = NaN;
            NewRow.Pa_AUCth_predFOOF = NaN;
            NewRow.Pa_AUCth_obs = NaN;
            NewRow.Pa_AUCth_diff = NaN;
            NewRow.Pa_AUCal_predFOOF = NaN;
            NewRow.Pa_AUCal_obs = NaN;
            NewRow.Pa_AUCal_diff = NaN;    
        end
        clear Ind_FOOOF Ind_currSubj SpectrumName


    % Occipital data %%%%%%%%%%%%%%%%%%%%%%
    % a) data availability
        Data_1Hz = Powerb_nsoctrls.Pow_ft_struc.powspctrm(:,Chind_occip,Powerb_nsoctrls.Pow_ft_struc.freq == 1);
        NewRow.Oc_PercCleanEpxCh = (numel(find(~isnan(Data_1Hz))))/(numel(Data_1Hz))*100;
        clear Data_1Hz
    % b) canonical power
        NewRow.Oc_CanLogPow_theta = mean(Powerb_nsoctrls.Logpowspec.Occip(1,Chind_FreqTheta),2);
        NewRow.Oc_CanLogPow_alpha = mean(Powerb_nsoctrls.Logpowspec.Occip(1,Chind_FreqAlpha),2);
    % c) FOOOF parameters and AUC
        % find the row for subj front Soc
        Ind_FOOOF = zeros(height(FOOOFdata),1);
        SpectrumName = strcat('/',Subj,'_Ab_NS_Oc');
        for ii = 1:height(FOOOFdata)
            if strcmp(FOOOFdata.Var1{ii},SpectrumName)
                Ind_FOOOF(ii) = 1;
            end
        end
        Ind_currSubj = find(Ind_FOOOF == 1);
        if ~isempty(Ind_currSubj)
            % FOOOF output
            NewRow.Oc_rSq = FOOOFdata.rSquared(Ind_currSubj);
            NewRow.Oc_fitError = FOOOFdata.fitError(Ind_currSubj);
            NewRow.Oc_intercept = FOOOFdata.intercept(Ind_currSubj);
            NewRow.Oc_slope = FOOOFdata.slope(Ind_currSubj);
            NewRow.Oc_peak1Freq = FOOOFdata.peak1Freq(Ind_currSubj);
            NewRow.Oc_peak1Ampl = FOOOFdata.peak1Amplitude(Ind_currSubj);
            NewRow.Oc_peak2Freq = FOOOFdata.peak2Freq(Ind_currSubj);
            NewRow.Oc_peak2Ampl = FOOOFdata.peak2Amplitude(Ind_currSubj);
            NewRow.Oc_peak3Freq = FOOOFdata.peak3Freq(Ind_currSubj);
            NewRow.Oc_peak3Ampl = FOOOFdata.peak3Amplitude(Ind_currSubj);
            % AUC
            % Predicted power
                offset = FOOOFdata.intercept(Ind_currSubj);
                exponent = FOOOFdata.slope(Ind_currSubj);
                freqs_cur = 1:.5:30;
                Predicted = zeros(length(freqs_cur),1);
                for ii = 1:length(Predicted)
                    freq_cur = freqs_cur(1,ii);
                    Predicted(ii,1) = offset - log10(0+freq_cur^(exponent));
                end
                clear ii
            % Observed power
                Freq_lb = find(Powerb_nsoctrls.Pow_ft_struc.freq == freqs_cur(1,1));
                Freq_ub = find(Powerb_nsoctrls.Pow_ft_struc.freq == freqs_cur(1,end));
                PowLog10 = log10(Powerb_nsoctrls.Abspowspec.Occip(1,Freq_lb:Freq_ub));
            % For theta
                Freq_th_lb = find(freqs_cur == 4);
                Freq_th_ub = find(freqs_cur == 7);
                Predicted = Predicted';
                NewRow.Oc_AUCth_predFOOF = trapz(Predicted(1,Freq_th_lb:Freq_th_ub));
                NewRow.Oc_AUCth_obs = trapz(PowLog10(1,Freq_th_lb:Freq_th_ub));
                NewRow.Oc_AUCth_diff = NewRow.Oc_AUCth_obs - NewRow.Oc_AUCth_predFOOF;
            % For alpha
                Freq_al_lb = find(freqs_cur == 8);
                Freq_al_ub = find(freqs_cur == 12);
                NewRow.Oc_AUCal_predFOOF = trapz(Predicted(1,Freq_al_lb:Freq_al_ub));
                NewRow.Oc_AUCal_obs = trapz(PowLog10(1,Freq_al_lb:Freq_al_ub));
                NewRow.Oc_AUCal_diff = NewRow.Oc_AUCal_obs - NewRow.Oc_AUCal_predFOOF;
           clear offset exponent Predicted PowLog10 Freq_lb Freq_ub Freq_th_lb Freq_th_ub 
           clear Freq_al_lb Freq_al_ub freq_cur
        else
            warning('Unable to find spectrum for current subject')
            NewRow.Oc_rSq = NaN;
            NewRow.Oc_fitError = NaN;
            NewRow.Oc_intercept = NaN;
            NewRow.Oc_slope = NaN;
            NewRow.Oc_peak1Freq = NaN;
            NewRow.Oc_peak1Ampl = NaN;
            NewRow.Oc_peak2Freq = NaN;
            NewRow.Oc_peak2Ampl = NaN;
            NewRow.Oc_peak3Freq = NaN;
            NewRow.Oc_peak3Ampl = NaN;
            NewRow.Oc_AUCth_predFOOF = NaN;
            NewRow.Oc_AUCth_obs = NaN;
            NewRow.Oc_AUCth_diff = NaN;
            NewRow.Oc_AUCal_predFOOF = NaN;
            NewRow.Oc_AUCal_obs = NaN;
            NewRow.Oc_AUCal_diff = NaN;    
        end
        clear Ind_FOOOF Ind_currSubj SpectrumName


        
         % add new data to tracker table and save 
            cd('XXX')
            if exist('SP_NSoc_EEGmeasures_perRegion.mat','file') == 2
                load SP_NSoc_EEGmeasures_perRegion.mat
            end
            TrackerNew = struct2table(NewRow);
            if ss == 1
                SP_NSoc_EEGmeasures_perRegion= TrackerNew;
            else
                SP_NSoc_EEGmeasures_perRegion(ss,:) = TrackerNew;
            end
            save('SP_NSoc_EEGmeasures_perRegion.mat','SP_NSoc_EEGmeasures_perRegion')
          % clear up
            clear SP_NSoc_EEGmeasures_perRegion NewRow TrackerNew
end