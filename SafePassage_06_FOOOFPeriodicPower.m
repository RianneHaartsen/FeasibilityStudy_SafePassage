%% Safe Passage; 6) Extraction of periodic power for statistical analyses in R

% This script calculates the average power from the periodic component for the
% frequency band of interest, for each region, for each condition.
% This is done for theta and alpha band, for each region. 
% Theta = 4:7 Hz
% Alpha = 8:12 Hz


% Output:
% 1 table with variables for social condition: SP_Soc_FOOOFPeaks_perRegion
% 1 table with variables for non-social condition: SP_NSoc_FOOOFPeaks_perRegion

% Created by Rianne Haartsen, July 2022

addpath('XXX/Safe_Passage')

%% Social

cd('XXX')
load SP_Soc_EEGmeasures_perRegion.mat

cd XXX
% load data table
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
    % For separate regions %%%%%%%%%%%%%%%%%%%%%%%%%%
     NewRow.ID = {Subj};    
            
    % Frontal data %%%%%%%%%%%%%%%%%%%%%%
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
            NewRow.Fr_fitError = FOOOFdata.fitError(Ind_currSubj);
            NewRow.Fr_intercept = FOOOFdata.intercept(Ind_currSubj);
            NewRow.Fr_slope = FOOOFdata.slope(Ind_currSubj);
            % Aperiodic power
                offset = FOOOFdata.intercept(Ind_currSubj);
                exponent = FOOOFdata.slope(Ind_currSubj);
                freqs_cur = 1:.5:30;
                Aperiodic = zeros(length(freqs_cur),1);
                for ii = 1:length(Aperiodic)
                    freq_cur = freqs_cur(1,ii);
                    Aperiodic(ii,1) = offset - log10(0+freq_cur^(exponent));
                end
                clear ii
                Aperiodic = Aperiodic';
            % Observed power
                Freq_lb = find(Powerb_soctrls.Pow_ft_struc.freq == freqs_cur(1,1));
                Freq_ub = find(Powerb_soctrls.Pow_ft_struc.freq == freqs_cur(1,end));
                PowLog10 = log10(Powerb_soctrls.Abspowspec.Front(1,Freq_lb:Freq_ub));
            % For theta
                Freq_th_lb = find(freqs_cur == 4);
                Freq_th_ub = find(freqs_cur == 7);
                Obs_th_mn = mean(PowLog10(1,Freq_th_lb:Freq_th_ub),2);
                Aper_th_mn = mean(Aperiodic(1,Freq_th_lb:Freq_th_ub),2);
                NewRow.Fr_PerPower_th = Obs_th_mn - Aper_th_mn;
            % For alpha
                Freq_al_lb = find(freqs_cur == 8);
                Freq_al_ub = find(freqs_cur == 12);
                Obs_al_mn = mean(PowLog10(1,Freq_al_lb:Freq_al_ub),2);
                Aper_al_mn = mean(Aperiodic(1,Freq_al_lb:Freq_al_ub),2);
                NewRow.Fr_PerPower_al = Obs_al_mn - Aper_al_mn; 
           clear offset exponent Aperiodic PowLog10 Freq_lb Freq_ub Freq_th_lb Freq_th_ub 
           clear Freq_al_lb Freq_al_ub freq_cur
        else
            warning('Unable to find spectrum for current subject')
            NewRow.Fr_fitError = NaN;
            NewRow.Fr_intercept = NaN;
            NewRow.Fr_slope = NaN;
            NewRow.Fr_PerPower_th = NaN;
            NewRow.Fr_PerPower_al = NaN;    
        end   
        clear Ind_FOOOF Ind_currSubj SpectrumName 

    % Central data %%%%%%%%%%%%%%%%%%%%%
    % find the row for subj 
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
            NewRow.Ce_fitError = FOOOFdata.fitError(Ind_currSubj);
            NewRow.Ce_intercept = FOOOFdata.intercept(Ind_currSubj);
            NewRow.Ce_slope = FOOOFdata.slope(Ind_currSubj);
            % Aperiodic power
                offset = FOOOFdata.intercept(Ind_currSubj);
                exponent = FOOOFdata.slope(Ind_currSubj);
                freqs_cur = 1:.5:30;
                Aperiodic = zeros(length(freqs_cur),1);
                for ii = 1:length(Aperiodic)
                    freq_cur = freqs_cur(1,ii);
                    Aperiodic(ii,1) = offset - log10(0+freq_cur^(exponent));
                end
                clear ii
                Aperiodic = Aperiodic';
            % Observed power
                Freq_lb = find(Powerb_soctrls.Pow_ft_struc.freq == freqs_cur(1,1));
                Freq_ub = find(Powerb_soctrls.Pow_ft_struc.freq == freqs_cur(1,end));
                PowLog10 = log10(Powerb_soctrls.Abspowspec.Cent(1,Freq_lb:Freq_ub));
            % For theta
                Freq_th_lb = find(freqs_cur == 4);
                Freq_th_ub = find(freqs_cur == 7);
                Obs_th_mn = mean(PowLog10(1,Freq_th_lb:Freq_th_ub),2);
                Aper_th_mn = mean(Aperiodic(1,Freq_th_lb:Freq_th_ub),2);
                NewRow.Ce_PerPower_th = Obs_th_mn - Aper_th_mn;
            % For alpha
                Freq_al_lb = find(freqs_cur == 8);
                Freq_al_ub = find(freqs_cur == 12);
                Obs_al_mn = mean(PowLog10(1,Freq_al_lb:Freq_al_ub),2);
                Aper_al_mn = mean(Aperiodic(1,Freq_al_lb:Freq_al_ub),2);
                NewRow.Ce_PerPower_al = Obs_al_mn - Aper_al_mn; 
           clear offset exponent Aperiodic PowLog10 Freq_lb Freq_ub Freq_th_lb Freq_th_ub 
           clear Freq_al_lb Freq_al_ub freq_cur
        else
            warning('Unable to find spectrum for current subject')
            NewRow.Ce_fitError = NaN;
            NewRow.Ce_intercept = NaN;
            NewRow.Ce_slope = NaN;
            NewRow.Ce_PerPower_th = NaN;
            NewRow.Ce_PerPower_al = NaN;    
        end   
        clear Ind_FOOOF Ind_currSubj SpectrumName 
        
        
     % Parietal data %%%%%%%%%%%%%%%%%%%%%%
        % find the row for subj 
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
            NewRow.Pa_fitError = FOOOFdata.fitError(Ind_currSubj);
            NewRow.Pa_intercept = FOOOFdata.intercept(Ind_currSubj);
            NewRow.Pa_slope = FOOOFdata.slope(Ind_currSubj);
            % Aperiodic power
                offset = FOOOFdata.intercept(Ind_currSubj);
                exponent = FOOOFdata.slope(Ind_currSubj);
                freqs_cur = 1:.5:30;
                Aperiodic = zeros(length(freqs_cur),1);
                for ii = 1:length(Aperiodic)
                    freq_cur = freqs_cur(1,ii);
                    Aperiodic(ii,1) = offset - log10(0+freq_cur^(exponent));
                end
                clear ii
                Aperiodic = Aperiodic';
            % Observed power
                Freq_lb = find(Powerb_soctrls.Pow_ft_struc.freq == freqs_cur(1,1));
                Freq_ub = find(Powerb_soctrls.Pow_ft_struc.freq == freqs_cur(1,end));
                PowLog10 = log10(Powerb_soctrls.Abspowspec.Pari(1,Freq_lb:Freq_ub));
            % For theta
                Freq_th_lb = find(freqs_cur == 4);
                Freq_th_ub = find(freqs_cur == 7);
                Obs_th_mn = mean(PowLog10(1,Freq_th_lb:Freq_th_ub),2);
                Aper_th_mn = mean(Aperiodic(1,Freq_th_lb:Freq_th_ub),2);
                NewRow.Pa_PerPower_th = Obs_th_mn - Aper_th_mn;
            % For alpha
                Freq_al_lb = find(freqs_cur == 8);
                Freq_al_ub = find(freqs_cur == 12);
                Obs_al_mn = mean(PowLog10(1,Freq_al_lb:Freq_al_ub),2);
                Aper_al_mn = mean(Aperiodic(1,Freq_al_lb:Freq_al_ub),2);
                NewRow.Pa_PerPower_al = Obs_al_mn - Aper_al_mn; 
           clear offset exponent Aperiodic PowLog10 Freq_lb Freq_ub Freq_th_lb Freq_th_ub 
           clear Freq_al_lb Freq_al_ub freq_cur
        else
            warning('Unable to find spectrum for current subject')
            NewRow.Pa_fitError = NaN;
            NewRow.Pa_intercept = NaN;
            NewRow.Pa_slope = NaN;
            NewRow.Pa_PerPower_th = NaN;
            NewRow.Pa_PerPower_al = NaN;    
        end   
        clear Ind_FOOOF Ind_currSubj SpectrumName 
        
    % Occipital data %%%%%%%%%%%%%%%%%%%%%%
        % find the row for subj 
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
            NewRow.Oc_fitError = FOOOFdata.fitError(Ind_currSubj);
            NewRow.Oc_intercept = FOOOFdata.intercept(Ind_currSubj);
            NewRow.Oc_slope = FOOOFdata.slope(Ind_currSubj);
            % Aperiodic power
                offset = FOOOFdata.intercept(Ind_currSubj);
                exponent = FOOOFdata.slope(Ind_currSubj);
                freqs_cur = 1:.5:30;
                Aperiodic = zeros(length(freqs_cur),1);
                for ii = 1:length(Aperiodic)
                    freq_cur = freqs_cur(1,ii);
                    Aperiodic(ii,1) = offset - log10(0+freq_cur^(exponent));
                end
                clear ii
                Aperiodic = Aperiodic';
            % Observed power
                Freq_lb = find(Powerb_soctrls.Pow_ft_struc.freq == freqs_cur(1,1));
                Freq_ub = find(Powerb_soctrls.Pow_ft_struc.freq == freqs_cur(1,end));
                PowLog10 = log10(Powerb_soctrls.Abspowspec.Occip(1,Freq_lb:Freq_ub));
            % For theta
                Freq_th_lb = find(freqs_cur == 4);
                Freq_th_ub = find(freqs_cur == 7);
                Obs_th_mn = mean(PowLog10(1,Freq_th_lb:Freq_th_ub),2);
                Aper_th_mn = mean(Aperiodic(1,Freq_th_lb:Freq_th_ub),2);
                NewRow.Oc_PerPower_th = Obs_th_mn - Aper_th_mn;
            % For alpha
                Freq_al_lb = find(freqs_cur == 8);
                Freq_al_ub = find(freqs_cur == 12);
                Obs_al_mn = mean(PowLog10(1,Freq_al_lb:Freq_al_ub),2);
                Aper_al_mn = mean(Aperiodic(1,Freq_al_lb:Freq_al_ub),2);
                NewRow.Oc_PerPower_al = Obs_al_mn - Aper_al_mn; 
           clear offset exponent Aperiodic PowLog10 Freq_lb Freq_ub Freq_th_lb Freq_th_ub 
           clear Freq_al_lb Freq_al_ub freq_cur
        else
            warning('Unable to find spectrum for current subject')
            NewRow.Oc_fitError = NaN;
            NewRow.Oc_intercept = NaN;
            NewRow.Oc_slope = NaN;
            NewRow.Oc_PerPower_th = NaN;
            NewRow.Oc_PerPower_al = NaN;    
        end   
        clear Ind_FOOOF Ind_currSubj SpectrumName 
        
        
        % add new data to tracker table and save 
            cd('XXX')
            if exist('SP_Soc_PeriodicPower.mat','file') == 2
                load SP_Soc_PeriodicPower.mat
            end
            TrackerNew = struct2table(NewRow);
            if ss == 1
                SP_Soc_PeriodicPower = TrackerNew;
            else
                SP_Soc_PeriodicPower(ss,:) = TrackerNew;
            end
            save('SP_Soc_PeriodicPower.mat','SP_Soc_PeriodicPower')
          % clear up
            clear SP_Soc_PeriodicPower NewRow TrackerNew
        
        
end





%% Non-Social

cd XXX
% load data table
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
    
    % Social trials %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
    % For separate regions %%%%%%%%%%%%%%%%%%%%%%%%%%
     NewRow.ID = {Subj};    
            
    % Frontal data %%%%%%%%%%%%%%%%%%%%%%
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
            NewRow.Fr_fitError = FOOOFdata.fitError(Ind_currSubj);
            NewRow.Fr_intercept = FOOOFdata.intercept(Ind_currSubj);
            NewRow.Fr_slope = FOOOFdata.slope(Ind_currSubj);
            % Aperiodic power
                offset = FOOOFdata.intercept(Ind_currSubj);
                exponent = FOOOFdata.slope(Ind_currSubj);
                freqs_cur = 1:.5:30;
                Aperiodic = zeros(length(freqs_cur),1);
                for ii = 1:length(Aperiodic)
                    freq_cur = freqs_cur(1,ii);
                    Aperiodic(ii,1) = offset - log10(0+freq_cur^(exponent));
                end
                clear ii
                Aperiodic = Aperiodic';
            % Observed power
                Freq_lb = find(Powerb_nsoctrls.Pow_ft_struc.freq == freqs_cur(1,1));
                Freq_ub = find(Powerb_nsoctrls.Pow_ft_struc.freq == freqs_cur(1,end));
                PowLog10 = log10(Powerb_nsoctrls.Abspowspec.Front(1,Freq_lb:Freq_ub));
            % For theta
                Freq_th_lb = find(freqs_cur == 4);
                Freq_th_ub = find(freqs_cur == 7);
                Obs_th_mn = mean(PowLog10(1,Freq_th_lb:Freq_th_ub),2);
                Aper_th_mn = mean(Aperiodic(1,Freq_th_lb:Freq_th_ub),2);
                NewRow.Fr_PerPower_th = Obs_th_mn - Aper_th_mn;
            % For alpha
                Freq_al_lb = find(freqs_cur == 8);
                Freq_al_ub = find(freqs_cur == 12);
                Obs_al_mn = mean(PowLog10(1,Freq_al_lb:Freq_al_ub),2);
                Aper_al_mn = mean(Aperiodic(1,Freq_al_lb:Freq_al_ub),2);
                NewRow.Fr_PerPower_al = Obs_al_mn - Aper_al_mn; 
           clear offset exponent Aperiodic PowLog10 Freq_lb Freq_ub Freq_th_lb Freq_th_ub 
           clear Freq_al_lb Freq_al_ub freq_cur
        else
            warning('Unable to find spectrum for current subject')
            NewRow.Fr_fitError = NaN;
            NewRow.Fr_intercept = NaN;
            NewRow.Fr_slope = NaN;
            NewRow.Fr_PerPower_th = NaN;
            NewRow.Fr_PerPower_al = NaN;    
        end   
        clear Ind_FOOOF Ind_currSubj SpectrumName 

    % Central data %%%%%%%%%%%%%%%%%%%%%
    % find the row for subj 
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
            NewRow.Ce_fitError = FOOOFdata.fitError(Ind_currSubj);
            NewRow.Ce_intercept = FOOOFdata.intercept(Ind_currSubj);
            NewRow.Ce_slope = FOOOFdata.slope(Ind_currSubj);
            % Aperiodic power
                offset = FOOOFdata.intercept(Ind_currSubj);
                exponent = FOOOFdata.slope(Ind_currSubj);
                freqs_cur = 1:.5:30;
                Aperiodic = zeros(length(freqs_cur),1);
                for ii = 1:length(Aperiodic)
                    freq_cur = freqs_cur(1,ii);
                    Aperiodic(ii,1) = offset - log10(0+freq_cur^(exponent));
                end
                clear ii
                Aperiodic = Aperiodic';
            % Observed power
                Freq_lb = find(Powerb_nsoctrls.Pow_ft_struc.freq == freqs_cur(1,1));
                Freq_ub = find(Powerb_nsoctrls.Pow_ft_struc.freq == freqs_cur(1,end));
                PowLog10 = log10(Powerb_nsoctrls.Abspowspec.Cent(1,Freq_lb:Freq_ub));
            % For theta
                Freq_th_lb = find(freqs_cur == 4);
                Freq_th_ub = find(freqs_cur == 7);
                Obs_th_mn = mean(PowLog10(1,Freq_th_lb:Freq_th_ub),2);
                Aper_th_mn = mean(Aperiodic(1,Freq_th_lb:Freq_th_ub),2);
                NewRow.Ce_PerPower_th = Obs_th_mn - Aper_th_mn;
            % For alpha
                Freq_al_lb = find(freqs_cur == 8);
                Freq_al_ub = find(freqs_cur == 12);
                Obs_al_mn = mean(PowLog10(1,Freq_al_lb:Freq_al_ub),2);
                Aper_al_mn = mean(Aperiodic(1,Freq_al_lb:Freq_al_ub),2);
                NewRow.Ce_PerPower_al = Obs_al_mn - Aper_al_mn; 
           clear offset exponent Aperiodic PowLog10 Freq_lb Freq_ub Freq_th_lb Freq_th_ub 
           clear Freq_al_lb Freq_al_ub freq_cur
        else
            warning('Unable to find spectrum for current subject')
            NewRow.Ce_fitError = NaN;
            NewRow.Ce_intercept = NaN;
            NewRow.Ce_slope = NaN;
            NewRow.Ce_PerPower_th = NaN;
            NewRow.Ce_PerPower_al = NaN;    
        end   
        clear Ind_FOOOF Ind_currSubj SpectrumName 
        
        
     % Parietal data %%%%%%%%%%%%%%%%%%%%%%
        % find the row for subj 
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
            NewRow.Pa_fitError = FOOOFdata.fitError(Ind_currSubj);
            NewRow.Pa_intercept = FOOOFdata.intercept(Ind_currSubj);
            NewRow.Pa_slope = FOOOFdata.slope(Ind_currSubj);
            % Aperiodic power
                offset = FOOOFdata.intercept(Ind_currSubj);
                exponent = FOOOFdata.slope(Ind_currSubj);
                freqs_cur = 1:.5:30;
                Aperiodic = zeros(length(freqs_cur),1);
                for ii = 1:length(Aperiodic)
                    freq_cur = freqs_cur(1,ii);
                    Aperiodic(ii,1) = offset - log10(0+freq_cur^(exponent));
                end
                clear ii
                Aperiodic = Aperiodic';
            % Observed power
                Freq_lb = find(Powerb_nsoctrls.Pow_ft_struc.freq == freqs_cur(1,1));
                Freq_ub = find(Powerb_nsoctrls.Pow_ft_struc.freq == freqs_cur(1,end));
                PowLog10 = log10(Powerb_nsoctrls.Abspowspec.Pari(1,Freq_lb:Freq_ub));
            % For theta
                Freq_th_lb = find(freqs_cur == 4);
                Freq_th_ub = find(freqs_cur == 7);
                Obs_th_mn = mean(PowLog10(1,Freq_th_lb:Freq_th_ub),2);
                Aper_th_mn = mean(Aperiodic(1,Freq_th_lb:Freq_th_ub),2);
                NewRow.Pa_PerPower_th = Obs_th_mn - Aper_th_mn;
            % For alpha
                Freq_al_lb = find(freqs_cur == 8);
                Freq_al_ub = find(freqs_cur == 12);
                Obs_al_mn = mean(PowLog10(1,Freq_al_lb:Freq_al_ub),2);
                Aper_al_mn = mean(Aperiodic(1,Freq_al_lb:Freq_al_ub),2);
                NewRow.Pa_PerPower_al = Obs_al_mn - Aper_al_mn; 
           clear offset exponent Aperiodic PowLog10 Freq_lb Freq_ub Freq_th_lb Freq_th_ub 
           clear Freq_al_lb Freq_al_ub freq_cur
        else
            warning('Unable to find spectrum for current subject')
            NewRow.Pa_fitError = NaN;
            NewRow.Pa_intercept = NaN;
            NewRow.Pa_slope = NaN;
            NewRow.Pa_PerPower_th = NaN;
            NewRow.Pa_PerPower_al = NaN;    
        end   
        clear Ind_FOOOF Ind_currSubj SpectrumName 
        
    % Occipital data %%%%%%%%%%%%%%%%%%%%%%
        % find the row for subj 
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
            NewRow.Oc_fitError = FOOOFdata.fitError(Ind_currSubj);
            NewRow.Oc_intercept = FOOOFdata.intercept(Ind_currSubj);
            NewRow.Oc_slope = FOOOFdata.slope(Ind_currSubj);
            % Aperiodic power
                offset = FOOOFdata.intercept(Ind_currSubj);
                exponent = FOOOFdata.slope(Ind_currSubj);
                freqs_cur = 1:.5:30;
                Aperiodic = zeros(length(freqs_cur),1);
                for ii = 1:length(Aperiodic)
                    freq_cur = freqs_cur(1,ii);
                    Aperiodic(ii,1) = offset - log10(0+freq_cur^(exponent));
                end
                clear ii
                Aperiodic = Aperiodic';
            % Observed power
                Freq_lb = find(Powerb_nsoctrls.Pow_ft_struc.freq == freqs_cur(1,1));
                Freq_ub = find(Powerb_nsoctrls.Pow_ft_struc.freq == freqs_cur(1,end));
                PowLog10 = log10(Powerb_nsoctrls.Abspowspec.Occip(1,Freq_lb:Freq_ub));
            % For theta
                Freq_th_lb = find(freqs_cur == 4);
                Freq_th_ub = find(freqs_cur == 7);
                Obs_th_mn = mean(PowLog10(1,Freq_th_lb:Freq_th_ub),2);
                Aper_th_mn = mean(Aperiodic(1,Freq_th_lb:Freq_th_ub),2);
                NewRow.Oc_PerPower_th = Obs_th_mn - Aper_th_mn;
            % For alpha
                Freq_al_lb = find(freqs_cur == 8);
                Freq_al_ub = find(freqs_cur == 12);
                Obs_al_mn = mean(PowLog10(1,Freq_al_lb:Freq_al_ub),2);
                Aper_al_mn = mean(Aperiodic(1,Freq_al_lb:Freq_al_ub),2);
                NewRow.Oc_PerPower_al = Obs_al_mn - Aper_al_mn; 
           clear offset exponent Aperiodic PowLog10 Freq_lb Freq_ub Freq_th_lb Freq_th_ub 
           clear Freq_al_lb Freq_al_ub freq_cur
        else
            warning('Unable to find spectrum for current subject')
            NewRow.Oc_fitError = NaN;
            NewRow.Oc_intercept = NaN;
            NewRow.Oc_slope = NaN;
            NewRow.Oc_PerPower_th = NaN;
            NewRow.Oc_PerPower_al = NaN;    
        end   
        clear Ind_FOOOF Ind_currSubj SpectrumName 
        
        
        % add new data to tracker table and save 
            cd('XXX')
            if exist('SP_NSoc_PeriodicPower.mat','file') == 2
                load SP_NSoc_PeriodicPower.mat
            end
            TrackerNew = struct2table(NewRow);
            if ss == 1
                SP_NSoc_PeriodicPower = TrackerNew;
            else
                SP_NSoc_PeriodicPower(ss,:) = TrackerNew;
            end
            save('SP_NSoc_PeriodicPower.mat','SP_NSoc_PeriodicPower')
          % clear up
            clear SP_NSoc_PeriodicPower NewRow TrackerNew
        
        
end