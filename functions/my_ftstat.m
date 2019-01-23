% Inputs
%   Y: (channels, times, conds, subjects)
%   t: times (ms)
%   urevent: ...
%   Nsuj: # of subjects
%   CONFIG: Params of the original FieldTrip function 
%       http://www.fieldtriptoolbox.org/tutorial/cluster_permutation_timelock

function [stat, times, ftdata, grandavg, cfg] = my_ftstat(Y,t,urevent,subjlist,CONFIG, CHANS)

Nsuj = length(subjlist);

%% eeglab2fieldtrip
    clear ftdata
    for su = 1:Nsuj
        for ico = 1:2;
            if ~isempty(CONFIG.EEGBLANK);
                EEG = CONFIG.EEGBLANK;
            else
                EEG = eegblank;
            end

            data = squeeze(Y(:,:,ico,su));

            EEG.setname = ['ave' num2str(ico)];
            EEG.nbchan  = size(data,1);
            EEG.trials  = size(data,3);
            EEG.pnts    = length(t);
            EEG.xmin    = round(min(t))/1000;  % in seconds
            EEG.xmax    = round(max(t))/1000;  % in seconds
            EEG.times   = t;            % in milliseconds
            EEG.srate   = round(mean(1000./diff(t)));
            EEG.data    = data;

            EEG.chanlocs= CHANS.chanlocs;
            EEG.chaninfo= CHANS.chaninfo;
            EEG.urchanlocs= CHANS.urchanlocs;
%             EEG         = eeg_checkset(EEG);  

            ftdata(su,ico) = eeglab2fieldtrip(EEG,'timelockanalysis','none' );
            ftdata(su,ico).time = t;
        end
    end
    elec    = ftdata(1,1).elec;
    times   = ftdata(1,1).time;

%     keyboard
%% average

    cfg = [];
    cfg.channel        = 'all'; % Nx1 cell-array with selection of channels (default = 'all'), see FT_CHANNELSELECTION for details
    cfg.latency        = 'all'; % [begin end] in seconds or 'all' (default = 'all')
    cfg.keepindividual = 'yes'; % 'yes' or 'no' (default = 'no')
    cfg.normalizevar   = 'N-1'; %'N' or 'N-1' (default = 'N-1')

    for ity=1:2;
        str = [];
        for su =1:Nsuj
            ftdata(su,ity).fsample = EEG.srate;
            ftdata(su,ity).label = elec.label';
            if ~isempty(urevent)
                ftdata(su,ity).dof = length(urevent{su,ity})*ones(size(ftdata(su,ity).avg)); % puedo estimarlo de antes... cantidad de trials que lo componen
            else
                ftdata(su,ity).dof = ones(size(ftdata(su,ity).avg));
            end
            ftdata(su,ity).dimord = 'chan_time';
            str=[str,'ftdata(',num2str(su),',',num2str(ity),'),'];
        end
        str(end)=[];    

        eval(['grandavg(',num2str(ity),')','=','ft_timelockgrandaverage(cfg,',str,')' ';']);
    end

    for ity=1:2;
        grandavg(ity).fsample   = ftdata(1,1).fsample;
        grandavg(ity).avg       = squeeze(mean(grandavg(ity).individual,1));
        grandavg(ity).elec      = ftdata(1,1).elec;
    end

%% define neighborhood
%     % Esto se parece mas a lo que tenia
%     cfg                 = [];
%     cfg.method        = 'distance';
%     cfg.elec            = elec;    
%     cfg.neighbourdist   = 0.35;
%     neighbours = ft_neighbourselection(cfg, grandavg(1));

%     % Esto me gusta mas ahora
%     cfg                 = [];
%     cfg.method          = 'triangulation';
%     cfg.elec            = elec;
%     neighbours = ft_neighbourselection(cfg, grandavg(1));

    % Nueva version
    % prepare_neighbours determines what sensors may form clusters
    cfg_neighb.method    = 'distance';
    neighbours       = ft_prepare_neighbours(cfg_neighb, grandavg(1));

%     % Para chequearlo
%     cfg                 = [];
%     cfg.elec            = elec;
%     cfg.neighbours      = neighbours;
%     ft_neighbourplot(cfg, grandavg(1));

%% Estadistica
    cfg                 = [];
    cfg.channel         = 'all'; %'EEG1010' % see CHANNELSELECTION
    cfg.latency         = 'all';

    cfg.method          = 'montecarlo'; % use the Monte Carlo Method to calculate the significance probability
    cfg.statistic       = 'depsamplesT';% use the independent samples T-statistic as a measure to evaluate 
                                        % the effect at the sample level
%     cfg.statistic       = 'ft_statfun_depsamplesT';% use the independent samples T-statistic as a measure to evaluate 
%                                         % the effect at the sample level

    % cfg.correctm        = 'fdr';
    cfg.correctm        = 'cluster';
    if isfield(CONFIG,'clusteralpha')
        cfg.clusteralpha= CONFIG.clusteralpha; 
    else
        cfg.clusteralpha= 0.01;         % alpha level of the sample-specific test statistic that will be used for thresholding
    end
    
    cfg.clusterstatistic= 'maxsum';     % test statistic that will be evaluated under the permutation distribution.
    if isfield(CONFIG,'minnbchan')
        cfg.minnbchan   = CONFIG.minnbchan;
    else
        cfg.minnbchan   = 2;            % minimum number of neighborhood channels that is required for a selected
                                        % sample to be included in the clustering algorithm (default=0).
    end
    cfg.tail            = 0;            % -1, 1 or 0 (default = 0); one-sided or two-sided test
    cfg.correcttail     = 'alpha';      % Note that if you want to run a two-sided test, 
                                        % you have to split the critical alpha value by 
                                        % setting cfg.correcttail = 'alpha'; i.e. this 
                                        % sets cfg.alpha = 0.025, corresponding to a 
                                        % false alarm rate of 0.05 in a two-sided test. 
                                        % The field cfg.alpha is not crucial.
    cfg.alpha           = 0.05;         % alpha level of the permutation test, this is no 
    cfg.clustertail     = 0;
    if isfield(CONFIG,'numrandomization')
        cfg.numrandomization= CONFIG.numrandomization;
    else
        cfg.numrandomization= 500;      % number of draws from the permutation distribution
    end
    
    % subject number
    subj = Nsuj;
    design = zeros(2,2*subj);
    for i = 1:subj
      design(1,i)       = i;
    end
    for i = 1:subj
      design(1,subj+i)  = i;
    end

    % condition number
    design(2,1:subj)        = 1;
    design(2,subj+1:2*subj) = 2;

    cfg.design          = design;
    cfg.uvar            = 1;            % "subject" is unit of observation
    cfg.ivar            = 2;            % "condition" is the independent variable
    fprintf('JK: Number of units = %d\n',   length(unique(cfg.design(cfg.uvar,:))))
    fprintf('JK: Number of levels = %d\n',  length(unique(cfg.design(cfg.ivar,:))))
    
    cfg.avgovertime     = 'no';
    cfg.avgoverchan     = 'no';

    cfg.elec            = elec;
    cfg.neighbours      = neighbours;

    stat = ft_timelockstatistics(cfg, grandavg(1), grandavg(2));
end

function EEG = eegblank()
    EEG.setname='';
    EEG.filename='';
    EEG.filepath='';
    EEG.subject='';
    EEG.group='';
    EEG.condition='';
    EEG.session=[];
    EEG.comments='';
    EEG.nbchan=0;
    EEG.trials=0;
    EEG.pnts=0;
    EEG.srate=1;
    EEG.xmin=0;
    EEG.xmax=0;
    EEG.times=[];
    EEG.data=[];
    EEG.icaact=[];
    EEG.icawinv=[];
    EEG.icasphere=[];
    EEG.icaweights=[];
    EEG.icachansind=[];
    EEG.chanlocs=[];
    EEG.urchanlocs=[];
    EEG.chaninfo=[];
    EEG.ref=[];
    EEG.event=[];
    EEG.urevent=[];
    EEG.eventdescription={};
    EEG.epoch=[];
    EEG.epochdescription={};
    EEG.reject=[];
    EEG.stats=[];
    EEG.specdata=[];
    EEG.specicaact=[];
    EEG.splinefile='';
    EEG.icasplinefile='';
    EEG.dipfit=[];
    EEG.history='';
    EEG.saved='no';
    EEG.etc=[];
end