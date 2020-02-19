function [data]=cleanTMS(cfg,data)


%
% Cleans TMS artefact from EEG-Data
%
%cfg.trials: 'all' or vecotor with trial indices
%cfg.mindistance: minimum distance between peaks (in points; default: length(cfg.time/ 2))
%cfg.ampsd: minimum standard deviation for to-be detected peaks (default:1)
%cfg.plottrial: show data and detected peaks (default: 'no')
%cfg.method: 'single' for Single-Pulse TMS (default) or 'rtms' for rTMS-Pulses
%cfg.interp: 'spline' (spline interpolation), 'NaN' (replacing with NaNs) or 'std' (default, fill up with noise)
%cfg.pretime: start latency of period to be replaced (in ms: default 5)
%cfg.posttime: stop latency to be replaced (in ms: default: 15)
%cfg.refwin: time window used to replace data (in ms: default -200 - -20), not needed for spline interpolation 
%-----
%output: x=data.trial ohne artefakt
%
%
%----- No longer in Use:
%cfg.input = data.trial; (now deprecated, jk 2012)
%cfg.time = data.time; (now deprecated, jk 2012)
%cfg.t_tms = where is the tms? (now deprecated, nw 2011)
%
% Hannah Schulz 2011
% Nathan Weisz 2011
% Julian Keil 2011, 2012

% Revision 10.01.2012: Cleaned up the code and added interp options
% Revision 27.04.2012: Added the NaN-Option and automatic peak detection to
% the highest peak in the Data (jk)

%%
% how many trials?
if ~isfield(cfg,'trials')
    trials=1:length(data.trial); % Vector for all trials 
else
    trials=cfg.trials; % vector for subset of trials
end

% Amplitude Standard Deviation Minimum
if ~isfield(cfg, 'ampsd')
    ampsd=1;
else
    ampsd=cfg.ampsd;
end

% Should the trials be plotted?
if isfield(cfg, 'plottrial')
    if strcmpi(cfg.plottrial,'yes')
        plottrial=1;
    else
        plottrial=0;
    end
else
    plottrial=0;
end

% Single or multiple Pulses
if isfield(cfg, 'method')
    if strcmpi(cfg.method,'single')
        singlepulse=1;
    else
        singlepulse=0;
    end
else
    singlepulse=1;
end

% How mmuch time before TM puls should be replaced
if isfield(cfg, 'pretime')
    pretime=cfg.pretime; % revised: shouldn't be divided by 1000, 
                         % otherwise, conversion to sample points
                         % doesn't work anymore (round(pretime/(1000/fs)), ca. line 210) 
else
    pretime=5; % 
end

% How much after?
if isfield(cfg, 'posttime')
    posttime=cfg.posttime; % revised: shouldn't be divided by 1000 (see above)
else
    posttime=15; 
end

% What data should be used to interpolate?
if isfield(cfg, 'refwin')
    refwin=abs(cfg.refwin)/1000;
else
    refwin=[200 20]/1000;
end
refwin=abs(refwin);

% How should the data be interpolated?
if isfield(cfg, 'interp')
    if strcmpi(cfg.interp,'spline')
        spl=1;
    elseif strcmpi(cfg.interp,'std')
        spl=0;
    elseif strcmpi(cfg.interp, 'NaN')
        spl=2;
    end
else
    spl=0;
end
%%
%t_tms = cfg.t_tms;
%%
for ltr=trials % Loop around trials defined in cfg.trials
    
    fprintf('Removing TMS from Trial: %i\n',ltr)
    
    y = data.trial{ltr}; % data in cell ltr 
    t = data.time{ltr};  % samples in cell ltr

%%
        if ~isfield(cfg, 'mindistance') % minimum samples between peaks in 'findpeaks': 1499.5
            mindistance=length(t)/2;    % could be set too high
        else
            mindistance=cfg.mindistance;% might be the reason TM pulses ain't removed
        end
        
%% find peaks
    z=mean(abs(y),1); % column-wise mean over electrodes for each sample


%% minimum amplitude set rest to zero

    zthresh=z.*(z >= mean(z)+ampsd*std(z)); % threshold: values 'ampsd' away from mean(z), rest is set to 0  
    
%%

    % 'findpeaks' only takes values above threshold and selects tallest
       [pks, locs] = findpeaks(zthresh,'MINPEAKDISTANCE',mindistance);%,'THRESHOLD',10); 

        if plottrial == 1
            ff=figure;
            plot(t,z);hold;plot(t(locs),pks,'g*')
            waitfor(ff)
        end

%%
        if singlepulse == 1        % does 'singlepulse' mean one pulse per trial?
            [s_pks l_pks]=max(pks);% find the largest peak in the data
            t_tms=t(locs(l_pks));  % set tms_time to time of largest peak
        else
            t_tms=t(locs); 
        end
%%
        if spl == 1 % Thut & Gross-Approach
            
              for ipulse=1:length(t_tms)   
                  
                  % 1. Replace TMS-Time with NaNs

                  x=y;
                  fs=1/(t(2)-t(1));

                  prewin=round(pretime/(1000/fs)); %10ms in samle points
                  postwin=round(posttime/(1000/fs)); %15 in samle points
                  maxind=nearest(t,t_tms(ipulse));
                  
                  x(:, (maxind-prewin):(maxind+postwin))=NaN;
                  x_int=x;


                  for chan=1:size(x,1) % Loop for Channels -> Can someone vectorize this?
                      
                      % 2. Define length of Data
                      time = 1:length(x(chan,:));

                      % 3. Define Mask for NaNs
                      mask = ~isnan(x(chan,:));

                      % 4. Interpolate all the NaNs
                      
                      x_int(chan,~mask) = interp1(time(mask),x(chan,mask),time(~mask),'spline');
                  
                  end % channel loop
                  data.trial{ltr}=x_int;
             
              end % t_tms
              
        elseif spl == 2 % Replace with NaNs-Approach
            
              for ipulse=1:length(t_tms)   
                  
                  % 1. Replace TMS-Time with NaNs

                  x=y;
                  fs=1/(t(2)-t(1));

                  prewin=round(pretime/(1000/fs)); %10ms in samle points
                  postwin=round(posttime/(1000/fs)); %15 in samle points
                  maxind=nearest(t,t_tms(ipulse));
                  
                  x(:, (maxind-prewin):(maxind+postwin))=NaN;

                  data.trial{ltr}=x;
             
              end % t_tms

              
        elseif spl == 0 % Replace with Noise-Approach
      
            x=y; % duplicate data
            fs=1/(t(2)-t(1)); % what does this do? sample rate in Hz?

            % the following converts pre-/posttime from ms in sample points
            prewin=round(pretime/(1000/fs)); %10ms in samle points 
            postwin=round(posttime/(1000/fs)); %15 in samle points

              for ipulse=1:length(t_tms) % # of iterations is currently set to 1 
                  
                           % preparations:
                           % defines indices for reference period (refwin limits and max peak) 
                           prestdind1=nearest(t,t_tms(ipulse)-refwin(1)); %zeiten ersetzen durch variablen
                           prestdind2=nearest(t,t_tms(ipulse)-refwin(2));
                           maxind=nearest(t,t_tms(ipulse));
                           
                           ss=std(y(:,prestdind1:prestdind2)'); % STD for every channel

                           repvec=[(maxind-prewin):(maxind+postwin)]; % Timevector (sample indices) for replacement
                           
                           % generate noise:
                           % random numbers * STD + mean of reference
                           % period (z-transform ?)
                           tmpx=(randn(size(y,1),length(repvec)).*repmat(ss,length(repvec),1)')...
                               +repmat(mean(y(:,prestdind1:prestdind2),2)',length(repvec),1)'; % evt offset berücksichtigen
            
                    % fill up the reference period with generated noise
                           x(:, (maxind-prewin):(maxind+postwin))=tmpx; 
                    data.trial{ltr}=x;
              end % t_tms
              
        end% if spline
     
end % trials

          