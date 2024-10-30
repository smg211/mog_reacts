function [state_out, signal_thresh, noise_thresh, min_upstate_dur, min_downstate_dur] ...
  = find_binary_state(cfg, dat)
% Given a 1-dimensional continuous signal, find when the signal crosses
% over the signal threshold and subsequently stays over the noise
% threshold for a minimum number of samples

doclose_checkplot   = ft_getopt(cfg, 'doclose_checkplot', true);
dosave_checkplot    = ft_getopt(cfg, 'dosave_checkplot', false);
checkplot_savpath   = ft_getopt(cfg, 'checkplot_savpath', {});
signal_thresh       = ft_getopt(cfg, 'signal_thresh', []);
noise_thresh        = ft_getopt(cfg, 'noise_thresh', signal_thresh);
end_up_when_decrease= ft_getopt(cfg, 'end_up_when_decrease', false);
min_upstate_dur     = ft_getopt(cfg, 'min_upstate_dur', 0); % samples, the minimum duration allowed for an above threshold period to qualify as an up state
min_downstate_dur   = ft_getopt(cfg, 'min_downstate_dur', 0); % samples, the minimum duration allowed for an above threshold period to qualify as a down state
docheck             = ft_getopt(cfg, 'docheck', false);
dat4checkplot       = ft_getopt(cfg, 'dat4checkplot', {});
dat4checkplot_str   = ft_getopt(cfg, 'dat4checkplot_str', {});
time4checkplot      = ft_getopt(cfg, 'time4checkplot', 1:length(dat));
n_plotwindows       = ft_getopt(cfg, 'n_plotwindows', 5);
plotwinsize         = ft_getopt(cfg, 'plotwinsize', 'auto');
allplotscontinuous  = ft_getopt(cfg, 'allplotscontinuous', false);
ylim_method         = ft_getopt(cfg, 'ylim_method', 'thresh'); % 'thresh' or 'abs'
i_winedges          = ft_getopt(cfg, 'i_winedges', []);

if docheck || dosave_checkplot
  if ~iscell(dat4checkplot)
    dat4checkplot = {dat4checkplot};
  end

  if length(dat4checkplot) ~= length(dat4checkplot_str)
    tmp = length(dat4checkplot_str);
    for i = 1:length(dat4checkplot) - length(dat4checkplot_str)
      dat4checkplot_str{tmp + i} = ['Associated Timeseries ' num2str(tmp + i)];
    end
  end
end

docontinue = 0;
retry = 0;
while ~docontinue
  abv_sig = dat > signal_thresh;
  i_bel2abv_sig = find(diff([0 abv_sig]) == 1);
  
  abv_nz = dat > noise_thresh;
  i_bel2abv_nz = find(diff([0 abv_nz]) == 1);
  i_abv2bel_nz = find(diff([0 abv_nz]) == -1);
  

  state_out = zeros(1, length(dat));
  i_up_end_last = 0;
  for i = 1:length(i_bel2abv_sig) % for each signal thresh crossing
    % find the previous noise thresh crossing
    i_up_start = i_bel2abv_nz(find(i_bel2abv_nz <= i_bel2abv_sig(i), 1, 'last'));
    
    if i_up_start > i_up_end_last % make sure this crossing is after the last up state ended
%       i_up_start = i_bel2abv_sig(i);
      
      % find the next time the signal goes below the noise threshold
      if any(i_abv2bel_nz > i_up_start)
        if end_up_when_decrease
          % the up state ends when the signal starts to decrease or becomes nan (note:
          % this was implemented for saccade detection under the assumption
          % that the eyes are accelerating or at a constant velocity
          % throughout the saccade)
          i_up_end = find((diff([0 dat]) < 0 | isnan(dat)) & 1:length(dat) > i_up_start, 1, 'first');
        else
          i_up_end = i_abv2bel_nz(find(i_abv2bel_nz > i_up_start, 1));
        end

        % if the end of this state is at least min_samp_dur away, then count
        % this as a state
        if i_up_end - i_up_start >= min_upstate_dur

          if ~end_up_when_decrease && any(i_bel2abv_sig > i_up_end)
            % if the start of the next up state is not at least min_downstate_dur
            % away from the end of this up state, then consider the next up state
            % and the time between the 2 up states as part of the same up state
            ix = find(i_bel2abv_sig > i_up_end, 1, 'first');
            i_next_state_start = i_bel2abv_nz(find(i_bel2abv_nz <= i_bel2abv_sig(ix), 1, 'last'));
            while (i_next_state_start - i_up_end < min_downstate_dur) && ... % while the next up state is less than min_down_dur away
                any(i_abv2bel_nz > i_next_state_start) % && ...
%                 (i_abv2bel_nz(find(i_abv2bel_nz > i_bel2abv_sig(ix), 1)) - i_bel2abv_sig(ix)) > min_upstate_dur % and the next up state is at least min_up_dur long
              i_up_end = i_abv2bel_nz(find(i_abv2bel_nz > i_bel2abv_sig(ix), 1));
              if isempty(i_up_end)
                i_up_end = length(state_out);
              end
              ix = ix + 1;
              
              if ix+1 > length(i_bel2abv_sig)
                break
              end
              
              i_next_state_start = i_bel2abv_nz(find(i_bel2abv_nz < i_bel2abv_sig(ix), 1, 'last'));
            end
          else
            ix = i;
          end

          if ix <= length(i_bel2abv_sig) % ix+1 <= length(i_bel2abv_sig)
            state_out(i_up_start:i_up_end) = 1;
          end
        end
      else
        i_up_end = length(abv_sig);
      end

      i_up_end_last = i_up_end;
    end
  end
  
  state_out(isnan(dat)) = 0;
  
  % Plot the result and check to see if it is acceptable
  if docheck || dosave_checkplot
    fig = figure('units','normalized','outerposition',[0 0 1 1]);
    if ~docheck && doclose_checkplot
      fig.Visible = 'off';
    end
    
    plot_dims(1) = min([n_plotwindows 5]);
    plot_dims(2) = ceil(n_plotwindows/5);
    
    if isempty(i_winedges)
      if strcmp(plotwinsize, 'auto')
        plotwinsize = 50*max([min_upstate_dur min_downstate_dur]);
      end
      if plotwinsize == 0
        plotwinsize = length(dat)/100;
      end
      % get the edges of the windows we are going to plot
    
      i_winedges = [round(linspace(1, length(dat)-plotwinsize, n_plotwindows))' ...
        round(linspace(plotwinsize, length(dat), n_plotwindows))'];
    end
    
    % get the edges of the up states
    i_up_edges = [find(diff([0 state_out]) == 1)' find(diff([state_out 0]) == -1)'];
    for w = 1:n_plotwindows
      if allplotscontinuous
        i_samp2plot = 1:length(dat);
      else
        i_samp2plot = i_winedges(w, 1):i_winedges(w, 2);
      end
      
      h_leg = [];
      subplot(plot_dims(1), plot_dims(2), w); hold on;

      % plot the data in black
      h_leg(end+1) = plot(time4checkplot(i_samp2plot), dat(i_samp2plot), '-k', 'LineWidth', 2);
      
      % if there is any additional timeseries to be plotted along with the
      % main input data, plot it in gray
      if ~isempty(dat4checkplot)
        for i = 1:length(dat4checkplot)
          h_leg(end+1) = plot(time4checkplot(i_samp2plot), dat4checkplot{i}(i_samp2plot), '-', 'Color', [0.5 0.5 0.5], 'LineWidth', 1);
        end
      end
      
      % plot the signal threshold in blue
      if retry
        h_leg(end+1) = plot(time4checkplot(i_samp2plot([1 end])), [signal_thresh_last signal_thresh_last], '--b', 'LineWidth', 1);
      end
      h_leg(end+1) = plot(time4checkplot(i_samp2plot([1 end])), [signal_thresh signal_thresh], '-b', 'LineWidth', 1);
      
      % plot the noise threshold in red
      if retry
        h_leg(end+1) = plot(time4checkplot(i_samp2plot([1 end])), [noise_thresh_last noise_thresh_last], '--r', 'LineWidth', 1);
      end
      h_leg(end+1) = plot(time4checkplot(i_samp2plot([1 end])), [noise_thresh noise_thresh], '-r', 'LineWidth', 1);
      
      % highlight the up state periods in green
      a = gca;
      ylims = a.YLim;
      if retry
        for i = 1:size(i_up_edges_last, 1)
          if i == 1
%             h_leg(end+1) = plot(i_up_edges_last(i, :), [signal_thresh([1 1])+0.5*signal_thresh]-0.1*diff(get(gca, 'YLim')), '-', 'Color', [0 0.5 0], 'LineWidth', 3);
            h_leg(end+1) = patch(time4checkplot(i_up_edges_last(i, [1 1 2 2])), ylims([1 2 2 1]), ...
                [0 0 1], 'FaceAlpha', 0.1, 'EdgeAlpha', 0);
          else
            if any(i_up_edges_last(i, :) >= i_samp2plot(1) & i_up_edges_last(i, :) <= i_samp2plot(end))
%               plot(time4checkplot(i_up_edges_last(i, :)), ...
%                 [signal_thresh([1 1])+0.5*signal_thresh]-0.1*diff(get(gca, 'YLim')), ...
%                 '-', 'Color', [0 0.5 0], 'LineWidth', 3);
              patch(time4checkplot(i_up_edges_last(i, [1 1 2 2])), ylims([1 2 2 1]), ...
                [0 0 1], 'FaceAlpha', 0.1, 'EdgeAlpha', 0);
            end
          end
        end
      end
      for i = 1:size(i_up_edges, 1)
        if i == 1
%           h_leg(end+1) = plot(i_up_edges(i, :), [signal_thresh([1 1])+0.5*signal_thresh], '-g', 'LineWidth', 3);
          h_leg(end+1) = patch(time4checkplot(i_up_edges(i, [1 1 2 2])), ylims([1 2 2 1]), ...
                [0 1 0], 'FaceAlpha', 0.1, 'EdgeAlpha', 0);
        else
          if any(i_up_edges(i, :) >= i_samp2plot(1) & i_up_edges(i, :) <= i_samp2plot(end))
%             plot(time4checkplot(i_up_edges(i, :)), [signal_thresh([1 1])+0.5*signal_thresh], ...
%               '-g', 'LineWidth', 3);
            patch(time4checkplot(i_up_edges(i, [1 1 2 2])), ylims([1 2 2 1]), ...
                [0 1 0], 'FaceAlpha', 0.1, 'EdgeAlpha', 0);
          end
        end
      end

      a = gca;
      a.XLim = time4checkplot(i_winedges(w, :));
      if ~isempty(dat4checkplot)
        if strcmp(ylim_method, 'abs')
          a.YLim = [min([dat [dat4checkplot{:}]]) max([dat [dat4checkplot{:}]])+0.05*range([dat [dat4checkplot{:}]])];
        elseif strcmp(ylim_method, 'thresh')
          a.YLim = [-signal_thresh-0.5*signal_thresh signal_thresh+0.5*signal_thresh];
        end
      else
        a.YLim = [min(dat) max(dat)+0.05*range(dat)];
      end
      
      if w == 1
        leg = legend(h_leg);
        if retry
          leg_string = {'Input Timeseries', 'Signal Thresh (previous)', 'Signal Thresh (current)', 'Noise Thresh (previous)', 'Noise Thresh (current)', 'Up State (previous)', 'Up State (current)'};
        else
          leg_string = {'Input Timeseries', 'Signal Thresh', 'Noise Thresh', 'Up State'};
        end
        if ~isempty(dat4checkplot)
          leg_string = [leg_string(1) dat4checkplot_str leg_string(2:end)];
        end
        leg.String = leg_string;
        leg.Location = 'NorthEast';
        
        if retry
          a.Title.String = {['Up Duration Thresh = ' num2str(min_upstate_dur_last) ' (previous) --> ' num2str(min_upstate_dur) ' (current) samples'], ...
            ['Down Duration Thresh = ' num2str(min_downstate_dur_last) ' (previous) --> ' num2str(min_downstate_dur) ' (current) samples']};
        else
          a.Title.String = {['Up Duration Thresh = ' num2str(min_upstate_dur) ' samples'], ...
            ['Down Duration Thresh = ' num2str(min_downstate_dur) ' samples']};
        end
      end
      a.YLim = ylims;
    end
    
    if docheck
      docontinue = ~input('Would you like to change the thresholds for binary state detection? \n true or false: ');
      if ~docontinue
        retry = 1;
        i_up_edges_last = i_up_edges;
        signal_thresh_last = signal_thresh;
        noise_thresh_last = noise_thresh;
        min_upstate_dur_last = min_upstate_dur;
        min_downstate_dur_last = min_downstate_dur;

        signal_thresh = input(['\n The previous signal threshold (units of raw data) was ' num2str(signal_thresh) '\n' ...
          'Enter the new signal threshold (or enter "same"): ']);
        if strcmp(signal_thresh, 'same')
          signal_thresh = signal_thresh_last;
        end

        noise_thresh = input(['\n The previous noise threshold (units of raw data) was ' num2str(noise_thresh) '\n' ...
          'Enter the new noise threshold (or enter "same"): ']);
        if strcmp(noise_thresh, 'same')
          noise_thresh = noise_thresh_last;
        end

        min_upstate_dur = input(['\n The previous upstate duration threshold was ' num2str(min_upstate_dur) '\n' ...
          'Enter the new upstate duration threshold (in # of samples) (or enter "same"): ']);
        if strcmp(min_upstate_dur, 'same')
          min_upstate_dur = min_upstate_dur_last;
        end

        min_downstate_dur = input(['\n The previous downstate duration threshold was ' num2str(min_downstate_dur) '\n' ...
          'Enter the new downstate duration threshold (in # of samples) (or enter "same"): ']);
        if strcmp(min_downstate_dur, 'same')
          min_downstate_dur = min_downstate_dur_last;
        end
      end
    else
      docontinue = 1;
    end
    if dosave_checkplot
      print(fig, checkplot_savpath, '-dpng');
    end
    if doclose_checkplot
      close(fig);
    end
  else
    docontinue = 1;
  end
end



