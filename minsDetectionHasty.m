function [AllMins,AllMinsLoc,Freq]=minsDetectionHasty(T,Y)
% Pulse minimums detection algorithm for the gene oscillator
% and bacterial growth time series
% (identify minimums)

    figflaglbl=0;% Set this to 1 if you want to plot debugging figures, 0 otherwise

% (2015/06/04) NOTE THAT I'M ASSUMING THE TIME SERIES I'M PASSING TO THE
% PEAK DETECTION ALGORITHM HAS ALREADY BEEN EQUILIBRATED AND SAMPLED
% UNIFORMLY

    TT = T;
    YY = Y;

    % Number of steps that we need in order to define the minimum peak distance
    minpeakdist=10;%(2015/06/07)10;% CAREFUL HERE!!!! DEPENDING ON THE TIME-SERIES THIS 
    		   % VALUE COULD INTRODUCE A BIAS! Set 'figflaglbl=1' (line 5) and double
		   % check everything is correct (this is, that maximums are properly
		   % detected) in the debugging figures.

    % Next, we will compute the distribution of frequencies
    AllPeaks = [];
    AllMins = [];

    dim=size(Y);% Number of species
    for n = 1:dim(2)
        dummydeletion = [];
        % In order to find the maxima
        if length(YY(:,n)) < minpeakdist
            disp('No oscillations!');
            Freq(n) = 0;
        else
            clearvars pksmin locsmin pks locs
            [pks,locs]=findpeaks(YY(:,n),'MINPEAKHEIGHT',median(YY(:,n)),'MINPEAKDISTANCE',minpeakdist);
            % In order to find the minima
            [pksmin,locsmin]=findpeaks(-YY(:,n),'MINPEAKHEIGHT',20*median(-YY(:,n)),'MINPEAKDISTANCE',minpeakdist);
            if figflaglbl == 1
                figure
                hold on
                plot(TT(:),YY(:,n))
                scatter(TT(locsmin),YY(locsmin,n),'filled','r')
                scatter(TT(locs),YY(locs,n),'filled','g')
                hold off
            end
        
            clearvars mins.mxloc mins.minleftloc peak.minrightloc
            if length(pks) > 1
                j=1;
                for i=1:length(pks)
                    % Find the closest match for locs(i) in locsmin vector
                    [idleft,idright] = locate(locs(i), locsmin);
                    if idleft > 0 && idright >0
                        % Creating the structure that will contain the location for the maxs
                        % as well as the left and right minimums locations
                        peak.mxloc(j)=locs(i);
                        peak.minleftloc(j)=locsmin(idleft);
                        peak.minrightloc(j)=locsmin(idright);
                        j=j+1;
                    end
                end
                AllPeaksORI{n}=peak.mxloc;
                
                deletion=[];
                % Now it's time to get rid of false peaks!
                for i=1:(length(peak.mxloc)-1)
                    if (peak.minleftloc(i) == peak.minleftloc(i+1)) && (peak.minrightloc(i) == peak.minrightloc(i+1))
                        % We keep the maximum peak and discard the other one
                        if YY(peak.mxloc(i),n) > YY(peak.mxloc(i+1),n)
                            %                 disp('deleting element of the peak structure')
                            %                 i+1
                            deletion=cat(1,deletion,i+1);
                        else
                            %                 disp('deleting element of the peak structure')
                            %                 i
                            deletion=cat(1,deletion,i);
                        end
                    end
                    % Check also that the difference between the max and its minimums
                    if (YY(peak.mxloc(i),n)-YY(peak.minleftloc(i),n)) <= 0.1*YY(peak.mxloc(i),n) || ...
                            (YY(peak.mxloc(i),n)-YY(peak.minrightloc(i),n)) <= 0.1*YY(peak.mxloc(i),n)
                        deletion=cat(1,deletion,i);
                    end

                end

                deletion=unique(deletion);% To avoid repetitions

                peak.mxloc(deletion')=[];
                peak.minleftloc(deletion')=[];
                peak.minrightloc(deletion')=[];

                
                if length(peak.mxloc) >= 1
                    % Finally, discard pulses with amplitudes <= 30% of the mean
                    % amplitude
                    for i=1:length(peak.mxloc)
                        Amps{n}(i)= max([YY(peak.mxloc(i),n)-YY(peak.minleftloc(i),n)...
                            YY(peak.mxloc(i),n)-YY(peak.minrightloc(i),n)]);
                    end
                    meanAmps(n) = mean(Amps{n});
                    deletion=[];
                    for i=1:length(peak.mxloc)
                        if Amps{n}(i) < 0.10*meanAmps(n) 
                            deletion=cat(1,deletion,i);
                        end
                    end

                    deletion=unique(deletion);% To avoid repetitions
                    % Storing rejected peaks for debugging purposes
                    dummydeletion=cat(1,dummydeletion,peak.mxloc(deletion')');
                    AllPeaksDel{n}=dummydeletion';

                    peak.mxloc(deletion')=[];
                    peak.minleftloc(deletion')=[];
                    peak.minrightloc(deletion')=[];

                    if figflaglbl == 1
                        figure
                        hold on
                        plot(TT(:),YY(:,n))
                        scatter(TT(peak.mxloc),YY(peak.mxloc,n),'filled','r')
                        hold off
                    end

                    % And now you have the peak structure ready to extract information from
                    % pulses.
                    for i=1:length(peak.mxloc)
                        AllPeaks{n}(i) = TT(peak.mxloc(i));
                        AllPeaksLoc{n}(i) = peak.mxloc(i);
                    end
                        
                    L(n) = length(peak.mxloc);%Then L(n) will give us the number of pulses per specie

                    Freq(n) = L(n)*3600.0 / (max(TT(:))-min(TT(:)));% In h^-1
                else
                    disp(sprintf('No pulses for specie %s\n',num2str(n)));
                    Freq(n) = 0;
                end
            else
                disp(sprintf('No pulses for specie %s\n',num2str(n)));
                Freq(n) = 0;
            end
            if length(pksmin) > 1
               j=1;
               for i=1:length(pksmin)
                   [idleftmin,idrightmin] = locate(locsmin(i), locs);
                   if idleftmin > 0 && idrightmin >0
                       
                       % Creating the structure that will contain the location
                       % for the mins
                       % as well as the left and right maximums locations
                       peak.minloc(j)=locsmin(i);
                       peak.maxleftloc(j)=locs(idleftmin);
                       peak.maxrightloc(j)=locs(idrightmin);
                       j=j+1;
                   end
               end
               AllMinsORI{n}=peak.minloc;
            
               deletion=[];
               % Now it's time to get rid of false minimums! These are
               % those minimums that are not surrounded by two maximums
               for i=1:(length(peak.minloc)-1)
                   if (peak.maxleftloc(i) == peak.maxleftloc(i+1)) && (peak.maxrightloc(i) == peak.maxrightloc(i+1))
                        % We keep the lowest minimum and discard the other one
                        if YY(peak.minloc(i),n) < YY(peak.minloc(i+1),n)
                            disp('deleting element of the peak structure')
                            i+1
                            deletion=cat(1,deletion,i+1);
                        else
                            disp('deleting element of the peak structure')
                            i
                            deletion=cat(1,deletion,i);
                        end
                   end
                end

                deletion=unique(deletion);% To avoid repetitions

                peak.minloc(deletion')=[];
                peak.maxleftloc(deletion')=[];
                peak.maxrightloc(deletion')=[];
               
               
               for i=1:length(peak.minloc)
                   AllMins{n}(i) = TT(peak.minloc(i));
                   AllMinsLoc{n}(i) = peak.minloc(i);
               end
               
            end
            
          end
            
        end
    
    if figflaglbl == 1
            n=1;
            figure
            hold on
            plot(TT(:),YY(:,n))
            scatter(TT(locsmin),YY(locsmin,n),'filled','r')
            scatter(TT(AllPeaksLoc{n}),YY(AllPeaksLoc{n},n),'filled','g')
            hold off
            
            figure
            hold on
            plot(TT(:),YY(:,n))
            scatter(TT(AllMinsLoc{n}),YY(AllMinsLoc{n},n),'filled','k')
            scatter(TT(AllPeaksLoc{n}),YY(AllPeaksLoc{n},n),'filled','g')
            hold off
    end

end

