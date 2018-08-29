function processingSimData_Portfolio(TTT,YYY,dataID)
% Processing oscillatory time-series. Using "noisy" simulated data for a synthetic gene oscillator
% in bacteria (simulating gene expression for a gene oscillator, and simulating bacterial
% cell growth having bacterial length as a reporter), this code detects local maximums 
% for both time-series and computes their periods as well as the 
% the following histograms/analysis:
%  - distribution of bacterial length periods (unidirectionalCoupling_LENperiodshistogram.pdf)
%  - distribution of gene oscillator periods (unidirectionalCoupling_OSCperiodshistogram.pdf)
%  - "Pulse triggered average" (PTA) (unidirectionalCoupling_PTAcellDivisionrawdata.pdf)
%    (each segment is colored in function of its "age" (blue corresponds to younger generations,
%    while red corresponds to older generations).
%  - PTA cell division raw data -> (unidirectionalCoupling_PTAcellDivision_OscMAXShistogramNormt.pdf)
%

% Simulation data contains several arrays from which we need: 
% TTT: time (in seconds)
% YYY: time-series for the different species/observables (YYY)
%      Specifically, YYY(:,1) corresponds to the gene oscillator reporter (in a.u.)
%      and YYY(:,3) corresponds to bacterial cell length (in a.u.)

    wup = 100*3600;% warm up time (CAREFUL HERE! IT NEEDS TO BE THE
    % SAME VALUE AS IN 'peakDetectionHasty.m' code!)
    % we need to dump the first data of both time-series as
    % the system needs to be equilibrated (in steady state).
    % for measurements to be meaningful.

    Tosc=TTT; % For convenience, gene oscillator time
    Tlen=TTT; % For convenience, bacterial length time
    Yosc=YYY(:,1); % gene oscillator time-series
    Ylen=YYY(:,3); % bacterial length time-series

    % Eliminating equilibration data from time-series
    idxosc = find(Tosc>wup);
    TTosc = Tosc(idxosc);
    YYosc = Yosc(idxosc,:);
    idxlen = find(Tlen>wup);
    TTlen = Tlen(idxlen);
    YYlen = Ylen(idxlen);
    % First of all, obtain an uniformly sampled time-series for cell
    % length and the gene oscillator
    interval = max([mean(diff(TTosc))*2 mean(diff(TTlen))*2]);
    epsilon = interval;
    t1 = max([TTosc(1) TTlen(1)]); t2 = min([TTosc(end) TTlen(end)]);
    TT = [t1+epsilon:interval:t2-epsilon]';
    YYoscnew = interp1(TTosc,YYosc,TT);
    YYlennew = interp1(TTlen,YYlen,TT);
    YYosc = YYoscnew;
    YYlen = YYlennew;

    % Second, detect maximums for cell length time series and the
    % gene oscillator time series. Only local maximums whose values
    % are above the threshold defined by 'peakThreshold' can be
    % candidates to be real maximums.
    peakThreshold=3;
    [AllPeaksOSC,AllPeaksLocOSC,FreqOSC]=peakDetectionHasty(TT,YYosc,peakThreshold);

    figure
    plot(TT/3600,YYosc(:,1),'b');
    hold on
    scatter(TT(AllPeaksLocOSC{1})/3600,YYosc(AllPeaksLocOSC{1}),'filled','g');
    hold off

    % Detecting maximums for cell length time-series (analogous as above)
    peakThreshold=0.5;
    [AllPeaksLEN,AllPeaksLocLEN,FreqLEN]=peakDetectionHasty(TT,YYlen,peakThreshold);

    figure
    plot(TT/3600,YYlen(:,1),'b');
    hold on
    scatter(TT(AllPeaksLocLEN{1})/3600,YYlen(AllPeaksLocLEN{1}),'filled','g');
    hold off

    % Detect minimums for the gene oscillator time series
    [AllMinsOSC,AllMinsLocOSC,dummyFreq]=minsDetectionHasty(TT,YYosc);

    % Constructing a new time series: only 1 column (vector index goes from 1
    % to length(TT)). An element is = 0 if there is no OSC MAX nor MIN, +1 if there
    % is an OSC MAXIMUM in that index, -1 if there is an OSC MIN in that index
    Xtremesosc = zeros(1,length(TT));

    dimMAXSosc=size(AllPeaksLocOSC{1});% Fdim(2) contains number of oscillator maximums
    dimMINSosc=size(AllMinsLocOSC{1});% Fdim(2) contains number of oscillator minimums

    if dimMAXSosc >= 1
        Xtremesosc(AllPeaksLocOSC{1})=1;
    end
    if dimMINSosc >= 1
        Xtremesosc(AllMinsLocOSC{1})=-1;
    end
    %
    % DOUBLE CHECKING WE HAVE ALTERNATE OSC MAXS AND MINS (so that we cannot have for
    % instance 2 minimums between maximums)
    %
    %%%%
    figure
    hold on
    plot(TT,YYosc)
    kkmax=find(Xtremesosc==1);
    kkmin=find(Xtremesosc==-1);
    scatter(TT(kkmax),YYosc(kkmax),'c','filled');
    scatter(TT(kkmin),YYosc(kkmin),'m','filled');
    hold off
    dummyXtremesosc=Xtremesosc;
    kkMAX=[AllPeaksLocOSC{1}' ones(length(AllPeaksLocOSC{1}),1)];
    kkMIN=[AllMinsLocOSC{1}' ones(length(AllMinsLocOSC{1}),1).*(-1)];
    kkota=cat(1,kkMAX,kkMIN);
    [~,srtdId]=sort(kkota(:,1),1);
    srtdkkota=kkota(srtdId,:); % "zipper" style (sorted final array containing indices for
    % maxs (+1) and mins (-1) (hence, location in
    % ascending order)
    % Now we need to detect when we have 2 consecutive maxs or mins
    repXtremes=find(diff(srtdkkota(:,2))==0);% indices affected: repXtremes and repXtremes+1

    while repXtremes
        % Once detected we need to act: choose only one of the consecutive maxs or
        % mins
        % Erase the selected repeated max or min from srtdkkota and
        % Xtremesosc arrays
        l=1;
        if srtdkkota(repXtremes(l),2)==1 % Repeated MAX
            if YYosc(srtdkkota(repXtremes(l),1),1) > YYosc(srtdkkota(repXtremes(l),1)+1,1)
                srtdkkota(repXtremes(l),:) = []; % Removing the false max
            else
                srtdkkota(repXtremes(l)+1,:) = [];
            end
        else % Repeated MIN
            if YYosc(srtdkkota(repXtremes(l),1),1) < YYosc(srtdkkota(repXtremes(l),1)+1,1)
                srtdkkota(repXtremes(l),:) = []; % Removing the false min
            else
                srtdkkota(repXtremes(l)+1,:) = [];
            end
        end
        repXtremes=find(diff(srtdkkota(:,2))==0);
    end
    % RECONSTRUCTING final version of Xtremesosc
    finalXtremesosc = zeros(1,length(TT));
    lensrtdkkota=size(srtdkkota);
    for l=1:lensrtdkkota(1)
        finalXtremesosc(srtdkkota(l,1))=srtdkkota(l,2);
    end
    figure
    hold on
    plot(TT,YYosc)
    finalkkmax=find(finalXtremesosc==1);
    finalkkmin=find(finalXtremesosc==-1);
    scatter(TT(finalkkmax),YYosc(finalkkmax),'m ','filled');
    scatter(TT(finalkkmin),YYosc(finalkkmin),'c','filled');
    hold off

    Xtremesosc=finalXtremesosc;
    %%%%

    % Constructing a new time series: only 1 column (vector index goes from 1
    % to length(TT)). An element is = 0 if there is no LEN MAX nor MIN, +1 if there
    % is a LEN MAXIMUM in that index, -1 if there is a LEN MINIMUM in that index
    % (index of LEN mins = index of LEN max +1)
    Xtremeslen = zeros(1,length(TT));

    dimMAXSlen=size(AllPeaksLocLEN{1});% Fdim(2) contains number of cell length maximums
    if dimMAXSlen >= 1
        Xtremeslen(AllPeaksLocLEN{1})=1;
        Xtremeslen(AllPeaksLocLEN{1}+1)=-1;
    end

    % Next, we can compute the Pulse Triggered Averages now (triggered by cell
    % length)
    %
    % SEGMENTING DATA VIA CELL LENGTH MAXIMUMS
    %

    if dimMAXSlen(2) >= 2 % we need at least 2 maximums to be able to construct a segment
        %
        % EXTRACTING DATA FROM SEGMENTS
        %
        % SEGMENTING DATA VIA CELL LENGTH MAXIMUMS
        for i=1:dimMAXSlen(2)-1
            Lsegdata(i).len = YYlen(AllPeaksLocLEN{1}(i):AllPeaksLocLEN{1}(i+1));% From max to max of cell LENGTH
            Lsegdata(i).osc = YYosc(AllPeaksLocLEN{1}(i):AllPeaksLocLEN{1}(i+1));% Idem

            Lsegdata(i).Oxtremes = Xtremesosc(AllPeaksLocLEN{1}(i):AllPeaksLocLEN{1}(i+1));
        end

        % ALIGNING SEGMENTS AT THE MINIMUM OF CELL LENGTH FOR EACH SEGMENT:
        % 1.- Find the minimum of cell length and keep the index (for each segment)
        % 2.- Split the segment into RIGHT data (reference is the minimum of cell length)
        %     and LEFT data (do this for fluorescence time series and cell lengths).
        % 3.- Due to different lengths between segments, we need to compute the
        %     average carefully (hence the use of the 'continue' function).
        for i=1:dimMAXSlen(2)-1
            % The minimum of cell lenght occurs (in PTA triggered by length) in
            % the segment index = 2 (as index = 1 is the first length maximum)
            dummyLidx=2;

            % Remember that Fdimstorage(i) contains the length of Lsegdata(i).mf
            % As segmenting has been done from len max to len max, the minimum
            % value will correspond to the real minimum (cell division time)

            Lsegdata(i).lenalignmin=dummyLidx;% Contains the index 'j' of Lsegdata(i).mf(j) from which
            % the segment should be aligned

            OSC.Loscavrgalmin(i).RIGHTvals=Lsegdata(i).osc(Lsegdata(i).lenalignmin:end);% Stores osc values
            % after (and including) the minimum of the segment (and for each
            % segment 'i')
            OSC.Loscavrgalmin(i).LEFTvals=Lsegdata(i).osc(Lsegdata(i).lenalignmin-1:-1:1);% Stores osc values
            % before (NOT including) the minimum of the segment (and for each
            % segment 'i')
            LEN.Llenavrgalmin(i).RIGHTvals=Lsegdata(i).len(Lsegdata(i).lenalignmin:end);% Idem for length data
            % (aligning length data taking as a reference the minimum of
            % cell length
            LEN.Llenavrgalmin(i).LEFTvals=Lsegdata(i).len(Lsegdata(i).lenalignmin-1:-1:1);% Idem for length data

            XTR.LOxtremealmin(i).RIGHTvals=Lsegdata(i).Oxtremes(Lsegdata(i).lenalignmin:end);
            XTR.LOxtremealmin(i).LEFTvals=Lsegdata(i).Oxtremes(Lsegdata(i).lenalignmin-1:-1:1);
        end

        OSC.Roscdummy=[];% len(Loscavrgalmin(i).RIGHTvals)=len(Llenavrgalmin(i).RIGHTvals)
        OSC.Loscdummy=[];% Idem for LEFTvals
        for i=1:dimMAXSlen(2)-1
            OSC.Roscdim=size(OSC.Loscavrgalmin(i).RIGHTvals');
            OSC.Roscdummy=cat(1,OSC.Roscdummy,OSC.Roscdim(2));% RIGHT
            OSC.Loscdim=size(OSC.Loscavrgalmin(i).LEFTvals');
            OSC.Loscdummy=cat(1,OSC.Loscdummy,OSC.Loscdim(2));% LEFT
        end

    end
    %
    % PLOTTING SEGMENTED RAW DATA
    %
    cc=jet(dimMAXSlen(2)-1); % Generating a colormap, EACH SEGMENT IS COLOURED IN FUNCTION OF
    % IT'S "AGE" (BLUE -> YOUNGEST, RED -> OLDEST)
    figname2=sprintf('%s_PTAcellDivisionrawdata',dataID);% Pulse Triggered Averages
    figtitle1='[Gene oscillator] (a.u.)';
    h2=figure(2);
    subplot(2, 1, 1)
    hold on
    for i=1:dimMAXSlen(2)-1
        plot([fliplr(-1*[1:length(OSC.Loscavrgalmin(i).LEFTvals)]) [0:length(OSC.Loscavrgalmin(i).RIGHTvals)-1]]*interval/60,[OSC.Loscavrgalmin(i).LEFTvals(end:-1:1)' OSC.Loscavrgalmin(i).RIGHTvals'],'color',cc(i,:),'LineStyle',':','LineWidth',2)
        segmentXTREMES=[fliplr(XTR.LOxtremealmin(i).LEFTvals) XTR.LOxtremealmin(i).RIGHTvals];
        kkidxMXS=find(segmentXTREMES==1);
        kksegmentTIME=[fliplr(-1*[1:length(OSC.Loscavrgalmin(i).LEFTvals)]) [0:length(OSC.Loscavrgalmin(i).RIGHTvals)-1]]*interval/60;
        kksegmentDATA=[OSC.Loscavrgalmin(i).LEFTvals(end:-1:1)' OSC.Loscavrgalmin(i).RIGHTvals'];
        scatter(kksegmentTIME(kkidxMXS)',kksegmentDATA(kkidxMXS)','k','filled');
    end
    title(figtitle1,'fontsize',12)
    set(gca,'fontsize',10)

    subplot(2,1,2)
    figtitle2='Length (a.u.)';
    hold on
    for i=1:dimMAXSlen(2)-1
        plot([fliplr(-1*[1:length(LEN.Llenavrgalmin(i).LEFTvals)]) [0:length(LEN.Llenavrgalmin(i).RIGHTvals)-1]]*interval/60,[LEN.Llenavrgalmin(i).LEFTvals(end:-1:1)' LEN.Llenavrgalmin(i).RIGHTvals'],'color',cc(i,:),'LineStyle',':','LineWidth',2)
    end
    xlabel('"Time" (min)','fontsize',16);
    title(figtitle2,'fontsize',12)
    set(gca,'fontsize',10)

    ha2 = axes('Position',[0 0 1 1],'Xlim',[0 1],'Ylim',[0 1],'Box','off',...
        'Visible','off','Units','normalized', 'clipping' , 'off');

    colormap(jet);
    ch = colorbar;
    yl=get(ch,'ylim');
    zl=get(gca,'zlim');
    set(ch,'ytick',yl);
    set(ch,'yticklabel',{'1st';strcat(num2str(dimMAXSlen(2)-1),'th gen')});

    xch=get(ch,'Position');
    xch(2)=0.03;
    xch(2)=0.03;
    xch(4)=0.9;
    set(ch,'Position',xch)

    saveas(h2,figname2,'png')
    saveas(h2,figname2,'pdf')

    %
    % COMPUTING THE HISTOGRAM OF ALL GENE OSCILLATOR MAXIMUMS IN EACH
    % SEGMENT
    %
    storeLoscMAXs=[];
    storeLoscMAXsNormt=[];
    noMAXscounter=0;% Counting the number of segments in which there is no maximum for the gene oscillator
    for i=1:dimMAXSlen(2)-1
        segmentDATA=[fliplr(OSC.Loscavrgalmin(i).LEFTvals') OSC.Loscavrgalmin(i).RIGHTvals'];
        segmentTIME=[fliplr(-1*[1:length(OSC.Loscavrgalmin(i).LEFTvals)]) [0:length(OSC.Loscavrgalmin(i).RIGHTvals)-1]]*interval/60;
        newsegmentTIME=segmentTIME(2:end);

        segmentTIMENormt=newsegmentTIME/max(newsegmentTIME);% t_1=0 AT THE LENGTH MINIMUM, t_2=1 AT
        % THE SECOND LENGTH MAXIMUM

        segmentXTREMES=[fliplr(XTR.LOxtremealmin(i).LEFTvals) XTR.LOxtremealmin(i).RIGHTvals];
        newsegmentXTREMES=segmentXTREMES(2:end);

        locsSEG=find(newsegmentXTREMES==1);
        if ~isempty(locsSEG)
            for k=1:length(locsSEG)
                storeLoscMAXsNormt=cat(1,storeLoscMAXsNormt,segmentTIMENormt(locsSEG(k)));
            end
        else
            noMAXscounter=noMAXscounter+1;
        end
    end

    % PLOTTING THE HISTOGRAM
    Roscmax=max(OSC.Roscdummy);% Maximum length for RIGHT vectors
    Loscmax=max(OSC.Loscdummy);% Idem for LEFT

    newbin=0:0.08:1+0.08;%0:0.05:1;%0:0.1:1; % 0:1/7.5:1+1/7.5;
    Losc_nNormt=histc(storeLoscMAXsNormt,newbin);

    figname1000012=sprintf('%s_PTAcellDivision_OscMAXShistogramNormt',dataID);% Pulse Triggered Averages;
    h1000012=figure(1000012);
    figtitle1000012 = sprintf('%s PTA by lenght, distr. of all GENE OSC MAX, N=%s',dataID,num2str(sum(Losc_nNormt)));
    bar(newbin,Losc_nNormt/sum(Losc_nNormt),'histc')
    title(figtitle1000012,'fontsize',12)
    xlabel('Norm. time (a.u.)','fontsize',12);
    set(gca,'fontsize',12)

    saveas(h1000012,figname1000012,'png')
    saveas(h1000012,figname1000012,'pdf')
    %
    % Rearranging the previous bar plot (PTA):
    % For the sake of clarity, phase has been redefined so that 0.5 corresponds to the
    % moment when cell achieves its maximum length (right before division).
    % In this way, it is easier to compare the three cases: No coupling, unidirectional
    % coupling, and bidirectional coupling.
    [mindist,minindex] = locate(0.5,newbin);
    if newbin(end) > 1
        kkota=newbin(1:end-1)-newbin(minindex);
        finalLoscNormt_n=Losc_nNormt(1:end-1);
    else
        kkota=newbin-newbin(minindex);
        finalLoscNormt_n=Losc_nNormt;
    end

    kk=find(kkota<0);
    kk3=find(kkota>=0);
    kkota3=kkota(kk3);
    kkota2=kkota(kk)+max(abs(kkota(kk)))+(kkota3(2)-kkota3(1))*length(kkota3);
    kkotaALL=[kkota2 kkota3];

    figname1000013=sprintf('%s_PTAcellDivision_OscMAXShistogramNormt_rearranged',dataID);
    h1000013=figure(1000013);
    figtitle1000013 = sprintf('%s PTA by Lenght, distr. of all GENE OSC MAX, N=%s, re-arranged',dataID,num2str(sum(finalLoscNormt_n)));
    %figtitle1000013 = sprintf('%s PTA by Lenght, distr. of all GENE OSC MAX, N=%s+%s, re-arranged',dataID,num2str(sum(finalLoscNormt_n)),num2str(noMAXscounter));
    bar(kkotaALL,finalLoscNormt_n/sum(finalLoscNormt_n),'histc');
    title(figtitle1000013,'fontsize',12)
    ylim([0 0.35]);
    xlabel('Norm. Time (a.u.)','fontsize',12);
    set(gca,'fontsize',12)

    saveas(h1000013,figname1000013,'png')
    saveas(h1000013,figname1000013,'pdf')

    % COMPUTING PERIODS DISTRIBUTIONS (FOR CELL LENGTH AND GENE OSCILLATOR)
    % Maximums for fluorescence
    storeOSCperiods=[];
    oscmaxsok=find(Xtremesosc==1);
    dimMAXSosc=size(oscmaxsok);
    if dimMAXSosc(2) >= 2
        for i=1:dimMAXSosc(2)-1
            dummyperiod=(TT(oscmaxsok(i+1))-TT(oscmaxsok(i)))/60;
            storeOSCperiods=cat(1,storeOSCperiods,dummyperiod);
        end
    end
    % Computing the corresponding histogram
    bins3=min(storeOSCperiods):5:max(storeOSCperiods);% 5
    OSCperiods_n=histc(storeOSCperiods,bins3);

    figname100003=sprintf('%s_OSCperiodshistogram',dataID);
    h100003=figure(100003);
    figtitle100003 = sprintf('OSC periods distribution, N=%s',num2str(sum(OSCperiods_n)));
    bar(bins3,OSCperiods_n/sum(OSCperiods_n),'histc')
    title(figtitle100003,'fontsize',12)
    xlabel('Periods (min)','fontsize',12);
    set(gca,'fontsize',12)

    saveas(h100003,figname100003,'png')
    saveas(h100003,figname100003,'pdf')

    disp('   ');
    disp(dataID);
    disp('   ');
    OSCdevResults=sprintf('OSC period: mean +- std = %s +- %s (in min)',num2str(mean(storeOSCperiods)), num2str(std(storeOSCperiods)));
    disp(OSCdevResults);
    disp('   ');

    % Maximums for cell length
    storeLENperiods=[];
    if dimMAXSlen(2) >= 2
        for i=1:dimMAXSlen(2)-1
            dummyperiod=(TT(AllPeaksLocLEN{1}(i+1))-TT(AllPeaksLocLEN{1}(i)))/60;
            storeLENperiods=cat(1,storeLENperiods,dummyperiod);
        end
    end
    % Computing the corresponding histogram
    bins4=min(storeLENperiods):5:max(storeLENperiods);% 1
    LENperiods_n=histc(storeLENperiods,bins4);

    figname100004=sprintf('%s_LENperiodshistogram',dataID);
    h100004=figure(100004);
    figtitle100004 = sprintf('LEN periods distribution, N=%s',num2str(sum(LENperiods_n)));
    bar(bins4,LENperiods_n/sum(LENperiods_n),'histc')
    xlim([20 220])
    title(figtitle100004,'fontsize',12)
    xlabel('Periods (min)','fontsize',12);
    set(gca,'fontsize',12)

    saveas(h100004,figname100004,'png')
    saveas(h100004,figname100004,'pdf')

    LENdevResults=sprintf('LEN period: mean +- std = %s +- %s (in min)',num2str(mean(storeLENperiods)), num2str(std(storeLENperiods)));
    disp(LENdevResults);
    disp('   ');

end
