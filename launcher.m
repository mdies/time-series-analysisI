% Launching time-series data processing for 2 different data sets:
% 1.- Unidirectional coupling (gene oscillator is coupled to bacterial cell
%                  cycle).
% 2.- Bidirectional coupling (gene oscillator is coupled to bacterial cell
%                  cycle AND, in turn, bacterial cell cycle is 'weakly'
%                  coupled to gene oscillator).
% See M.Dies, L.Galera-Laporta, and J.Garcia-Ojalvo, Integrative Biology
% (2016) for further details.
listOfSims={'unidirectionalCouplingTimeSeries';...
    'bidirectionalCouplingTimeSeries'};

for i=1:length(listOfSims)
% Simulation data contains several arrays from which we need: 
% TTT: time (in seconds)
% YYY: time-series for the different species/observables (YYY)
%      Specifically, YYY(:,1) corresponds to the gene oscillator reporter (in a.u.)
%      and YYY(:,3) corresponds to bacterial cell length (in a.u.)
    load(listOfSims{i});
    dataID = erase(listOfSims{i},"TimeSeries");
    processingSimData_Portfolio(TTT,YYY,dataID);
    close all;
end