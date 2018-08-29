# Time series analysis (I)

In order to study the interaction between a natural oscillator, such as bacterial cell cycle, and a [synthetic
one](https://www.nature.com/articles/nature07389), I propose the following time series analysis. 
Using simulated data for both oscillators (see [M.Dies, L.Galera-Laporta, and J.Garcia-Ojalvo, _Integrative Biology_, 2015](http://pubs.rsc.org/en/content/articlelanding/2016/ib/c5ib00262a/unauth#!divAbstract) for details), I 
want to characterize the entrainment between these two oscillators and also quantify their individual periods.

Codes developed in here detect local maximums 
for both time-series, compute their periodicity, and perform a _Pulse triggered average_ (PTA) analysis
in order to establish the degree of entraiment between these two oscillators. This PTA analysis consists in the following steps (see [M.Dies, L.Galera-Laporta, and J.Garcia-Ojalvo, _Integrative Biology_, 2015](http://pubs.rsc.org/en/content/articlelanding/2016/ib/c5ib00262a/unauth#!divAbstract) for further details).
First I defined a phase that accounted for the progress of the system through a cycle, and correspondingly assigned it to each point of the time series data of the two oscillators. I defined a cycle as the segment of data going from one minimum of cell length to the following maximum (thus spanning the entire cell life), and fixed this phase to be 0 at the beginning of the cycle (when the cell is born) and 1 at the end of the cycle (just before the cell divides). Gene oscillator time series was also segmented according to cell length cycles, and phases were assigned correspondingly. Next, I proceed to quantify the existing phase shift between the two oscillators by asking when gene oscillator maxima occur within a cell length cycle.

![PhaseAssignment](https://github.com/mdies/time-series-analysisI/blob/master/PhaseAssignment.pdf)

For the sake of clarity, in PTA final histograms phase was redefined so that 0.5 corresponds to the moment when cell achieves its maximum length (right before division). The data shows that only in the bidirectional case the distribution of gene oscillator maxima is clearly unimodal (and centred around the division time).

## Installation / execution
* Clone this repo to your computer.
* Open MATLAB and cd into `time-series-analysisI` directory.
* To run the program execute `launcher` in the MATLAB prompt.

The program will analyze the time series for two different cases:
* Unidirectional coupling: gene oscillator is coupled to bacterial cell cycle.
* Bidirectional coupling: gene oscillator is coupled to bacterial cell cycle and, in turn, bacterial 
cell cycle is 'weakly' coupled to gene oscillator.

## Results
After execution, for each one of the two scenarios the program returns on the MATLAB prompt 
the computed periods
for gene oscillator and for bacterial cell cycle, and will create the following graphs:

  - distribution of bacterial length periods (uni/bidirectionalCoupling\_LENperiodshistogram.pdf/png)
  - distribution of gene oscillator periods (uni/bidirectionalCoupling\_OSCperiodshistogram.pdf/png)
  - "Pulse triggered average" (PTA raw data) (uni/bidirectionalCoupling\_PTAcellDivisionrawdata.pdf/png)
    (each segment is colored in function of its "age" (blue corresponds to younger generations,
    while red corresponds to older generations)).
  - PTA (uni/bidirectionalCoupling\_PTAcellDivision\_OscMAXShistogramNormt.pdf/png)
  - PTA (uni/bidirectionalCoupling\_PTAcellDivision\_OscMAXShistogramNormt\_rearranged.pdf/png)

Rearranged PTA graphs (for which phase has been redefined) show that only the bidirectional coupling 
leads to a clear unimodal peaked distribution,
meaning that these two oscillators become significantly co-entrained only in the bidirectional case.
The developed PTA and phase asignment methods for processing the time-series in this problem, allow us 
to gain meaningful insight into data.

