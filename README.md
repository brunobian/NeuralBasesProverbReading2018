# Data and analysis code from Bianchi et al. (Frontiers in Human Neuroscience, under review)
These are the code and data to build the figures from Bianchi, Shalom and Kamienkowski (under review). 

The code is written almost entirely in Matlab, and only uses some functions from [FieldTrip] for Fig. S2B. Pre-analysis (not included) also uses [EEGLAB]. The code for the LMM-CBPT (https://github.com/brunobian/LMM-CBP) combines Matlab and R functions. In particular, it uses lme4 for Linear Mixed Models. For further details please refer to Bianchi et al. (under review).

This work is part of a larger [project](http://reading.liaa.dc.uba.ar/) aimed to understand cognitive and neural processes of prediction during reading. We combined eye movement, EEG and online experiments with computational modelling in order to disentangle the different sources of predictibility using natural stimuli, such as proverbs or short stories. We are part of the [Applied Artificial Intelligence Laboratory][LIAA] @[Computer Science Department][DC], [School of Exact and Natural Science][FCEyN], [University of Buenos Aires][UBA], Argentina, and the [Neuroscience Lab][LNUTDT] @Torcuato Di Tella University, Buenos Aires, Argentina.

# How to cite us
#### Please, if you like it / use it cite us:
Bianchi B, Shalom DE, and Kamienkowski JE, “Predicting known sentences: neural basis of proverb reading using nonparametric statistical testing and mixed-effects models” (under review)
#### And let us know!!
Bruno Bianchi (bbianchi (arroba) dc (dot) uba (dot) ar)
Diego E. Shalom (diegoshalom (arroba) gmail (dot) com)
Juan E. Kamienkowski (juank (arroba) dc (dot) com (dot) ar)

# Required toolboxes
* EEGLAB (version XXXX): https://sccn.ucsd.edu/eeglab/index.php
* Fieldtrip (version XXXX): http://www.fieldtriptoolbox.org/
* LMM-CBP (version XXXX): https://github.com/brunobian/LMM-CBP

EEGLAB and Fieldtrip run entirely over Matlab, LMM-CBP run from Matlab and but it uses R (VERSION) and lme4 package (VERSION, cita) 

# Data
#### exp.mat
Data from stimuli used on the experiment. This file contains a MatLab structure array with 120 records, one per sentence. For further details please refer to Bianchi et al. (under review). 

#### EEG.mat
This file contains 3 matlab structures: erp, win, CHANS:
* CHANS contains information from the biosemi cap used for the EEG recording. 
* win contains information about windows used for averaging within epochs (i.e. time limits and electrodes). 
* erp contains all the information from the EEG preanalysis. 
For further details please refer to Bianchi et al. (under review).

#### EEG_remef.csv
For Fig 5, the 'Position*Sentences type' partial effect was removed in the R environment, using the R function remef() (Hohenstein & Kliegl, 2013; R package version 1.0.6.9000, https://github.com/hohenstein/remef) (see run_remef.R). ERP amplitudes were then calculated after removing this effect and exported into a csv. This csv is loaded in MatLab and analized in figures/fig5.m. For further details please refer to Bianchi et al. (under review).

# Figures
figures/\*.m contains all the files with the codes to create all the Figures from Bianchi et al. (under review). They uses EEGLAB, Fieltrip and LMM-CBP (See 'required packages') toolboxes and niceBars2.m (also included). Thus, it's necessary to add fuctions/ folder to the path, as well as those packages depending on which figure you  want to run. It's also necessary to add data/ folder to the path or include in the same folder the necessary summary data to create the figures,





