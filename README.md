# Data and analysis code from Bianchi et al. (Frontiers in Human Neuroscience, under review)
The code for novel EEG analysis technique to study multiple co-varibles over the whole spatio-temporal ERP signal can be found at https://github.com/brunobian/LMM-CBP. It combines Linear-Mixed Models and Cluster-Based Permutation procedures into a single framework, more details can be found in Bianchi, Shalom and Kamienkowski (under review).

These are the code and data to build the figures from Bianchi, Shalom and Kamienkowski (under review). 

The code is written almost entirely in Matlab, and only uses some functions from [FieldTrip] for Fig. S2B. Pre-analysis (not included) also uses [EEGLAB]. The code for the [LMM-CBP](https://github.com/brunobian/LMM-CBP) combines Matlab and R functions. In particular, it uses lme4 for Linear Mixed Models. For further details please refer to Bianchi et al. (under review).

This work is part of a larger [project](http://reading.liaa.dc.uba.ar/) aimed to understand the cognitive and neural processes of prediction during reading. We combined eye movement, EEG and online experiments with computational modelling in order to disentangle the different sources of predictibility using natural stimuli, such as proverbs or short stories. We are part of the [Applied Artificial Intelligence Laboratory](http://liaa.dc.uba.ar/) [LIAA] @[Computer Science Department](http://dc.uba.ar/) [DC], [School of Exact and Natural Science](http://exactas.uba.ar) [FCEyN], [University of Buenos Aires](http://www.uba.ar) [UBA], Argentina, and the [Neuroscience Lab](https://www.utdt.edu/ver_contenido.php?id_contenido=10518&id_item_menu=20132) [LNUTDT] @Torcuato Di Tella University, Buenos Aires, Argentina.


# How to cite us
#### Please, if you like it / use it cite us:
Bianchi B, Shalom DE, and Kamienkowski JE, “Predicting known sentences: neural basis of proverb reading using nonparametric statistical testing and mixed-effects models” (under review)
#### And let us know!!
* Bruno Bianchi (bbianchi (at) dc (dot) uba (dot) ar)
* Diego E. Shalom (diegoshalom (at) gmail (dot) com)
* Juan E. Kamienkowski (juank (at) dc (dot) com (dot) ar)

# Required toolboxes
* [MatLab] EEGLAB (version v14.1.1): https://sccn.ucsd.edu/eeglab/index.php 
* [MatLab] Fieldtrip (version 3e7ad536c, 20170827): http://www.fieldtriptoolbox.org/ 
* [R] lm4 (version 1.1-18-1): https://cran.r-project.org/web/packages/lme4/lme4.pdf 
* [MatLab, Bash, R] LMM-CBP (version 0.1): https://github.com/brunobian/LMM-CBP 

EEGLAB and Fieldtrip run entirely over Matlab. LMM-CBP also run from Matlab, but it connects with R (version 3.4.4) and lme4 package (version 1.1-18-1, Bates, et al. (2015) “Fitting linear mixed-effects models using lme4” Journal of Statistical Software 67, 1–48). This conection is made through bash scripts.

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
For Fig 5, the 'Position*Sentences type' partial effect was removed in the R environment, using the R function remef() (Hohenstein & Kliegl (2013) “remef: Remove partial effects” R package version 1.0.6.9000, https://github.com/hohenstein/remef) (see run_remef.R). ERP amplitudes were then calculated after removing this effect and exported into a csv. This csv is loaded in MatLab and analyzed in figures/fig5.m. For further details please refer to Bianchi et al. (under review).

#### LMM-CBP data (not included in this repository)
The data to run the LMM-CBP is not included in this repository because it weights ~4.9Gb. Please, find it in our project website: http://reading.liaa.dc.uba.ar/ 

# Figures' code
figures/\*.m contains all the files with the codes to create all the Figures from Bianchi et al. (under review). They use EEGLAB, FielTrip and LMM-CBP (See 'required packages') toolboxes and niceBars2.m (also included). Thus, it's necessary to add functions/ folder to the path, as well as those packages depending on which figure you want to run. It's also necessary to add data/ folder to the path or include in the same folder the necessary summary data to create the figures,
