# Master_Thesis
Code for the pipeline described in the Master Thesis

Flow of the pipeline: 

1. script_23_12_20.m: For every patient, filter and interval: 

	1. Get a table for ripples and fast ripples with the columns: (1) Snippet ID; (2) channel; (3) location index; (4) # of reconstrctions; (5) index of the atom; (6) coefficients of the atom; (7) reconstrction error; (8) v_factor; (9) inclusion AMP; (10) inclusion REC; (11) inclusion REC; (12) threshold REC
  
	2. Get all the spikes: (1) index and (2) channel


2. Get_countsAbdatoms.m: For every patient and detector version: stores the counts of detected EoI and the atom features for every interval and channel 



3. Get_classifications.m: For every patient and detector version: get the predicted channels, the temporal consistency across intervals, the classification
and the rates



4. Calculate_performance.m: Calculate the performance (sensitivity, specificity, NPV, PPV, accuarcy) for every detector version using the classifications 
of all the patients 



Some parts of the code are adapted from:
Besheli, B. F. et al. A sparse representation strategy to eliminate pseudo-HFO events from intracranial EEG for seizure onset zone localization. J Neural Eng 19, 10.1088/1741-2552/ac8766 (2022).

