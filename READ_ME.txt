Title: Larger anticipatory postural adjustments relate to later and larger reactive steps in people with PD
READ_ME File by: Keith Lohse
Date: 2018-05-01

In this folder you will find a script file for recreating the statistical analyses for the paper, "Larger 
anticipatory postural adjustments relate to later and larger reactive steps in people with PD" by Dan 
Peterson, Keith Lohse, and Martina Macini. 

The folder also includes two data files, "data_PETERSON_MASTER_050118.csv" and "data_PETERSON_MASTER_no15.csv". 
The ...050118.csv can be used for the initial analyses, whereas the ...no15.csv can be used for the final 
analyses. The difference between these two data files is the exclusion of participant #15, who was found to be 
overly influential in several analyses based on the Cook's distance statistic. 

Within each data file, however, the constinuent variables are the same:
"subject" = anonymous numerical ID for each participant, we convert this to a factor in the script
"direction" = categorical variable for the direction of perturbation
"trial" = the relative trial number in each direction (1-25)
"trial_num" = the absolute trial number across both directions (1-50)
"FOGstatus" = a numerical variable indicated if the participant experience freezing of gait (1) or not (-1)
"COMdisp" = center of mass displacement	
"AP_MOS" = anterior-posterior margin of stability	
"ML_MOS" = mediolateral margin of stability	
"stepWidth" = step width for the protective step in meters	
"StepLength" = step length for the protective step in meters	
"AP_COMatFO" = anterior-posterior center of mass at foot off	
"ML_COMatFO" = mediolateral center of mass at foot off	
"Sla" = step latency in seconds	
"APA" = anticipatory postural adjustment in meters	
"COP_pos" = center of pressure position at the perturbation	
"StepLeg" = categorical variable indicating the stepping leg (left of right)
