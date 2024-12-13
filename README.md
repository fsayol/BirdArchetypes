# README file 

In order to run analysis, you first need to run the “0_Setup.R” script.

You can modify the folder where you run to run all the analysis. In such folder, you just need to have the “data” folder with the raw data (Data S1-S4).

The “figures” and the “output” folders, if not created before, will be generated with this initial script, inside the specified repository folder (“RepoFolder”). A “data_tmp” will also be automatically created inside the “data” folder to store some temporary data generated during the analysis.

All the other scripts should be run in consecutive order, as some of them are depending on the results of previous scripts, which will be stored in the output or temporal data folders (data_tmp).
