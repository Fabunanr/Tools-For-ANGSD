# Tools-For-ANGSD

These tools will help visualize the pestPG output files from ANGSD. One visualization tool uses the shiny package from R and I just modified it from @arundurvasula (https://github.com/arundurvasula/angsd-wrapper/wiki/Shiny-Graphing). The other is a python script that will help pull out regions with specific thetaW and Tajima's D values from two different pestPG files.

## Shiny visualization

This is a modified version of @arundurvasula's script that now allows three different files to be visualized in the main panel. 
In this case, each graph is titled as either high, middle, or low. This is specific for my use.

### Requirements
R

#### R packages:

Shiny

genomeIntervals

## Python script for specific thetaW and Tajima's D

This python script allows you to obtain regions from the pestPG file that have specific Tajima's D values in two different pestPG files. The custumization is easy as it will ask which files to use, chromosome, above or below tajima's D value, above or below thetaW value.

### Requirements
any version of python.

#### How to use:
place the script (.py file) into a directory that contains the two pestPG files needed. Then type "python script.py" and the terminal should guide you through the rest of the steps.
