#iSite Tools
All files should be copied into the same directory as the iSite executable.
##Queueing jobs
###iSite.pbs
Customize PBS directives for your machine.
###iSite.py
 Customize parameters then run this file to begin queueing jobs.

##Collecting Data
I can't guarantee that changing any parameters in these files won't break anything. Meddle at your own risk (Or make them more generic and submit a pull request!).

###average.rb
Averages a single parameter set's runs and places the averages in a single file.
###average\_dist.rb
Averages a single parameter set's domain distribution and places the averages in a single file.
###collect.rb
Collects averages from each parameter set and places them in a single file. Should be sorted by fusion/fission values due to naming conventions.
###collect\_data.sh
Run this to aggregate results. This really just runs the above tools multiple times. End results will be in `results/result`.
