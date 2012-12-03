#!/usr/bin/env python
import subprocess
import time
SEEDGRAPH='newseed.r100'
PROBSUBFUNC='.82'
PROBASYM='.5'
PROBHOMO='.539'
ENDORDER='2647'
#Number of iterations for single program instance
ITERATIONS='100'
OUTDIR='results'
OUTFILE='result'
PRINTRESULTS='noPrintResult'
GRAPHFILE='g'
SRCDIR='/home/skirkpatrick/workspace/iSite'
#0%-20% inclusive, step by 1%
for PROBFUSION in [x*.01 for x in range(21)]:
    for PROBFISSION in [y*.01 for y in range(21)]:
        args = 'SEEDGRAPH=' + SEEDGRAPH +                           \
               ',PROBSUBFUNC=' + PROBSUBFUNC +                      \
               ',PROBASYM=' + PROBASYM +                            \
               ',PROBHOMO=' + PROBHOMO +                            \
               ',PROBFUSION=' + '{:.2f}'.format(PROBFUSION) +       \
               ',PROBFISSON=' + '{:.2f}'.format(PROBFISSION) +      \
               ',ENDORDER=' + ENDORDER +                            \
               ',ITERATIONS=' + ITERATIONS +                        \
               ',OUTDIR=' + OUTDIR + '.{:.2f}.{:.2f}'.format(PROBFUSION, PROBFISSION) + \
               ',OUTFILE=' + OUTFILE + '.{:d}'.format(iteration) +  \
               ',PRINTRESULTS=' + PRINTRESULTS +                    \
               ',GRAPHFILE=' + GRAPHFILE +                          \
               ',SRCDIR=' + SRCDIR
        subprocess.call('qsub -v ' + args + ' iSite.pbs', shell=True)
user = subprocess.check_output('id -u -n', shell=True).decode().strip()
num_jobs_cmd = 'qstat | grep ' + user + ' | grep -v \' C \' | wc -l'
while True:
    num_jobs = subprocess.check_output(num_jobs_cmd, shell=True).decode().strip()
    if (int(num_jobs) < 100):
        break
    print('Remaining jobs: ' + num_jobs)
    time.sleep(10)
