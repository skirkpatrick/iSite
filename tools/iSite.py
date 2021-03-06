#!/usr/bin/env python
import subprocess
import time
SEEDGRAPH='newseed.r100.83self'
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
DISTFILE='dist'
SRCDIR='/home/skirkpatrick/workspace/iSite'
iteration=1
#0%-20% inclusive, step by 1%
for PROBFUSION in [x*.01 for x in range(4, 15)]:
    for PROBFISSION in [y*.01 for y in range(11)]:
        args = 'SEEDGRAPH=' + SEEDGRAPH +                           \
               ',PROBSUBFUNC=' + PROBSUBFUNC +                      \
               ',PROBASYM=' + PROBASYM +                            \
               ',PROBHOMO=' + PROBHOMO +                            \
               ',PROBFUSION=' + '{:.2f}'.format(PROBFUSION) +       \
               ',PROBFISSION=' + '{:.3f}'.format(PROBFISSION) +     \
               ',ENDORDER=' + ENDORDER +                            \
               ',ITERATIONS=' + ITERATIONS +                        \
               ',OUTDIR=' + OUTDIR + '.{:.2f}.{:.3f}'.format(PROBFUSION, PROBFISSION) + \
               ',OUTFILE=' + OUTFILE + '.{:d}'.format(iteration) +  \
               ',PRINTRESULTS=' + PRINTRESULTS +                    \
               ',GRAPHFILE=' + GRAPHFILE +                          \
               ',DISTFILE=' + DISTFILE +                            \
               ',SRCDIR=' + SRCDIR
        subprocess.call('qsub -v ' + args + ' iSite.pbs', shell=True)
user = subprocess.check_output('id -u -n', shell=True).decode().strip()
num_jobs_cmd = 'qstat | grep ' + user + ' | grep -v \' C \' | wc -l'
while True:
    num_jobs = subprocess.check_output(num_jobs_cmd, shell=True).decode().strip()
    if (int(num_jobs) < 1):
        break
    print('Remaining jobs: ' + num_jobs)
    time.sleep(10)
