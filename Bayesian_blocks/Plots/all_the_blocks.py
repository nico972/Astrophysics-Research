import astropy.io.fits
from astropy.io import fits
import matplotlib.pyplot as plt
import numpy as np
import array


       
ObsID = ['14702',
 '15040',
 '14703',
 '15651',
 '15654',
 '14946',
 '15041',
 '15042',
 '14945',
 '15043',
 '14944',
 '15044',
 '14943',
 '14704',
 '15045',
 '16508',
 '16211',
 '16212',
 '16213',
 '16214',
 '16210',
 '16597',
 '16215',
 '16216',
 '16217',
 '16218']

plot_step = 300

files = []
files2 = []
i=0
for i in range(0,len(ObsID)):
    files.append(ObsID[i] + "_sgra_events_2-8keV_1.25asec.fits")

i=0
for i in range(0,len(ObsID)):
    files2.append(ObsID[i] + "_sgra_bkgsub_lc_2-8keV_1.25asec.fits")
    
ffile = []
data = []
time=[]
vertical = []
time_starts  = []
tstart = []
tstop = []

ffile2 = []
data2 = []

i=0
for i in range(0,len(files)):
    ffile.append(fits.open(files[i]))
    data.append(ffile[i][1].data)
    tstart.append(ffile[i][1].header['TSTART'])
    tstop.append(ffile[i][1].header['TSTOP'])

i=0
for i in range(0,len(files)):
    ffile2.append(fits.open(files2[i]))
    data2.append(ffile2[i][1].data)

order = np.argsort(tstart)

real = []

i=0
for i in range(0,len(ObsID)):
    real.append(ObsID[i])

output = []
blocks = []
rates = []
bining = []
bstart = []
bend = []
time = []
time_true = []
step=[]
time_stops = []

y = []
err = []
vertical = []

import pwkit
from pwkit import bblocks

i=0
for i in range(0,len(ObsID)):
    time_true.append(data[i].field('TIME'))
    if i > 0:
        time.append(time_true[i]-tstart[i])
        vertical.append(time[i][-1]+vertical[i-1])
    else:
        time.append(time_true[i]-tstart[i])
        vertical.append(time[i][-1])
    
    y.extend(data2[i].field('COUNT_RATE'))
    err.extend(data2[i].field('COUNT_RATE_ERR'))

cool = 300 *(len(y))

timing = np.arange(150,cool,300)

i=0
j=0


for i in range(0,len(ObsID)):
    end = tstop[i] - tstart[i]
    time_starts.append(np.arange(0,end,plot_step))
    step.append([300]*len(time_starts[i]))
    time_stops.append(time_starts[i] + step[i])
    output.append(bblocks.tt_bblock(time_starts[i], time_stops[i], time[i], p0=0.05))

'''
i=0
for i in range(0,len(ObsID)):
    blocks.append(output[i].blockstarts)
    rates.append(output[i].rates)
    bining.append(output[i].widths)
    bstart.append(output[i].ledges)
    bend.append(output[i].redges)
'''
i=0
for i in range(0,len(ObsID)):
    if i > 0:
        blocks.extend(output[i].blockstarts+blocks[-1])
        rates.extend(output[i].rates)
        bining.extend(output[i].widths)
        bstart.extend(output[i].ledges+bend[-1])
        bend.extend(output[i].redges+bend[-1])
    else:
        blocks.extend(output[i].blockstarts)
        rates.extend(output[i].rates)
        bining.extend(output[i].widths)
        bstart.extend(output[i].ledges)
        bend.extend(output[i].redges)
        
#plt.errorbar(timing, y,err,ecolor='g',alpha=0.1)

j=0

    
'''
for j in range(0,len(rates)):
    if j < len(rates)-1:
        plt.bar(bstart[j]-bstart[0],rates[j], width=(bstart[j+1]-bstart[j]),
                alpha=1,color='r')
    else:
       plt.bar(bstart[j]-bstart[0],rates[j], width=(bend[j]-bstart[j]),
                alpha=1,color='r')

'''

flare_block_number = [ [2,3] , ['nan'], ['nan'], [2]]

'''
flares=[]
flares.append(1500+3062.2011278/2)
flares.append(14086.69902754 + 2413.30097246/2)
flares.append(7382.38815305 + 1041.72473171/2)
flares.append(34500 + 833.64574081/2)
'''
fig, ax = plt.subplots(2,1,sharex=True, sharey=False)
#fig.tight_layout()
ax[0].errorbar(timing, y,err,ecolor='g')
#ax[0].vlines(flares,0,0.6,colors='r',linestyles='dashed')
ax[0].vlines(vertical,0,0.6,colors='r',linestyles='dashed')

j=0
for j in range(0,len(rates)):
    if j < len(rates)-1:
        ax[1].bar(bstart[j]-bstart[0],rates[j], width=(bstart[j+1]-bstart[j]),
                alpha=1,color='r')
    else:
       ax[1].bar(bstart[j]-bstart[0],rates[j], width=(bend[j]-bstart[j]),
                alpha=1,color='r')

plt.xlabel("Time(s)")
for ax in ax:
    ax.set_ylabel('Count rates (count/s)')

                
plt.show()

    