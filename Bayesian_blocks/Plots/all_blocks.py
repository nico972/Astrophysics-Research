import astropy.io.fits
from astropy.io import fits
import matplotlib.pyplot as plt
import numpy as np
import array

ObsID = ["14702","14703","14704","14943","14944"
       ,"14945","14946","15040","15041","15042","15043","15044",
       "15045","15651","15654","16210","16211","16212","16213",
       "16214","16215","16216","16217","16218","16508","16597"]
       
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
i=0
for i in range(0,len(ObsID)):
    files.append(ObsID[i] + "_sgra_events_2-8keV_1.25asec.fits")
    
ffile = []
data = []
time=[]
vertical = []
time_starts  = []
tstart = []
tstop = []

i=0
for i in range(0,len(files)):
    ffile.append(fits.open(files[i]))
    data.append(ffile[i][1].data)
    tstart.append(ffile[i][1].header['TSTART'])
    tstop.append(ffile[i][1].header['TSTOP'])
    
order = np.argsort(tstart)

real = []

for i in range(0,len(ObsID)):
    real.append(ObsID[order[i]])

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

import pwkit
from pwkit import bblocks

i=0
for i in range(0,len(order)):
    time_true.append(data[order[i]].field('TIME'))
    if order[i] > 0:
        time.append(time_true[i]-tstart[order[i]])
    else:
        time.append(time_true[i]-tstart[order[i]])


i=0
j=0

 
i = 0 
j=0

for i in range(0,len(ObsID)):
    end = tstop[order[i]] - tstart[order[i]]
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
        bstart.extend(output[i].ledges+bstart[-1])
        bend.extend(output[i].redges+bend[-1])
    else:
        blocks.extend(output[i].blockstarts)
        rates.extend(output[i].rates)
        bining.extend(output[i].widths)
        bstart.extend(output[i].ledges)
        bend.extend(output[i].redges)
        

j=0

    

for j in range(0,len(rates)):
    if j < len(rates)-1:
        plt.bar(bstart[j]-bstart[0],rates[j], width=(bstart[j+1]-bstart[j]),
                alpha=0.5,color='r')
    else:
       plt.bar(bstart[j]-bstart[0],rates[j], width=(bend[j]-bstart[j]),
                alpha=0.5,color='r')
                
plt.show()

    