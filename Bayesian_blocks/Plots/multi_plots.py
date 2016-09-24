import astropy.io.fits
from astropy.io import fits
import matplotlib.pyplot as plt 
import numpy as np

ObsID = ["14702","14703","14704","14941","14942","14943","14944"
       ,"14945","14946","15040","15041","15042","15043","15044",
       "15045","15651","15654","16210","16211","16212","16213",
       "16214","16215","16216","16217","16218","16508","16597"]
       
       

#For the first 7 
#ObsID = ["14702","14703","14704","14941","14942","14943","14944"
 #      ,"14945","14946","15040","15041","15042"]

#This is the original bining of the data. 
plot_step = 300 

files = []
i=0
for i in range(0,len(ObsID)):
    files.append(ObsID[i] + "_sgra_bkgsub_lc_2-8keV_1.25asec.fits")

ffile = []
data = []

time=[]
time1 = []
time2 = []
y = []
err = []
tstart = []
tstop = []
start = []
end = []
vertical = []


'''
WHat I did here is that I deleted the gaps between the observations
and I stacked each observation next to each other.
'''

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
#delete the gaps

import pwkit
from pwkit import bblocks

i=0
for i in range(0,len(order)):
    if i > 0: 
        time.extend(data[order[i]].field('TIME')-tstart[order[i]]+time[-1])
        vertical.append(time[-1])
    else: 
        time.extend(data[order[i]].field('TIME')-tstart[order[i]])
    
    #start.append(data[i][0].field('TIME')-tstart[0])
    y.extend(data[order[i]].field('COUNT_RATE'))
    err.extend(data[order[i]].field('COUNT_RATE_ERR'))

    end = tstart[order[-1]] - tstart[order[0]]

    time_starts = np.arange(0,end,plot_step)

    array_step = []
    for i in range(0,len(time_starts)):
        array_step.append(plot_step)
    
    time_stops = time_starts + array_step 

    output = bblocks.tt_bblock(time_starts, time_stops, time, p0=0.05)





#Most imoprtant outputs. 
blocks = output.blockstarts
rates = output.rates
bining = output.widths
output.final = output.finalp0

plt.xlabel('time (s)')
plt.ylabel('COUNT_RATE (counts/s)')
plt.title('All the observations')

bstart = output.ledges
bend = output.redges

#Plot every single bar with a width being the difference between
#...the present block and the start of the next block
j=0

for j in range(0,len(rates)-1):
    plt.bar(bstart[j]-bstart[0],rates[j], width=(bstart[j+1]-bstart[j]),
                alpha=0.5,color='r')


#Plot the last block that doesn't plot in the for loop.             
plt.bar(bstart[-1]-bstart[0],rates[-1],width = x[-1]-bstart[-1],
            alpha=0.5,color='r')
            

'''
#SECOND !!!!!!!!!!!



#for the second 7 
ObsID_2 = ["15043","15044","15045","15651","15654"
       ,"16210","16211","16212","16213"]

files_2 = []
i=0
for i in range(0,len(ObsID_2)):
    files_2.append(ObsID_2[i] + "_sgra_bkgsub_lc_2-8keV_1.25asec.fits")

ffile_2 = []
data_2 = []

time_2=[]
y_2 = []
err_2 = []
tstart_2 = []
tstop_2 = []
start_2 = []
end_2 = []
vertical_2 = []


WHat I did here is that I deleted the gaps between the observations
and I stacked each observation next to each other.


i=0
for i in range(0,len(files_2)):
    ffile_2.append(fits.open(files_2[i]))
    data_2.append(ffile_2[i][1].data)
    tstart_2.append(ffile_2[i][1].header['TSTART'])
    tstop_2.append(ffile_2[i][1].header['TSTOP'])
    #delete the gaps
    if i > 0: 
        time_2.extend(data_2[i].field('TIME')-tstart_2[i]+time_2[-1])
        vertical_2.append(time_2[-1])

    else: 
        time_2.extend(data_2[i].field('TIME')-tstart_2[i])
    
    #start.append(data[i][0].field('TIME')-tstart[0])
    y_2.extend(data_2[i].field('COUNT_RATE'))
    err_2.extend(data_2[i].field('COUNT_RATE_ERR'))
#SECOND !!!!!!!!!!!



#for the second 7 
ObsID_3 = ["16214","16215","16216","16217"
       ,"16218","16508","16597"]

files_3 = []
i=0
for i in range(0,len(ObsID_3)):
    files_3.append(ObsID_3[i] + "_sgra_bkgsub_lc_2-8keV_1.25asec.fits")

ffile_3 = []
data_3= []

time_3=[]
y_3 = []
err_3 = []
tstart_3 = []
tstop_3 = []
start_3 = []
end_3 = []
vertical_3 = []

flares=[]
flares.append(1500+3062.2011278/2)
flares.append(14086.69902754 + 2413.30097246/2)
flares.append(7382.38815305 + 1041.72473171/2)
flares.append(34500 + 833.64574081/2)




WHat I did here is that I deleted the gaps between the observations
and I stacked each observation next to each other.



i=0
for i in range(0,len(files_3)):
    ffile_3.append(fits.open(files_3[i]))
    data_3.append(ffile_3[i][1].data)
    tstart_3.append(ffile_3[i][1].header['TSTART'])
    tstop_3.append(ffile_3[i][1].header['TSTOP'])
    #delete the gaps
    if i > 0: 
        time_3.extend(data_3[i].field('TIME')-tstart_3[i]+time_3[-1])
        vertical_3.append(time_3[-1])

    else: 
        time_3.extend(data_3[i].field('TIME')-tstart_3[i])
    
    #start.append(data[i][0].field('TIME')-tstart[0])
    y_3.extend(data_3[i].field('COUNT_RATE'))
    err_3.extend(data_3[i].field('COUNT_RATE_ERR'))

fig, ax = plt.subplots(3,1,sharex=True, sharey=True)
#fig.tight_layout()
ax[0].errorbar(time, y,err,ecolor='g')
ax[0].vlines(vertical,0,0.05,colors='m',linestyles='dashed')
ax[0].vlines(flares,0,0.6,colors='r')

ax[1].errorbar(time_2, y_2,err_2,ecolor='g')
ax[1].vlines(vertical_2,0,0.05,colors='m',linestyles='dashed')

ax[2].errorbar(time_3, y_3,err_3,ecolor='g')
ax[2].vlines(vertical_3,0,0.05,colors='m',linestyles='dashed')
'''
plt.errorbar(time, y,err,ecolor='g')

plt.show()