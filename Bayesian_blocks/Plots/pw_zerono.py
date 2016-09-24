OBS_ID = raw_input("What's the Obs_ID? ")
file_events = OBS_ID + "_sgra_events_2-8keV_1.25asec.fits"
file_bin = OBS_ID + "_sgra_bkgsub_lc_2-8keV_1.25asec.fits"
p_0 = raw_input("Value of p0: ")
version = raw_input("Version normal or bootstrap(or b) ?: ")

p0 = float(p_0)

import astropy.io.fits
from astropy.io import fits
import matplotlib.pyplot as plt 
import numpy as np

#Allow to load column data in text or any other type of file
import spinmob as s



#Open the fits file with the time-tagged event array
#You have to have these files in the same directory as you are. 
ffile = fits.open(file_events)

#open the fits file with the binned data
ffile2 = fits.open(file_bin)
data2 = ffile2[1].data
 
#This is the original bining of the data. 
plot_step = 300 

#this is to see the headers
head = ffile[0].header  
head2 = ffile2[1].header  

tstart = head['TSTART']
tstop = head['TSTOP']
deadtime_cor = head2['DTCOR']
exposure = head2['EXPOSURE']

date_start = head['DATE-OBS']
date_end = head['DATE-END']

ObsID = head['OBS_ID']
filename = ObsID + '.txt'

#This is the names in the different headers of the certain arrays we want to 
#access. 
plot_counts = 'COUNTS'
plot_rate = 'COUNT_RATE'
ploterr = 'COUNT_RATE_ERR'
plotx = 'TIME'


#Load the data of the particular arrays you want from the fits file
time_bin = data2.field(plotx)
#...and substart this array of time to the start time to start at 0 instead of 
#...a random time
x = time_bin - tstart
y = data2.field(plot_rate)
err = data2.field(ploterr)


#This is importing the package to output Peter Williams algorithm version
import pwkit
from pwkit import bblocks

data = ffile[1].data
plot_t = 'TIME'

# extract the array data for the header 'TIME' 
times_events = data.field(plot_t)


#For plotting purposes
end = tstop - tstart
time_starts = np.arange(0,end,plot_step)

array_step = []
for i in range(0,len(time_starts)):
    array_step.append(plot_step)
    
time_stops = time_starts + array_step 

#Again substract the start time to start at 0 instead of some random time 
times = times_events - tstart

#THIS IS THE COMMAND FROM THE KIT THAT OUTPUT THE BLOCKS DATA
if version == 'normal':
    output = bblocks.tt_bblock(time_starts, time_stops, times, p0=p0)
elif version == 'bootstrap'or version == 'b':
    output= bblocks.bs_tt_bblock(times,time_starts,time_stops,p0=p0,nbootstrap=512)
else:
    print "Choose between normal and bootstrap"

#Most imoprtant outputs. 
blocks = output.blockstarts
rates = output.rates
bining = output.widths
output.final = output.finalp0


plt.xlabel('time (s)')
plt.ylabel('COUNT_RATE (counts/s)')
plt.title('ObsID '+OBS_ID)

'''
#This is if we want to use the start of the blocks from 'blockstarts' output
bstart_2 = []
j = 0
while j<len(blocks):
    bstart_2.append(times[blocks[j]])
    j=j+1;
'''
#This is using the start of the blocks from 'ledges' output
bstart = output.ledges
bend = output.redges

#Plot the lightcurve
plt.errorbar(x,y,err,color='g',fmt='-o',ecolor='g')

if version == 'normal':
    boot = ''

#Plot every single bar with a width being the difference between
#...the present block and the start of the next block
    j=0

    for j in range(0,len(rates)-1):
        plt.bar(bstart[j]-bstart[0],rates[j], width=(bstart[j+1]-bstart[j]),
                alpha=0.5,color='r')


#Plot the last block that doesn't plot in the for loop.             
    plt.bar(bstart[-1]-bstart[0],rates[-1],width = x[-1]-bstart[-1],
            alpha=0.5,color='r')

elif version == 'bootstrap' or version == 'b':
    berr = output.bsrstds
    boot = 'bs'
 
#...the present block and the start of the next block
    j=0
    for j in range(0,len(rates)-1):
        plt.bar(bstart[j]-bstart[0],rates[j], width=(bstart[j+1]-bstart[j]),
                alpha=0.5,color='r')
    i=0
    bmid = []
    for i in range(0,len(bstart)-1):
        bmid.append(bstart[i]+(bstart[i+1] - bstart[i])/2)
    bmid.append(bstart[-1]+(bend[-1]-bstart[-1])/2)

    w=3

#Plot the error bars for the bblocks
    plt.errorbar(bmid,rates,berr,fmt='o',linewidth=w,color='m')

#Plot the last block that doesn't plot in the for loop.
    plt.bar(bstart[-1]-bstart[0],rates[-1],width = x[-1]-bstart[-1],
            alpha=0.5,color='r')
    
    bg_rate = np.amin(rates) 
    place = np.argmin(rates)
    valid = bg_rate + 2 * berr[place]
    if bg_rate ==0: 
        bg_rate = np.amin(np.delete(rates,np.argmin(rates)))
        place = np.where(rates==bg_rate)
        valid = bg_rate + 2 * berr[place[0][0]]
    
  


    rate_down_err = rates - berr
    is_fare = []
    i=0

    flare = np.where(rate_down_err>valid)
    flare = flare[0]
    
else:
    print "You need to choose between normal and bootstrap"

duration = bend - bstart


print "The flares are:" 
i=0 
for i in range(0,len(flare)):
    print "block" + str(flare[i]+1)
    


f = open(filename,'w')
#f.write(print_array(bstart,bend,rates))
f.write("ObsID = %s"%OBS_ID)
f.write("\nstart date = %s" %date_start)
f.write("\nend date = %s" %date_end)
f.write("\ntstart = %d" %tstart)
f.write("\ntstop = %d" %tstop)
f.write("\ndt_cor= %d"%deadtime_cor)
f.write("\nexposure = %d"%exposure)
f.write("\nbg_rate = %d"%bg_rate)
f.write("\nduration= ")
f.write(str(duration))

f.write("\n\nFor p0 = %f "%p0)
f.write("and version = %s \n"%version)


f.write('bstart= ')
f.write(str(bstart))
f.write('\nbstop= ')
f.write(str(bend))
f.write("\nrates= ")
f.write(str(rates))
 
if version == 'bootstrap' or version == 'b':
    f.write('\nberr= ')
    f.write(str(berr))
    i = 0
    f.write("\n\nThe flares are:")
    for i in range(0,len(flare)):
        f.write("\nblock" + str(flare[i]+1))
    
    

f.close()

plt.ion()


name_pdf = "%s_"%OBS_ID + "%s_"%boot + "%s"%p_0 +".pdf"

plt.savefig(name_pdf)

plt.close()



