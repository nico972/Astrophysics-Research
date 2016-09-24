tstart = 4.8672271671264* 10**8
tstop  =  4.8674335508876 * 10**8

import astropy.io.fits
from astropy.io import fits
import matplotlib.pyplot as plt 
import numpy as np

#Allow to load column data in text or any other type of file
import spinmob as s

#That's loading the sizes of the cells 
#... and the number of counts in each cell
# it is not necessary.
#It's only if peter williams input is actually the width of the cells
#...the same way isis views it


#Open the fits frile with the time-tagged event array
#You have to have these files in the same directory as you are. 
ffile = fits.open('14703_sgra_events_2-8keV_1.25asec.fits')


#open the fits file with the binned data
ffile2 = fits.open('14703_sgra_bkgsub_lc_2-8keV_1.25asec.fits')
data2 = ffile2[1].data
 
#This is the original bining of the data. 
plot_step = 300 

#this is to see the headers
head = ffile[0].header  
head2 = ffile2[1].header


#This is the names in the different headers of the certain arrays we want to 
#access. 
plot_counts = 'COUNTS'
plot_rate = 'COUNT_RATE'
ploterr = 'COUNT_RATE_ERR'
plot_tstart = 'TSTART'
plotx = 'TIME'


#THIS IS FOR PLOTTING THE LIGHTCURVE...
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

#This is the ouput command that would output the infos about the bblocks
#...bblocks = bayesian blocks

output= bblocks.bs_tt_bblock(times,time_starts,time_stops,p0=0.95,nbootstrap=512)

#Most imoprtant outputs. 
blocks = output.blockstarts
berr = output.bsrstds
brate = output.rates
bining = output.widths
output.final = output.finalp0


plt.xlabel('time (s)')
plt.ylabel('COUNT_RATE (counts/s)')


#Plot the lightcurve
plt.errorbar(x,y,err,color='g',fmt='-o',ecolor='g')
#plt.errorbar(x[3:-1]-900,y[3:-1],err[3:-1],fmt='-o',color='g')

'''
#Create x_component array with elements being the start time of each block
#... coming from the output
#...considering the fact that we are starting at 0 also
bstart_2 = []
j = 0
while j<len(blocks):
    if blocks[j]==0:
        j= j+1
    bstart_2.append(times[blocks[j]])
    j=j+1;
'''

bstart = output.ledges
bend = output.redges

#Plot every single bar with a width being the difference between
#...the present block and the start of the next block
j=0
for j in range(0,len(brate)-1):
    plt.bar(bstart[j]-bstart[0],brate[j], width=(bstart[j+1]-bstart[j]),
            alpha=0.5,color='r')        

#for the errors
i=0
bmid = []
for i in range(0,len(bstart)-1):
    bmid.append(bstart[i]+(bstart[i+1] - bstart[i])/2)
bmid.append(bstart[-1]+(bend[-1]-bstart[-1])/2)

w=3

#Plot the error bars for the bblocks
plt.errorbar(bmid,brate,berr,fmt='o',linewidth=w,color='m')

#Plot the last block that doesn't plot in the for loop.             
plt.bar(bstart[-1]-bstart[0],brate[-1],width = x[-1]-bstart[-1],
        alpha=0.5,color='r')

