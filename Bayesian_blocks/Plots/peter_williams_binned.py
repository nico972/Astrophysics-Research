import astropy.io.fits
from astropy.io import fits
import matplotlib.pyplot as plt 
import numpy as np

import spinmob as s

'''
d = s.data.load("peter_cells.txt")
cell_size = d[0]
cell_pops = d[1]
'''

ffile = fits.open('14942_sgra_events_2-8keV_1.25asec.fits')

ffile2 = fits.open('14942_sgra_bkgsub_lc_2-8keV_1.25asec.fits')
data2 = ffile2[1].data

 
plot_step = 300 

head = ffile[0].header  
head2 = ffile2[1].header  

plot_counts = 'COUNTS'
plot_rate = 'COUNT_RATE'
ploterr = 'COUNT_RATE_ERR'
plotx = 'TIME'

tstart = 5.1504030562620 * 10**8
tstop  =  5.1509243712907 * 10**8

x = data2.field(plotx) - tstart
y = data2.field(plot_rate)
err = data2.field(ploterr)

time_bin = data2.field(plotx) - tstart
counts = data2.field(plot_counts)
width_array = np.empty(len(counts))
i=0
for i in range(0,len(counts)):
    width_array[i] = 300



import pwkit
from pwkit import bblocks

data = ffile[1].data
plot_t = 'TIME'


# extract the array data for the header 'TIME'
times = data.field(plot_t)


'''
array_step = []
for i in range(0,len(time_starts)):
    array_step.append(plot_step)
 
time_starts = np.arange(tstart,tstop,300)  
time_stops = time_starts + array_step    

temp = []
for i in range(0,len(time_starts)):
    temp.append(tstart)
    
starts = time_starts - temp
stops = time_stops - temp
'''

end = tstop - tstart
time_starts = np.arange(0,end,plot_step)


array_step = []
for i in range(0,len(time_starts)):
    array_step.append(plot_step)
    
time_stops = time_starts + array_step 

times_plot = times - tstart

#This is a command of the pwkit website that works

#pwkit.bblocks.tt_bblock(tstarts, tstops, times, p0=0.05))
#output = bblocks.tt_bblock(time_starts, time_stops, times_plot, p0=0.99)


#BINNED MODE !!!!!!!!!!!!!
output_bin= bblocks.bin_bblock(width_array, counts, p0=0.05)

#Using Peter cells
'''
output_bin = bblocks.bin_bblock(cell_size,cell_pops,p0=0.95)
'''
blocks = output_bin.blockstarts
rates = output_bin.rates
bining = output_bin.widths

plt.xlabel('time (s)')
plt.ylabel('COUNT_RATE (counts/s)')

plt.errorbar(x,y,err,fmt='-o',ecolor='g')

x_component = []

j = 0
for j in range(0,len(blocks)):
    x_component.append(time_bin[blocks[j]])

j=0
for j in range(0,len(rates)-1):
    plt.bar(x_component[j],rates[j], width=(x_component[j+1]-x_component[j]),
            alpha=0.3,color='r')
            
#plt.bar(x_component[-1],rates[-1],width = times_plot[-1]-x_component[-1],
    #    alpha=0.3,color='r')