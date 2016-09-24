import matplotlib.pyplot as plt 
import spinmob as s
import numpy as np

#On isis !!!!!!!!

'''
.load bblocks_examp.sl_2

start = ev.lo_t-tstart;
end = ev.hi_t-tstart;
writecol("figure.txt",start,end,ev.rate);

bstart = results.lo_t-tstart;
bend = results.hi_t-tstart;
writecol("bayesian.txt",bstart,bend,results.rate,results.err);

writecol("title.txt",plot_step,ncp_prior);
'''

'''
    .load bblocks_examp.sl_2
    
    start = ev.lo_t;
    end = ev.hi_t;
    writecol("figure.txt",start,end,ev.rate);
    
    bstart = results.lo_t;
    bend = results.hi_t;
    writecol("bayesian.txt",bstart,bend,results.rate,results.err);
    
    writecol("title.txt",plot_step,ncp_prior);
    '''

#Load the data by column for the different text files
d = s.data.load("figure.txt")

data = s.data.load("bayesian.txt")

title = s.data.load("title.txt")
step = title[0] #1st column -> plot_step for plotting
ncp = title[1] #2nd column -> significance 

plot_step = step[0]

percent = 0

#This should be improved but it setting the associating the ncp values to 
#... the percent values. 
if ncp == 1.1:
    percent = 68
if ncp == 3.0 : 
    percent = 95
if ncp == 4.6: 
    percent = 99
if ncp == 7.0: 
    percent = 99.9

start = d[0]  #start times of the lightcurve starting at 0
end = d[1]
rate = d[2]   #rates of the lightcurve for each constant bin of plot_step
#rate_err = d[3]

bstart = data[0] #start times of the blocks 
bend = data[1] 
brate = data[2] #rates of each block during the duraction of the block
berr = data[3] #error of each block

bining = []
baybin = []
#Create a bin array where the elements are the bin of lightcurve
#...therefore array of a constant value of plot_Step 
i=0 
for i in range(0,len(start)):
    bining.append(end[i]-start[i])
#CReate a bin array where the elements are the bin of each block
i=0 
for i in range(0,len(bstart)):
    baybin.append(bend[i]-bstart[i])

#With the command you kind keep editing the figure that just diplayed
plt.ion()
   
#plt.hist(start,bins=300, width=300,weights=rate,color='b')

#Plot the lightcurve as a blue histogram with constant bin = plot_step
plt.bar(start,rate,width=plot_step,color='b')

#plt.bar(bstart,brate, width=baybin,color='b',alpha=0.4,edgecolor='red',linewidth=1.5)
w=3

#plot the bayesian blocks in a steps-post mode (looks like histogram)
plt.plot(bstart,brate,linestyle='steps-post',linewidth=w,color='r')
#the last block doesn't show. so we plot it ourselves here
plt.plot((bstart[-1], bend[-1]), (brate[-1], brate[-1]), 'k-',color='r',linewidth=w)


i=0
bmid = []
for i in range(0,len(bstart)-1):
    bmid.append(bstart[i]+(bstart[i+1] - bstart[i])/2)
bmid.append(bstart[-1]+(bend[-1]-bstart[-1])/2)

#Plot the error bars for the bblocks
plt.errorbar(bmid,brate,berr,fmt='o',linewidth=w,color='m')

plt.ylabel("Count rate")
plt.xlabel("Time (s)")
plt.title("significance=%lf" %percent + "%" + " and bin=%lf" %step)

plt.show()