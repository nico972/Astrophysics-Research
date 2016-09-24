This is the python version used to plot Bayesian Blocks (bblocks) on top of a light-curve. 

This file should be adaptable to any light-curve. Some basics features need to be set depending on which observation the data comes from. 

1) Change the names of both fits file imported. 'ffile' imports the fits file containing the time-tagged event array and the 'ffile1' imports the fits file containing the binned time array and counts per bin array. 

2) Change the start time 'tstart' and the stop time 'tstop' according to what the first file of the specific observation says. 



-----------------------
The commands of the kit:
-----------------------
The kit has three commands that can be found at http://pwkit.readthedocs.io/en/latest/science/pwkit-bblocks/
We are using the 2 last which use the time-tagged event array

---------------
Inconsistency ?
---------------
One possible inconsistency I don't yet understand is the fact that using the output (from the pwkit command) 'blockstarts' doesn't give the same blocks widths results than using the output 'ledges' which gives the time of the left edge of the block which in my opinion should be the same. For some reason 'ledge' seems to make more sense than 'blockstart'.

---------------------------
Important things to mention: 
---------------------------
1) The blocks were plotted one by one by using 'plt.bar' in a for loop because I didn't find a way to plot the blocks using an histogram command in python.
2) The widths of the blocks were plotted by using the start of the blocks and assuming that the width is the length between the next block and the present block. We didn't use the width array from the output 'width' (from the pwkit command)
3) For the 'pw_bs.py' python file where we use the third command of the kit Peter Williams mentioned on his website (mentioned above) that the uncertainties are not very good. I suppose we can consider them as a help but we can't rigorously use them. 