
%  SITAR - S-lang/ISIS Timing Analysis Routines 
% ***  Version 1.2.0 *** Sept. 18, 2013  ***

% (Mostly) Alphabetical List of Routines in this File

% epfold_rate_use 
% fstat           
% glob_opt_use    
% lngamma       
% log_post_XXX
% mdc_use         
% msg
% null_it         
% pfold_rate_use  
% readasm_use     
% reb_rat_use     
% sitar_avg_cpd               
% sitar_avg_psd   
% sitar_bin_events            
% sitar_define_psd            
% sitar_epfold_rate           
% sitar_global_optimum        
% sitar_lags                  
% sitar_lbin_cpd              
% sitar_lbin_psd              
% sitar_make_data_cells       
% sitar_pfold_event            
%  _sitar_event_profile       
%  _sitar_pfold_event          
% sitar_pfold_rate            
%  _sitar_pfold_profile       
%  _sitar_pfold_rate          
% sitar_readasm               
% sitar_rebin_rate            

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% If you've got the s-lang GNU Scientific Module loaded,
% uncomment the following *before* running.  Necessary for statistics
% of epoch fold (dummy function used if left uncommented), and for
% additional functionality of the Bayesian Blocks routine (i.e.,
% `event mode' with alpha & beta > 0, **which won't crash, but won't
% give real results, either, if used without GSL**).

%require("gsl");

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

provide("sitar");

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

private define null_it(a,i)
{
   variable j = Char_Type[ length(a) ];
   j[i] = 1;
   return a[ where(j==0) ];
}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

variable fp = stderr;

private define mdc_use()
{
   () = fprintf(fp,  "%s\n", %{{{
`
  cell = sitar_make_data_cells(tt,type,max_delt,frame,tstart,tstop);
    
  Inputs: 
    tt          : The times of the detected events
    type        : Describe cells as 'events' (type=1 or 2), with the
                  size (i.e., normalized duration) being the distance
                  from halfway from the previous event to halfway to
                  the next event (type=1), or from the current event
                  to right before the subsequent event (type=2).
                  Alternatively, events can be assigned to bins of
                  uniform size (type=3).
    max_delt    : Events closer together than max_delt are grouped
                  together in a single cell.
    frame       : The output cell sizes are in units of frame, or this
                  is the bin size for type=3.
    tstart/tstop: start/stop times for making the output cells.
    
  Outputs: 
    cell.pops   : Array with the number of events per cell  
    cell.size   : Event mode: size (i.e. duration) of a cell in units of
                  frame.  Deadtime/efficiency should be folded into this
                  definition.  E.g., size= total "good time" / total "good
                  time" for an individual frame.  (Here, we use the equivalent
                  for uniform deadtime, size = total time / total frame time)
                  Binned mode: size is mean (over whole observation) counts
                  per bin.  (Again, here is where you would put in efficiency 
                  factors, *if* efficiency/deadtime is not uniform bin to bin)
    cell.lo/hi_t: cell start and stop times
    cell.dtcor  : Array, currently set to unity, for storing 'dead 
                  time' or 'efficiency' corrections for the cells.
                  Must be set, after the fact, by the user.
                  Currently only used in sitar_global_optimum to
                  give output block rates corrected for dead time    
`);  %}}}    
   return;
}
  

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

public define sitar_make_data_cells( )
{ 
%     Based upon Jeff Scargle's MATLAB routines to perform a 
%     Bayesian blocks decomposition of an ordered array.
%
%     Paper: "Studies in Astronomical Time Series Analysis. VI.
%             Optimum Partition of the Interval: Bayesian Blocks,
%	      Histograms, and Triggers", J.D. Scargle et al., 2003
%             (to be submitted).
%		
%     Web Site:  http://astrophysics.arc.nasa.gov/~jeffrey

   switch ( _NARGS )
   {
      case 0:
      mdc_use( );
      return;
   }
   {
      variable tt,type,max_delt,frame,tstart,tstop;
      (tt,type,max_delt,frame,tstart,tstop) = ( );
   }

   %  Initially, one datum per cell.

   variable cell = struct{ pops, size, lo_t, hi_t, dtcor };
   tt = tt[ where( tt >= tstart and tt <= tstop ) ];
   cell.pops = ones( length(tt) );

   variable mintt = min(tt);
   variable maxtt = max(tt);
   variable difftt = shift(tt,1) - tt;
   difftt = difftt[ [0:length(difftt)-2] ];

   variable ii_close = where( difftt < max_delt );
   variable ii_start, ii_beyond, ii_end, ii_clump;
   variable  clump_pop, clump_tt, ind_null;

   if ( type !=3 ) while( length(ii_close) )
   {
      ii_start = ii_close[0];
      ii_beyond = where( shift(ii_close,1) - ii_close > 1 );

      if ( length(ii_beyond) == 0 )
      {
         ii_end = ii_start + length(ii_close);
      }
      else
      {
         ii_end = ii_start + ii_beyond[0] + 1;
      }

      ii_clump = [ii_start:ii_end];

      clump_pop = sum( cell.pops[ii_clump] );
      clump_pop = typecast(clump_pop,Integer_Type);
      clump_tt = mean( tt[ii_clump] );

      % Put the members of a clump into one cell

      cell.pops[ii_start] = clump_pop;
      tt[ii_start] = clump_tt;

      % Null the cells evacuated by this operation

      if ( ii_end > ii_start+1 )
      {
         ind_null = [ii_start+1:ii_end];
      }
      else
      {
         ind_null = [ii_end:ii_end];
      }
      cell.pops = null_it(cell.pops,ind_null);
      tt = null_it(tt,ind_null);

      difftt = shift(tt,1) - tt;
      difftt = difftt[ [0:length(difftt)-2] ];
      ii_close = where( difftt < max_delt );
   }

   % Now define the cell sizes and locations.

   variable dt = shift(tt,1) - tt;
   variable ndt = length(dt)-1;
   dt = dt[[0:ndt-1]];
   variable cstart = 1, cstop = ndt-1;

   switch (type)
   { 
    case 1:

    % type = 1: Set to midpoints

      cell.size = 0.5 * (dt[[0:ndt-2]] + dt[[1:ndt-1]]);
      cell.lo_t = tt[[1:ndt-1]]-dt[[0:ndt-2]]/2.;
      cell.hi_t = tt[[1:ndt-1]]+dt[[1:ndt-1]]/2.;

      % If tstart is < first event, add in the first event.  Design
      % feature for sparse lightcurves, where a long wait time until
      % the first event might mean something.

      if ( tstart < mintt )
      {
         cell.size = [ tt[0]-tstart+dt[0]/2., cell.size ];
         cell.lo_t = [ tstart, cell.lo_t ];
         cell.hi_t = [ tt[0]+dt[0]/2., cell.hi_t ];
         cstart = 0;
      }

      % If tstop is > last event, add in the last event.  Design
      % feature for sparse lightcurves, where a long drought at the
      % end might mean something.

      if ( tstop > maxtt )
      {
         cell.size = [ cell.size, tstop-tt[ndt]+dt[ndt-1] ];
         cell.lo_t = [ cell.lo_t, tt[ndt-1]+dt[ndt-1]/2. ];
         cell.hi_t = [ cell.hi_t, tstop ];
         cstop = ndt;
      }

      % Now make cell.pops match the above.

      cell.pops = cell.pops[ [cstart:cstop] ];
   }

   {
    case 2:

    % type =2: Intervals instead of midpoints. 

      cell.size = dt;
      cell.lo_t = tt[ [0:ndt-1] ];
      cell.hi_t = tt[ [1:ndt] ];
      cstart = 0;
      cstop = ndt-1;

      % If tstop is > last event, add in the last event.  Design
      % feature for sparse lightcurves, where a long drought at the
      % end might mean something.

      if ( tstop > maxtt )
      {
         cell.size = [ cell.size, tstop-tt[ndt] ];
         cell.lo_t = [ cell.lo_t, tt[ndt] ];
         cell.hi_t = [ cell.hi_t, tstop ];
         cstop = ndt;
      }

      % I'm not entirely happy with any choices for the endpoints for
      % sparse lightcurves.  But one has to choose *something*.  For
      % tstart < mintt, append that time to the first bin.  Thought
      % for alternative and type=2: add a bin with 0 photons going
      % from tstart to right before the first photon?  (Should be OK
      % with the statistic?)

      if ( tstart < mintt )
      {
         cell.size[0] = cell.size[0] + ( tt[0]-tstart );
         cell.lo_t[0] = tstart;
      }
      cell.pops = cell.pops[ [cstart:cstop] ];
   }

   {
    case 3:

    % type = 3: Bin using ISIS functions

      variable nbins = typecast( (tstop-tstart)/frame, Integer_Type );
      ( cell.lo_t, cell.hi_t ) = linear_grid( tstart, tstop, nbins);

      cell.pops = histogram( tt, cell.lo_t, cell.hi_t );
      cell.dtcor = ones( length(cell.pops) );
      cell.size = @cell.dtcor;

      % Based upon my thoughts of applying this to gratings data, I've
      % decided upon the following modification.  This forces the
      % "rate per cell size" to be typically less than one.  Or, in
      % terms of the marginalized likelihood, this puts a sum of
      % "model counts per block" in the denominator (presuming that
      % your zeroth order model is `constant rate over the whole
      % observation').

      cell.size = cell.size * ( sum(cell.pops)/length(cell.size) );
      return cell;
   }

   {

    % If it's not type 1, 2, or 3, something is wrong! 

      error("Type not 1, 2, or 3 in 'sitar_make_data_cells'");
   }

   % The return for Type 1 or 2

   cell.size = cell.size/frame;
   cell.dtcor = ones( length(cell.pops) );	
   return cell;
}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#ifnexists GSL_EUNDRFLW
% If GSL is not defined, define our own
private define lngamma(z)
{
   variable c = [76.18009173,-86.50532033,24.01409822,-1.231739516,
                 1.20858003e-3,-5.36382e-6];
   variable csum=1.;
   foreach(c)
   {  
      csum+=()/z;
      z+=1.;
   }
   csum = log( sqrt(2.*PI)*csum ) - (z-1.5) + (z-6.5)*log(z-1.5);
   return csum;
}
#endif

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#ifexists GSL_EUNDRFLW
   gsl_set_error_disposition (GSL_EUNDRFLW, 0);
#endif

#ifnexists beta_inc
   private define beta_inc(x,y,z)
   {
      return 10.;
   }
#endif

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% log_post_XXX functions
%
% PURPOSE:
%     Used in Bayesian decomposition of *time tagged event* data,
%     or *binned* data into blocks with statistically uniform "rate".
%     Note: this does not have to be used on a lightcurve, but instead
%     could be used, e.g., on a spectrum to identify potential lines.  
%     The only requirement is that the data be sequential, with an
%     independent variable (e.g., time, wavelength, etc.), and a dependent
%     variable (e.g., counts, rate).
%
%     This function calculates the log of the probability function (or 
%     fitness function) for a given decomposition.  Called as a subroutine
%     from "global_optimum".
%
% CALLING SEQUENCE:
%     log_prob = log_post_XXX(data_pops, data_size, ncp_prior, alpha, beta);
%
% INPUTS:
%     data_pops  : Array with cell populations (i.e., number of "events" 
%                  per cell)
%     data_size  : Array with cell sizes (i.e., number of "ticks" 
%                  per cell; can contain window function/efficiency)
%     ncp_prior  : log parameter for number of change points.  It is
%                  assumed that the prior probability of # of change
%                  points goes as gamma^N_cp, then:
%	              ncp_prior = -log(gamma)
%                  MAN suggests: ncp_prior = log(# of bins in lightcurve/
%                                             maximum # of desired cells)
%                  ncp_prior=7 is ~ 10^-3 significance for each new block
%     alpha, beta: Parameters that affect choices on priors.  
%                  See glob_opt_use.
%
% OUTPUTS:
%     log_prob   : The logarithm of the probability function for a given
%                  decomposition of the lightcurve.
% 
% MODIFICATION HISTORY:
%     Based upon Jeff Scargle's MATLAB routines to perform a 
%     Bayesian blocks decomposition of an ordered array.
%
%     Paper: "Studies in Astronomical Time Series Analysis. VI.
%             Optimum Partition of the Interval: Bayesian Blocks,
%	      Histograms, and Triggers", J.D. Scargle et al., 2003
%             (to be submitted).
%		
%     Web Site:  http://astrophysics.arc.nasa.gov/~jeffrey
%
% Note: no error checking, since this subroutine shouldn't be called 
% directly, only through "sitar_global_optimum();"

private define log_post_bin( data_pops, data_size, ncp_prior, alpha, beta)
{
   variable dpop_plus_alpha = data_pops + alpha;
   variable addenda = alpha*log( beta ) - lngamma( alpha ) - ncp_prior;
   variable log_prob = lngamma( dpop_plus_alpha ) 
                       - dpop_plus_alpha*log( data_size + beta ) + addenda;
   return log_prob;   
}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

private define log_post_max( data_pops, data_size, ncp_prior, alpha, beta )
{
   variable log_prob = data_pops*( log( (data_pops+1.e-12)/data_size ) - 1) - ncp_prior;
   return log_prob;   
}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

private define log_post_evt_rate( data_pops, data_size, ncp_prior, alpha, beta )
{
   variable log_prob;

   variable arg = data_size - data_pops;
   variable ii = where( arg > 0 );
   variable addenda = -log( beta - alpha ) - ncp_prior;

   % Slightly different tack than Scargle here.  His MATLAB code
   % sets this to ~0.  I set this M=N in Eq. (23) of Scargle (1998,
   % ApJ, 504, pp. 405-418).  Neither is strictly correct, but I
   % don't know if it makes a substantial difference.  (If it does,
   % you should be using the binned version instead.)  

   % The following represent my own thoughts for this case.

   variable ealpha = @data_size, ebeta = @ealpha;
   ealpha[*] = exp(-alpha);
   ebeta[*] = exp(-beta);
      
   % This is only an approximate expression, to make up for my
   % lack of integration skills
         
   log_prob = data_size * log(1.-ebeta) + log(beta)
              - log( beta - alpha );

   % This, however, I believe, is correct.  It *will not*
   % evaluate correctly if the GSL module is not initialized

   log_prob[ii] = lngamma( arg[ii] ) +
                  lngamma( data_pops[ii] + 1 ) -
                  lngamma( data_size[ii] + 1 ) +
                  log( 
                      beta_inc(arg[ii],data_pops[ii]+1,ealpha[ii]) -
                      beta_inc(arg[ii],data_pops[ii]+1,ebeta[ii])  
                      );
                  
   log_prob += addenda;

   return log_prob;   
}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

private define log_post_evt_mean( data_pops, data_size, ncp_prior, alpha, beta )
{
   variable log_prob;

   variable arg = data_size - data_pops;
   variable ii = where( arg > 0 );
   variable dsize_plus_one_over_alpha = data_size + 1./alpha;
   variable dsize_plus_one_over_alpha_plus_one = __tmp(dsize_plus_one_over_alpha) + 1;
   variable addenda = -log(alpha) - ncp_prior;

   log_prob = lngamma( dsize_plus_one_over_alpha ) 
            - lngamma( dsize_plus_one_over_alpha_plus_one );

   log_prob[ii] = lngamma( data_pops[ii] + 1 ) 
                + lngamma( arg[ii] + 1./alpha )
                - lngamma( dsize_plus_one_over_alpha_plus_one[ii] );

   log_prob += addenda;

   return log_prob;   
}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

private define log_post_evt_alpha( data_pops, data_size, ncp_prior, alpha, beta)
{
   variable log_prob;

   variable arg = data_size - data_pops;
   variable ii = where( arg > 0 );
   variable dsize_plus_alpha_plus_two = data_size + alpha + 2;
   variable addenda = log( alpha+1 ) - ncp_prior;

   log_prob = lngamma( data_size + 1 ) + lngamma( alpha + 1 )
            - lngamma( dsize_plus_alpha_plus_two );

   log_prob[ii] = lngamma( data_pops[ii] + 1 ) 
                + lngamma( arg[ii] + alpha + 1 )
                - lngamma( dsize_plus_alpha_plus_two[ii] );

   log_prob += addenda;

   return log_prob;   
}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

private define log_post_evt( data_pops, data_size, ncp_prior, alpha, beta )
{
   variable log_prob;

   variable arg = data_size - data_pops;
   variable ii = where( arg > 0 );
   variable dsize_plus_two = data_size + 2.;

   log_prob = lngamma( data_size + 1. ) 
            - lngamma( dsize_plus_two );

   log_prob[ii] = lngamma( data_pops[ii] + 1 ) 
                + lngamma( arg[ii] + 1.  )
                - lngamma( dsize_plus_two[ii] );
 
   log_prob -= ncp_prior;

   return log_prob;   
}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


private define glob_opt_use()
{
   () = fprintf(fp,  "%s\n", %{{{
`
   ans = sitar_global_optimum( cell, ncp_prior, type 
                               [; first, mean_prior, rate_prior, alpha=#, beta=#] );
    
   Inputs:
     cell      : A structure, with cell.pops, cell.size, cell.lo_t
                 cell.hi_t,cell.dtcor (see sitar_make_data_cells)
     ncp_prior : Parameter for prior on number of 'blocks' 
     type      : Identical to type from sitar_make_data_cells

   Qualifiers:
     first     : If set, the code runs in 'trigger' mode, and returns upon
                 the first sign of a change, so long as it is greater
                 than 'first' cells from the beginning of the
                 lightcurve

                EVENT MODE PRIORS (type=1 or 2)-

                 Default prior is that p_1, the probability of one or
                 more photons in a frame, is uniformly dsitributed
                 from 0->1.

     alpha     : >=0 implies p_1 prior is (1+alpha) (1-p_1)^alpha 

     mean_prior: if set, prior on *rate* is exp(-Lambda/Lambda_0),
                 where \Lambda_0 is the mean rate from the entire
                 observation.  Supercedes any choice on alpha.

     rate_prior: if set, prior is chosen to be uniform between
     rate_low,   rates ranging from rate_low -- rate_high, with
     rate_high   defaults of rate_low = 1/3 mean rate and 
                 rate_high = 3 times mean rate. Supercedes any
                 choice of alpha or mean_prior.

                BINNED MODE PRIORS AND POSTERIORS (type=3)-

     alpha,    : The prior on the bin rate (sort of) goes as
     beta        (Lambda)^(alpha-1) Exp(-beta*Lambda), where
                 Lambda=rate, and the default is alpha=beta=1

     max_like  : Use maximum likelihood, log(Pmax) = N log(N/T) -N,
                 as the posterior test function.  Overrides any
                 choice of alpha & beta for binned mode.

   Outputs:
     results.cpt  : Array of change point locations for the maximum 
                    likelihood solution (indices are specific to the
                    cell input);
     results.last,: (Diagnostic purposes only) Arrays of the
            .best   location of the last change point and the
                    associated maximum log probability
     results.cts  : Counts in each block
     results.rate : Rate in each block
     results.err  : Poisson error for the block rate
     results.lo_t,: Times of lower (>=) and upper (<) block
            .hi_t   boundaries.
`); %}}}
   return;
}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

public define sitar_global_optimum()
{
%     Used in Bayesian decomposition of *time tagged event* data,
%     or *binned* data into blocks with statistically uniform "rate".
%     Note: this does not have to be used on a lightcurve, but instead
%     could be used, e.g., on a spectrum to identify potential lines.  
%     The only requirement is that the data be sequential, with an
%     independent variable (e.g., time, wavelength, etc.), and the 
%     dependent variable (i.e., counts)
%
%     This is the main routine to determine the optimal partitioning of 
%     the sequential data.
% 
% MODIFICATION HISTORY:
%     Based upon Jeff Scargle's MATLAB routines to perform a 
%     Bayesian blocks decomposition of an ordered array.
%
%     Paper: "Studies in Astronomical Time Series Analysis. VI.
%             Optimum Partition of the Interval: Bayesian Blocks,
%	      Histograms, and Triggers", J.D. Scargle et al., 2003
%             (to be submitted).
%		
   % Some really basic error checking for the input

   variable cell, ev_type, ncp_prior, first, alpha, beta, log_prob, lo, ab_default;
   switch ( _NARGS )
   {
    case 3:
      ( cell, ncp_prior, ev_type ) = ();
   }
   {
      glob_opt_use();
      if(_NARGS != 0)
      {
         () = fprintf(fp, "\n Incorrect number of arguments. \n");
         _pop_n(_NARGS);
      }
      return;
   }

   if ( is_struct_type(cell) )
   {
      variable fields = get_struct_field_names(cell);
      if ( length( where( fields == "pops" or fields == "size" or
                          fields == "lo_t" or fields == "hi_t" or
                          fields == "dtcor" ) ) != 5 )
      {
         () = fprintf(fp, "\n Fields on cell structure not properly set.\n");
         return;
      }
   }
   else
   {
      () = fprintf(fp, "\n 'cell' was not a structure\n");
      return;
   }

   if (
        typeof(cell.pops) != Array_Type or 
        typeof(cell.size) != Array_Type or 
        typeof(cell.lo_t) != Array_Type or 
        typeof(cell.hi_t) != Array_Type or 
        typeof(cell.dtcor)  != Array_Type 
      ) 
   {
       
      () = fprintf(fp,"\n Error: All fields of structure 'cell' must be Array_Type.\n");
      return;
   }

   variable num_cells = length(cell.pops);

   if ( num_cells < 3 ) 
   {
      () = fprintf(fp, "\n Error: length of lightcurve is %d. Needs >= 3.\n", num_cells);
      return;
   }

   if (
        num_cells != length(cell.size) ||
        num_cells != length(cell.lo_t) ||
        num_cells != length(cell.hi_t) ||
        num_cells != length(cell.dtcor) 
      ) 
   {
      () = fprintf(fp, "\n Error: All arrays in cell structure must have same length\n");  
      return;
   }

   if ( ev_type != 1 && ev_type != 2 && ev_type != 3 ) 
   {
      () = fprintf(fp, "\n Error: 'ev_type' set to %d, but must be 1, 2, or 3", ev_type);
      return;
   }

   first=0;
   if(qualifier_exists("first")) first=1;

   variable lpost_fun;

   if(ev_type == 3 && qualifier_exists("max_like"))
   {
      alpha = 1.;
      beta = 1.;
      lpost_fun = &log_post_max;
   }
   else if(ev_type == 3)
   {
      alpha = qualifier("alpha",1.);
      beta = qualifier("beta",1.);
      lpost_fun = &log_post_bin;
   }
   else if(qualifier_exists("rate_prior"))
   {
      alpha = qualifier("rate_low",sum(cell.pops)/sum(cell.size)/3.);
      beta = qualifier("rate_high",3.*sum(cell.pops)/sum(cell.size));
      lpost_fun = &log_post_evt_rate;
   }
   else if(qualifier_exists("mean_prior"))
   {
      alpha = qualifier("rate_low",sum(cell.pops)/sum(cell.size));
      beta = 1.;
      lpost_fun = &log_post_evt_mean;
   }
   else if(qualifier_exists("alpha"))
   {
      alpha = qualifier("alpha",1);
      beta = 1.;
      lpost_fun = &log_post_evt_alpha;
   }
   else
   {
      alpha = 1.;
      beta = 1.;
      lpost_fun = &log_post_evt;
   }

   if ( alpha < 0. || beta < 0. )
   {
      () = fprintf(fp, "\n Error: qualifier values must be >= 0.\n");
      return;
   }

   variable ans = struct{cpt,last,best,cts,rate,err,lo_t,hi_t};
   variable iw, intvl;
   variable best=[0.], last=[0];

   % Long ago there was a reason I did it this way ...

   best[0] =  max( @lpost_fun( cell.pops[[0:0]],cell.size[[0:0]],
                               ncp_prior, alpha, beta ) );

   % Now onto the meat of the job. Slightly longer than the one line
   % in MATLAB.

   variable besttest;
   variable i, irange;
   for( i=1; i<= num_cells-1 ; i++ )
   {
      irange = [i:0:-1];
      besttest =  [0.,best] +
                  reverse(
                            @lpost_fun( 
                                        cumsum( cell.pops[irange] ),
                                        cumsum( cell.size[irange] ),
                                        ncp_prior, alpha, beta
                                      )
                         );

      best = [best, max(besttest)];
      last = [last, min( where( besttest == best[i] ) )];

      % Is this trigger mode?

      if ( first and last[i] ) 
      {
         ans.best = best;
         ans.last = last;

         ans.cpt = [ 0, last[i] ];
	 ans.cts = Double_Type[2];
	 ans.rate = @ans.cts;
	 ans.err = @ans.cts;
	 ans.lo_t = @ans.cts;
	 ans.hi_t = @ans.cts;

         variable ia = Array_Type[2];
         ia[0] = [ 0: last[i]-1 ];
         ia[1] = [ last[i]: i ];
         i = 0;
	 loop(2)
         {
            ans.cts[i] = sum( cell.pops[ ia[0] ] );
            ans.err[i] = ( 1. + sqrt(0.75 + ans.cts[i] ) ); 

            ans.lo_t[i] = cell.lo_t[ min( ia[0] ) ];
            ans.hi_t[i] = cell.hi_t[ max( ia[0] ) ];

            intvl = sum( ( cell.hi_t[ ia[0] ] - cell.lo_t [ ia[0] ] ) 
                         * cell.dtcor[ ia[0] ] );
            ans.rate[i] = ans.cts[i] / intvl;
            ans.err[i] = ans.err[i] / intvl;

	    i++;
         }
         ans.cts = typecast(ans.cts,Integer_Type);
         return ans;
      }
   }

   % Find the optimum partition

   variable index=last[num_cells-1];
   ans.cpt = [index];
   index = last[index-1];	

   % Long ago there was a reason why I did it this way...

   while ( index )
   {
      ans.cpt = [ index, ans.cpt ];
      index = last[index-1];
   }

   % Always return the first changepoint index as 0, if not already 0
   % (as opposed to MATLAB, the latter happens because of how we
   % initialized 'last').

   if ( ans.cpt[0] )
   {
      ans.cpt = [ 0, ans.cpt ];
   }
   ans.best=best;
   ans.last=last;

   % Make those useful things for plotting, etc.

   variable cpt_length = length(ans.cpt);
   ans.rate = Double_Type[cpt_length];
   ans.err = @ans.rate;
   ans.cts = @ans.rate;
   ans.lo_t = @ans.rate;
   ans.hi_t = @ans.rate;

   ans.lo_t[0] = cell.lo_t[0];
   ans.hi_t[cpt_length-1]=cell.hi_t[length(cell.hi_t)-1];
   i = 0;
   loop ( cpt_length-1 )
   {
      iw = [ ans.cpt[i]:ans.cpt[i+1] - 1 ];

      ans.cts[i] = sum( cell.pops[iw] );
      ans.err[i] = ( 1. + sqrt(0.75 + ans.cts[i] ) ); 

      ans.hi_t[i] = cell.hi_t[ ans.cpt[i+1] - 1 ];
      ans.lo_t[i+1] = cell.lo_t[ ans.cpt[i+1] ];
%
      intvl = sum( ( cell.hi_t[iw] - cell.lo_t [iw] ) * cell.dtcor[iw] );
      ans.rate[i] = ans.cts[i] / intvl;
      ans.err[i] = ans.err[i] / intvl;

      i++;
   }

   % Mop up that endpoint

   iw = [ ans.cpt[i]:length(cell.hi_t) - 1 ];

   ans.cts[i] = sum( cell.pops[iw] );
   ans.err[i] = ( 1. + sqrt(0.75 + ans.cts[i] ) ); 

   intvl = sum( ( cell.hi_t[iw] - cell.lo_t [iw] ) * cell.dtcor[iw] );
   ans.rate[i] = ans.cts[i] / intvl;
   ans.err[i] = ans.err[i] / intvl;

   ans.cts = typecast(ans.cts,Integer_Type);
   return ans;
}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

private define readasm_use()
{
   () = fprintf(fp,  "%s\n", %{{{
`
   event = sitar_readasm( file [; tstart=#, tstop=#, maxchi2=#, chnl=#, mjd, jd] );
    
     Variables omitted or set to 0 take on default values.
     Variables in [] are optional, but are order specific.
    
   Inputs:
     file      : ASM filename, e.g. 'xa_cygx1_d1.lc' or 'xa_cygx1_d1.col'

   Qualifiers:
     tstart    : Start time to be read (units of MET or MJD)
     tstop     : Stop time to be read (units of MET or MJD)
     maxchi2   : The maximum reduced chi2 for an 'ASM solution' which
                 will be accepted to consider a data point
     chnl      :  =0 for total band (source.lc files = default)
                 !=0 for three ASM colors + total (source.col files) 
     mjd       : If set changes from default of Mission Elapsed Time (days)
                 to MJD = Julian Date - 2,400,000.5
     jd       : If set changes from default of Mission Elapsed Time (days)
                 to Julian Date.  Supersedes mjd flag.
    
   Outputs: 
     event.time: Time of each event
     event.rate: Total ASM count-rate
     event.err : Uncertainty in count rate
     event.chi2: Chi2 of ASM solution
    
   Optional Outputs:
     event.[ch1,ch2,ch3]_rate: ASM count rate in each channel
     event.[ch1,ch2,ch3]_err : ASM count rate error in each channel
     event.[ch1,ch2,ch3]_chi2: Chi2 for ASM solution in each channel
                               (Overrides event.chi2, which won't be output)
`); %}}} 
   return;
}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

public define sitar_readasm( )
{
   variable file, event;
   variable tstart  = qualifier("tstart",-1.e32);
   variable tstop   = qualifier("tstop",1.e32);
   variable maxchi2 = qualifier("maxhi2",1.e6);
   variable chnl    = qualifier("chnl",0);
   variable mjd     = qualifier_exists("mjd");
   variable jd      = qualifier_exists("jd");

   if (jd ) mjd=1;

   switch (_NARGS)
   {
    case 1:
      ( file ) = ();
   }
   {
      if(_NARGS != 0)
      {
         () = fprintf(fp, "\n Incorrect number of arguments. \n");
         _pop_n(_NARGS);
      }
      readasm_use();
      return;
   }

   if ( orelse{tstart == NULL}{tstart == 0.})
   {
      tstart = -1.e32;
   }
   if ( orelse{tstop == NULL}{tstop == 0.}{tstop<tstart})
   {
      tstop = 1.e32;
   }
   if ( orelse{maxchi2 == NULL}{maxchi2 <= 0.})
   {
      maxchi2 = 1.e6;
   }

   file = strtrans(file," ","");

   variable t, rate, err, chi2, minchan, maxchan, band;

   (t,rate,err,chi2,minchan,maxchan,band) = 
      fits_read_col(file,"TIME","RATE","ERROR","RDCHI_SQ","MINCHAN","MAXCHAN", "BAND");

   variable mjdrefi = fits_read_key(file,"MJDREFI");
   variable mjdrefr = fits_read_key(file,"MJDREFF");
   mjdrefr = mjdrefr + mjdrefi;

   if ( mjd )     % Convert Mission Elapsed time to MJD if Flag Set
   {
      t += mjdrefr;
   }
   if ( jd )
   {
      t += 2400000.5;
   }

   variable mint = min(t), maxt = max(t);
   if( maxt < tstart or mint > tstop )
   {
      () = fprintf(fp, "\n No times specified within [tstart,tstop]\n");
      () = fprintf(fp, "   [tstart,tstop] = [%e,%e]\n",tstart,tstop);   
      () = fprintf(fp, "   [min(t),max(t)] = [%e,%e]\n",mint,maxt); 
      return NULL;
   }

   variable mincolor = [ 410, 410,1189,1861];
   variable maxcolor = [4750,1189,1861,4750];
   variable ndx;

   if ( chnl != 0 )
   {
      event = struct{time, rate, ch1_rate, ch2_rate, ch3_rate, 
                     err, ch1_err, ch2_err, ch3_err, 
                     ch1_chi2, ch2_chi2, ch3_chi2};

      variable ndxa=where(band=="A");
      variable ndxb=where(band=="B");
      variable ndxc=where(band=="C");

      variable ta = t[ndxa], tb = t[ndxb], tc = t[ndxc];
      variable iw = where( ta == tb and tb == tc );

      if( length(iw) != length(ta) )
      {
         () = fprintf(fp, "\n Color channels not of equal length\n");
	return NULL;
      }

      event.time   = t[ndxa];

      event.ch1_rate =  rate[ndxa];
      event.ch1_err  = err[ndxa];
      event.ch1_chi2 =  chi2[ndxa];

      event.ch2_rate =  rate[ndxb];
      event.ch2_err  = err[ndxb];
      event.ch2_chi2 =  chi2[ndxb];

      event.ch3_rate =  rate[ndxc];
      event.ch3_err  = err[ndxc];
      event.ch3_chi2 =  chi2[ndxc];

      ndx = where( event.ch1_chi2 <= maxchi2 and
                   event.ch2_chi2 <= maxchi2 and
                   event.ch3_chi2 <= maxchi2 and
		   event.time >= tstart and event.time <=tstop );
      if( length(ndx) == 0 )
      {
        () = fprintf(fp, "\n No solutions with chi2 < %f found\n", maxchi2);
        return NULL;
      }

      event.time   =   event.time[ndx];
      event.ch1_rate = event.ch1_rate[ndx];
      event.ch1_err  =  event.ch1_err[ndx];
      event.ch1_chi2 = event.ch1_chi2[ndx];

      event.ch2_rate = event.ch2_rate[ndx];
      event.ch2_err  =  event.ch2_err[ndx];
      event.ch2_chi2 = event.ch2_chi2[ndx];

      event.ch3_rate = event.ch3_rate[ndx];
      event.ch3_err  =  event.ch3_err[ndx];
      event.ch3_chi2 = event.ch3_chi2[ndx];

      event.rate = event.ch1_rate + event.ch2_rate + event.ch3_rate;
      event.err = sqrt( event.ch1_err^2 + 
                        event.ch2_err^2 + 
                        event.ch3_err^2 );
   }
   else
   {
      event = struct{time, rate, err, chi2};

      ndx = where(band=="SUM");

      if( length(ndx) == 0 )
      {
         () = fprintf(fp, "\n Empty lightcurve!\n");
        return NULL;
      }
      event.time =     t[ndx];
      event.rate =  rate[ndx];
      event.err  = err[ndx];
      event.chi2 =  chi2[ndx];

      ndx = where( event.chi2 <= maxchi2 and
		   event.time >= tstart and event.time <= tstop );

      if( length(ndx) == 0 )
      {
        () = fprintf(fp, "\n No solutions with chi2 < %f found\n", maxchi2);
        return NULL;
      }

      event.time = event.time[ndx];
      event.rate = event.rate[ndx];
      event.err  =  event.err[ndx];
      event.chi2 = event.chi2[ndx];
   }

   return event;
}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


private define reb_rat_use()
{
   () = fprintf(fp,  "%s\n", %{{{
`
   lc = sitar_rebin_rate(t,dt,rate [,err] [;tstart=#,tstop=#,minbin=#,
                                            user_bin=array,weight,delgap]);
    
     Variables in [] are optional, but are order specific.
    
   Inputs:
     t          : Discrete times at which rates are measured
     dt         : Width of evenly spaced time bins in rebinned lightcurve
                  (Ignored if a user_bin is input.)
     rate       : (Presumed GTI & Exposure Corrected) rates

   Optional Input:
     err        : Error on the rates

   Qualifiers:
     tstart     : Beginning of first output time bin
     tstop      : End of last output time bin
     minbin     : Minimum required events per bin, else the bin is
                  set to 0, or ignored. (Default=1.)
     user_bin.lo: Lower time bounds for user defined grid
     user_bin.hi: Upper time bounds for a user defined grid
                  (Overrides dt, tstart, and tstop inputs)
     delgap     : If set, delete empty bins from the lightcurve,
                  otherwise, set them to 0 (default)
     weight     : If set and err input, the mean is weighted by 
                  the error.  I.e., mean = (\sum rate/err^2)/(\sum 1/err^2),
                  and the output variance becomes 
                  (\sum (mean-rate/err^2)^2)/(\sum 1/err^2)^2
    
   Outputs:
     lc.rate    : Mean rate of resulting rebinning
     lc.bin_lo  : Lower time bounds of rebinned lightcurve
     lc.bin_hi  : Upper time bounds of rebinned lightcurve
     lc.num     : Number of events going into a bin
     lc.var     : Variance of the rate in a time bin. (Only calculated
                  where lc.num >= 2, otherwise set to zero.) 
    
   Optional Outputs:
     lc.err     : Rate errors combined in quadrature, i.e., it's 
                  sqrt{\sum{ err^2 }}/N, where N is the number of points
                  in that particular bin (i.e., lc.num).  Only computed if 
                  err is input [otherwise, just use sqrt(lc.var)].
`);  %}}}    
   return;
}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

public define sitar_rebin_rate( )
{
   % Bin a rate lightcurve to a new scale

   variable t, rate; 
   variable err=NULL, dt, custf=0, istck;
   variable lc;

   variable tstart = qualifier("tstart",NULL);
   variable tstop  = qualifier("tstop",NULL);
   variable minbin = qualifier("minbin",1);
   variable cust   = qualifier("user_bin",NULL);
   variable wght   = qualifier_exists("weight");
   variable gap    = qualifier_exists("delgap");

   if ( gap !=0 and minbin == 0 ) { minbin = 1; };

   switch (_NARGS)
   {
    case 3:
      ( t, dt, rate ) = ( );
   }
   {
    case 4:
      ( t, dt, rate, err ) = ( );
   }
   {
      reb_rat_use();
      _pop_n(_NARGS);
      return;
   }

   if ( cust != NULL )
   {
      if ( is_struct_type(cust) == 1 && struct_field_exists(cust,"lo")
                                     && struct_field_exists(cust,"hi") )
      { 
         custf = 1;
         dt = 0;
      }
      else 
      {
         reb_rat_use();
         () = fprintf(fp,  "%s\n", %{{{
`
  Custom user_bin is incorrect.  Input structure requires: 
     variable user_bin = struct{ lo, hi };
`); %}}}
         return;
      }
   }

   if (length(t) != length(rate) )
   {
      reb_rat_use();
      () = fprintf(fp, "\n Incommensurate array lengths for time and rate\n\n");
      return;
   }

   if ( dt == NULL ) dt = 0;

   variable mint = min(t);
   variable maxt = max(t);

   if ( tstart == NULL ) tstart = mint;
   if ( tstop == NULL )  tstop = maxt;

   if( maxt < tstart or mint > tstop or mint > maxt)
   {
      reb_rat_use( );
      () = fprintf(fp, "\n Data outside of bounds of start & stop times:\n");
      () = fprintf(fp, "      [tstart,tstop] = [%e,%e]\n",tstart,tstop);
      () = fprintf(fp, "      [min_t, max_t] = [%e,%e]\n\n",mint,maxt);
      return;
   }

   lc = struct{ rate, bin_lo, bin_hi, num, var };
   variable errf = 0;

   if( err != NULL && length(err) == length(rate) )
   {
      lc = struct_combine(lc,"err");
      errf = 1;
   }
   else if ( err != NULL )
   {
      reb_rat_use();
      () = fprintf(fp, "\n Incommensurate array length for error array\n\n");
      return;
   }

   if( custf == 0 && dt > 0)
   {
      variable nfrac = ( tstop - tstart )/dt;
      variable nbins = int(nfrac);
      lc.bin_lo = [0:nbins-1]*dt+tstart;
      lc.bin_hi = [1:nbins]*dt+tstart;
      if( nfrac-nbins > 0 )
      {
         lc.bin_lo = [lc.bin_lo,lc.bin_hi[nbins-1]];
         lc.bin_hi = [lc.bin_hi,tstop]; 
      }
   }
   else if ( custf )
   {
      lc.bin_lo = @cust.lo;
      lc.bin_hi = @cust.hi;
   }
   else
   {
      reb_rat_use();
      () = fprintf(fp, "\n New time bin or custom bin array must be input.\n\n");
      return;
   }

   variable rev;

   lc.num = histogram(t, lc.bin_lo, lc.bin_hi, &rev);
   lc.rate = Double_Type[length(lc.num)]; 
   lc.var = @lc.rate;
   if ( errf ) lc.err = @lc.rate;

   variable i=0, mom, mome;
   loop( length(lc.num) )
   {
      if ( lc.num[i] )
      {
         if( errf !=0 && wght !=0 )
         {
            mom = moment( rate[rev[i]]/err[rev[i]] );
            mome = moment( 1/err[rev[i]] );
         }
         else
         {
            mom = moment( rate[rev[i]] );
            mome = @mom;
            mome.ave = 1.;
            mome.var = 1.;
         }

         lc.rate[i] = mom.ave/mome.ave;
	 lc.var[i] = mom.var/(mome.ave)^2;

	 if ( errf )
	 {
	    lc.err[i] = sqrt(sum( (err[rev[i]])^2 ))/lc.num[i];
	 }
      }
      i++;
   }
   
   variable idx;
   if ( minbin )
   {
      if ( gap )
      {
         idx = where( lc.num >= minbin );

	 if ( length(idx) )
         {
            lc.rate = lc.rate[idx];
	    lc.bin_lo = lc.bin_lo[idx];
	    lc.bin_hi = lc.bin_hi[idx];
	    lc.var = lc.var[idx];
	    lc.num = lc.num[idx];

            if ( errf ) lc.err = lc.err[idx];
         }
	 else
	 {
            reb_rat_use( );
            () = fprintf(fp, "\n Minimum events requirement yields null lightcurve\n\n");
	    return;
	 }
      }
      else
      {
         idx = where( lc.num < minbin );

         if ( length(idx) )
         {
            lc.rate[idx] = 0.;
            lc.var[idx] = 0.;

            if ( errf )
            {
               lc.err[idx] = 0.;
            }
         }
      }
   }

   return lc;
}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

private define pfold_rate_use()
{
   () = fprintf(fp,  "%s\n", %{{{
`
   profile  = sitar_pfold_rate(t,rate,p [;pdot=#,pddot=#,nphs=#,tzero=#]);
    
     Variables omitted or set to 0 take on default values.
     Variables in [] are optional, but are order specific.
    
   Inputs:
     t             : Times at which rates are measured
     rate          : Lightcurve rate
     p             : Period on which to fold the lightcurve

   Qualifiers:
     pdot          : Period derivative
     pddot         : Period second derivative
     nphs          : Number of phase bins in output folded lightcurve.
                     Default=20.
     tzero         : Time of zero phase. Default = t[0] 
    
   Outputs:
     profile.bin_lo: Start value of phase bin
     profile.bin_hi: Stop value of phase bin
     profile.mean  : Mean rate in phase bin
     profile.var   : Variance of rate in phase bin
     profile.sdm   : Standard deviation of mean rate in a phase bin
     profile.num   : Number of data points in phase bin
`); %}}}    
   return;
}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

private define _sitar_pfold_profile(nphs)
{
   variable profile = struct{bin_lo, bin_hi, mean, var, sdm, num};
   (profile.bin_lo, profile.bin_hi) = linear_grid(0,1,nphs);
   profile.num  = Integer_Type[nphs];
   profile.mean = Double_Type[nphs];
   profile.var = Double_Type[nphs];
   profile.sdm = Double_Type[nphs];
   return profile;
}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

private define _sitar_pfold_rate( t, rate, p, pdot, pddot, nphs, profile)
{
   % Convert time to phase

   variable phs;

   switch ( 2*(pddot!=0)+(pdot!=0) )
   {
    case 0:
      phs = p;
   }
   {
    case 1:
      phs = p + pdot*t/2. ; 
   }
   {
      phs = p + ( pdot/2. + pddot*t/6. )*t ; 
   }

   phs = ( t mod phs ) / phs;

   variable ndx = where( phs < 0. );
   if( length(ndx) > 0 ) 
   {
      phs[ndx] = phs[ndx]+1.;
   }

   variable rev;

   profile.num = histogram(phs, profile.bin_lo, profile.bin_hi, &rev);

   variable i, j, mom;
   _for i (0, nphs-1, 1 )
   {
      j = rev[i];
      if ( length(j) )
      {
         mom = moment( rate[j] );
         profile.mean[i] = mom.ave;
         profile.var[i] = mom.var;
         profile.sdm[i] = mom.sdom;
      } 
   }

   return profile;
}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

public define sitar_pfold_rate( )
{
   variable t, rate, p, profile;

   variable pdot  = qualifier("pdot",0);
   variable pddot = qualifier("pddot",0);
   variable nphs  = qualifier("nphs",20);
   variable tzero = qualifier("tzero",-1.e32);

   switch (_NARGS)
   {
    case 3:
      ( t, rate, p ) = ( );
   }
   {
      pfold_rate_use();
      _pop_n(_NARGS);
      return;
   }

   variable ltime = length(t), lrate = length(rate);
   if( ltime < 2 || lrate < 2 || ltime != lrate )
   {
      pfold_rate_use();
      () = fprintf(fp, "\n Time and rate must be of equal length > 2\n");
      return;
   }

   if ( orelse{nphs == NULL }{ nphs < 1 }{nphs > ltime} )
   {
      nphs = min([20,ltime]);
   }
   else
   {
      nphs = typecast(nphs,Integer_Type);
   }
   if ( orelse{ tzero == NULL }{ tzero < -0.9e18 } )
   {
      tzero = t[0];
   }

   t = t - tzero;

   profile = _sitar_pfold_profile(nphs);
   profile = _sitar_pfold_rate( t, rate, p, pdot, pddot, nphs, profile);

   return profile;
}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

private define pfold_event_use()
{
    () = fprintf(fp,  "%s\n", %{{{
`
   profile  = sitar_pfold_event(t,p, gti_lo, gti_hi [;pdot=#,pddot=#,nphs=#,tzero=#]);
    
   Inputs:
     t             : Times at which rates are measured
     p             : Period on which to fold the lightcurve
     gti_lo        : Array of lower boundaries of good time intervals
     gti_hi        : Array of upper boundaries of good time intervals

   Qualifiers:
     pdot          : Period derivative
     pddot         : Period second derivative
     nphs          : Number of phase bins in output folded lightcurve. (Default=20)
     tzero         : Time of zero phase. (Default = t[0])
    
   Outputs:
     profile.bin_lo: Start value of phase bin
     profile.bin_hi: Stop value of phase bin
     profile.cts   : Counts in phase bin
     profile.var   : Square root of counts in phase bin
     profile.expos : Integrated exposure in phase bin
`); %}}}    
   return;
}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

private define _sitar_pfold_event_profile(nphs)
{
   variable profile = struct{bin_lo, bin_hi, cts, var, expos};
   (profile.bin_lo, profile.bin_hi) = linear_grid(0,1,nphs);
   profile.cts  = Integer_Type[nphs];
   profile.var = Double_Type[nphs];
   profile.expos = Double_Type[nphs];
   return profile;
}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

private define _sitar_pfold_event( t, p, pdot, pddot, nphs, profile, 
                                                  gti_strt, gti_stop )
{
   % Convert time to phase

   variable per, per_gti_strt, per_gti_stop, phs, 
            phs_gti_strt, phs_gti_stop, phs_lo, phs_hi;

   switch ( 2*(pddot!=0)+(pdot!=0) )
   {
    case 0:
      per = p;
      per_gti_strt = p;
      per_gti_stop = p;
   }
   {
    case 1:
      per = p + ( pdot/2. )*t ; 
      per_gti_strt = p + ( pdot/2. )*gti_strt ; 
      per_gti_stop = p + ( pdot/2. )*gti_stop ; 
   }
   {
      per = p + ( pdot/2. + pddot*t/6. )*t ; 
      per_gti_strt = p + ( pdot/2. + pddot*gti_strt/6. )*gti_strt ; 
      per_gti_stop = p + ( pdot/2. + pddot*gti_stop/6. )*gti_stop ; 
   }

   phs = ( t mod per ) / per;

   variable ndx = where( phs < 0 );
   if( length(ndx) ) 
   {
      phs[ndx] = phs[ndx]+1.;
   }

   profile.cts = histogram(phs, profile.bin_lo, profile.bin_hi);
   profile.var = sqrt(profile.cts);

   % The (absolute) start & stop phases of each Good Time Interval
   phs_gti_strt = gti_strt / per_gti_strt;
   phs_gti_stop = gti_stop / per_gti_stop;

   % Loop through each GTI

   variable i, j, phs_gti_min, phs_gti_max, tot_phs, tlo, thi, expos, 
            jfrst, jlast, jprf_strt, jprf_last, isum, phs_jfrst, phs_jlast;

   _for i (0,length(phs_gti_strt)-1,1)
   {
      % The first integer phase before the start of the Good Time
      % Intervals
      if(phs_gti_strt[i]<0)
      {
         phs_gti_min = int(phs_gti_strt[i]-1);
      }
      else
      {
         phs_gti_min = int(phs_gti_strt[i]);
      }
      % The first integer phase after the end of the Good Time
      % Intervals
      if(phs_gti_stop[i]<0)
      {
         phs_gti_max = int(phs_gti_stop[i]);
      }
      else
      {
         phs_gti_max = int(phs_gti_stop[i]+1);
      }

      % The total number of periods bracketing the Good Time Intervals
      tot_phs = phs_gti_max - phs_gti_min;
      (phs_lo,phs_hi) = linear_grid(phs_gti_min,phs_gti_max,tot_phs*nphs);

      % Convert the phases to times (including period derivatives to
      % leading order - perhaps slight overkill), and get bin exposure

      switch ( 2*(pddot!=0)+(pdot!=0) )
      {
       case 1:
         tlo = p*phs_lo;
         thi = p*phs_hi;
      }
      {
       case 2:
         tlo = p*phs_lo * (1 + pdot*phs_lo/2);
         thi = p*phs_hi * (1 + pdot*phs_hi/2);
      }
      {
         tlo = p*phs_lo * (1 + (pdot/2+pddot*p*phs_lo/3)*phs_lo);
         thi = p*phs_hi * (1 + (pdot/2+pddot*p*phs_hi/3)*phs_hi);
      }
      expos = thi-tlo;

      % The first full phase bin in the GTI
      jfrst = min(where(phs_lo>=phs_gti_strt[i]));
      % The last full phase bin in the GTI
      jlast = max(where(phs_hi<=phs_gti_stop[i]));

      % Where does the first phase bin start in the profile?
      phs_jfrst = phs_lo[jfrst];

      if(phs_jfrst<0)
      {
         phs_jfrst = (phs_jfrst mod 1) + 1;
      }
      else
      {
         phs_jfrst = (phs_jfrst mod 1);
      }

      % Where does this phase start in the output profile?  (We use
      % the profile hi, rather than the profile lo bin, to minimize
      % roundoff errors.)
      jprf_strt = min(where(profile.bin_hi>phs_jfrst));

      % Loop through all of the phase bins
      _for j (0,nphs-1,1)
      {
         isum = [jfrst+j:jlast:nphs];
         if(length(isum))
         {
            profile.expos[((jprf_strt+j) mod nphs)] = sum(expos[isum]);
         }
      }

      % The end points might have some missed exposure
      % Did we miss some exposure at the start?
      if(phs_gti_strt[i] < phs_lo[jfrst])
      {
         profile.expos[jprf_strt-1] += thi[jfrst-1]-gti_strt[i];
      }

      % Did we miss some exposure at the end?
      if(phs_gti_stop[i] > phs_hi[jlast])
      {
         phs_jlast = phs_hi[jlast];
   
         if(phs_jlast<0)
         {
            phs_jlast = (phs_jlast mod 1) + 1;
         }
         else
         {
            phs_jlast = (phs_jlast mod 1);
         }

         % Next bin center has to be greater than this (end bin, upper
         % value) phase.  I'm not sure if the try/catch is needed -
         % this is in case the last bin wraps around to 0 phase.

         try{ jprf_last = min(where((profile.bin_hi+profile.bin_lo)/2 >phs_jlast)); }
         catch AnyError: { jprf_last = 0; };

         profile.expos[jprf_last] += gti_stop[i]-thi[jlast];
      }
   }
   return profile;
}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

public define sitar_pfold_event( )
{
%
%     Fold events on a given period (plus first and second derivatives)
%
   variable profile;
   variable t, p, pdot, pddot, nphs, tzero, gti_lo, gti_hi;

   switch (_NARGS)
   {
    case 4:
      ( t, p, gti_lo, gti_hi ) = ( );
   }
   {
      pfold_event_use();
      _pop_n(_NARGS);
      return;
   }

   pdot = qualifier("pdot",0);
   pddot = qualifier("pddot",0);
   nphs = qualifier("nphs",0);
   tzero = qualifier("tzero",-1.e18);

   if ( orelse{ nphs == NULL }{ nphs < 1 } )
   {
      nphs = 20;
   }
   else
   {
      nphs = typecast(nphs,Integer_Type);
   }
   if ( orelse{ tzero == NULL }{ tzero < -0.9e18 } )
   {
      tzero = t[0];
   }

   t -= tzero;
   gti_lo -= tzero;
   gti_hi -= tzero;

   profile = _sitar_pfold_event_profile(nphs);
   profile = _sitar_pfold_event( t, p, pdot, pddot, nphs, profile, gti_lo, gti_hi);

   return profile;
}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

private define fstat(v,dfn,dfd)
{
   % Modelled after IDL f_pdf, which computes the probabilty (p) such
   % that: Probability(X <= v) = p, where X is a random variable from
   % the F distribution with (dfn) and (dfd) degrees of freedom.
   % beta_inc, is taken from the Gnu Scientific Library

   variable x, logq=10.;
   if ( v > 0 )
   {  
      logq = beta_inc( dfd/2.0, dfn/2.0, dfd/( dfd + dfn*v ) );
      logq = log10( logq );
      return logq;
   }
   else
   {
      return 1.;
   }
}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

private define epfold_rate_use()
{
   () = fprintf(fp,  "%s\n", %{{{
`
   fold  = sitar_epfold_rate(t,rate [,pstart,pstop] [;nphs=#,nsrch=#,loggrid,prds=array]);
    
     Variables omitted or set to 0 take on default values.
     Variables in [] are optional, but are order specific.
    
   Inputs:
     t         : Times at which rates are measured
     rate      : Lightcurve rate

   (Semi-)Optional Inputs:
     pstart    : Start value of periods to search
     pstop     : Stop value of periods to search

     pstart & pstop  ** or ** prds=[array] must be input

   Qualifiers:
     nphs      : Number of phase bins in folded lightcurve. (Default=20)
     nsrch     : nsrch = number of periods to search. (Default=20)
     loggrid   : If input, the period search grid is logarithmic
     prds      : Alternative array of periods to search
                 prds supersedes pstart & pstop inputs
    
   Outputs:
     fold.prds : Array of trial periods
     fold.lstat: Array of L-statistic for trial period
     fold.lsig : Array of log10 of single trial significance levels.
                 (Note: requires use of the GNU Scientific Library 
                 package, otherwise, defaults to a value of 1.)
     fold.nphs : Array of number of phases (<= nphs) at a given
                 period that had data.
`); %}}}   
   return;
}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

public define sitar_epfold_rate( )
{
   variable fold = struct{prds,lstat,lsig,nphs};
   variable t, rate, pstart=NULL, pstop=NULL;

   variable nphs  = qualifier("nphs",20);
   variable nsrch = qualifier("nsrch",20);
   variable prds  = qualifier("prds",NULL);
   variable lgrid = qualifier_exists("loggrid");

   switch (_NARGS)
   {
    case 2:
      ( t, rate ) = ( );
   }
   {
    case 4:
      ( t, rate, pstart, pstop ) = ( );
   }
   {
      epfold_rate_use();
      _pop_n(_NARGS);
      return;
   }

   variable ltime = length(t), lrate = length(rate);
   if ( ltime < 2 || lrate < 2 || ltime != lrate )
   {
      epfold_rate_use();
      () = fprintf(fp, "\n Time and rate must be of equal length > 2 \n");
      return;
   }

   % Going with input periods
   if ( prds == NULL )
   {  
      if ( orelse{ pstart == NULL }{ pstop == NULL }{ pstart<=0 }{ pstop<=0 } )
      {
         epfold_rate_use();
         () = fprintf(fp, "\n Start and stop periods must be set to positive values.\n");
         return;
      }
      else
      {
         variable pers = [pstart,pstop];
         pstart = min(pers);
         pstop = max(pers);
      }

      if ( orelse{ nphs == NULL }{ nphs < 1 }{ nphs > ltime } )
      {
         nphs = min([20,ltime]);
         () = fprintf(fp, "\n Setting nphs=20.\n");
      }
      else
      {
         nphs = typecast(nphs,Integer_Type);
      }

      if ( orelse{ nsrch == NULL }{ nsrch <= 0 } )
      {
         nsrch = 20;
         () = fprintf(fp, "\n Setting nsrch=20.\n");
      }
      else
      {
         nsrch = typecast(nsrch,Integer_Type);
      }

      variable b;
      if ( lgrid )
      {
         (,b) = linear_grid(0.,log10(pstop/pstart),nsrch-1);
         prds = pstart*[1.,10.^b];
      }
      else
      {
         (,b) = linear_grid(pstart,pstop,nsrch-1);
         prds = [pstart,b];
      }
   }
   else
   {
      if ( typeof(prds) != Array_Type )
      {
         epfold_rate_use();
         () = fprintf(fp,"\n Optional periods array 'prds' must be 'Array_Type'");
	 return;
      }
   }

   variable i, prof, lprds=length(prds);
   variable dltime=1.*ltime-nphs;
   variable factor = dltime / ( nphs -1. );

   % No loss of generality with subtracting the mean

   rate = rate - mean(rate);
   fold.lstat = @Double_Type[lprds];
   fold.lsig = @fold.lstat;
   fold.prds = prds;
   factor =  ( ltime - nphs ) / ( nphs -1. );

   prof = _sitar_pfold_profile(nphs);

   _for i (0, lprds-1, 1 )
   {
      prof = _sitar_pfold_rate(t,rate,prds[i],0,0,nphs,prof);

      fold.lstat[i] = factor * ( sum( prof.num * prof.mean^2 ) ) /
                      ( sum( prof.var * ( prof.num - 1. ) ) );

      fold.lsig[i] = fstat( fold.lstat[i], nphs-1, dltime );
   }

   fold.nphs = nphs;

   return fold;
} 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% If the GSL module is loaded, it's fft will overwrite ISIS's, and we
% have to use a different normalization.

private variable gsl_fft = 0;

#ifexists _gsl_fft_complex

gsl_fft = 1;

#endif

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

public define sitar_avg_cpd()
{
%     Take two evenly spaced lightcurves (presumed counts vs. time),
%     and calculate the PSD and CPD in segments of length l, averaged
%     over the whole lightcurve.  Segments with data gaps are skipped.
%     Deadtime corrections to Poisson noise are *not* made. Poisson
%     noise is *not* subtracted from average PSD.

   variable j,s,f,l,k,alen,ntrans;
   variable transa,transb,stransa,stransb;
   variable freq,atrans=0,btrans=0,ctrans=0.+0.i,ntot=0;
   variable mcounta=0,mcountb=0;
   variable dt, wdt, nwdt, wwdt;

   variable counts_a, counts_b, len;

   variable tbin = qualifier("tbin",1.);
            tbin = qualifier("dt",tbin);
   variable bin_time = qualifier("times",NULL);
   variable donorm = qualifier("norm",1);

   switch(_NARGS)
   {
    case 3:
      (counts_a, counts_b, len) = ();
   }
   {
     () = fprintf(fp, "%s\n", %{{{
`
   (f,psd_a,psd_b,cpd,navg,avg_a,avg_b) = sitar_avg_cpd(cnts_a,cnts_b,l [;dt=#,times=array,norm=#]);

     Take an evenly spaced lightcurve (presumed counts vs. time), and
     calculate the PSD and CPD in segments of length l, averaged over
     the whole lightcurve.  Segments with data gaps are skipped.
     Deadtime corrections to Poisson noise are *not* made. Poisson
     noise is *not* subtracted from average PSDs.

   Inputs:
     cnts_a/b: Arrays of total counts in each time bin
     l       : Length of individual PSD segments (use a power of 2!!!)

   Qualifiers:
     dt    : Length of evenly spaced bins (default == 1)
     times   : Times of measurements (otherwise presumed to have no gaps)
     norm    : Determine PSD normalization convention.
                =1, 'Leahy normalization'. Poisson noise PSD == 2, in 
                    absence of deadtime effects).
               !=1, 'rms' or 'Belloni-Hasinger' or 'Miyamoto' normalization, 
                    i.e., PSD == (rms)^2/Hz & Poisson noise== 2/Rate.

   Outputs:
     f       : PSD frequencies ( == 1/Input Time Unit)
     psd_a/b : Average PSDs 
     cpd     : Normalized Cross Power Spectral Density (complex array)
     navg    : Number of data segments going into the averages
     avg_a/b : Average number of counts per segment of length l
`);  %}}}
      return;
   }

   len = typecast(len,Long_Type);

   alen = len/2;

   % Lengths of the counts and time arrays, and the number of
   % transforms we will expect to be able to do

   variable lcounts_a, lcounts_b, ltime;

   lcounts_a = length(counts_a);
   lcounts_b = length(counts_b);

   if(bin_time != NULL)
   {
      ltime = length(bin_time);
   }
   else
   {
      bin_time = [0:lcounts_a-1]*tbin;
      ltime = lcounts_a;
   }

   if(lcounts_a != ltime or lcounts_b != lcounts_b) 
   { 
      () = fprintf(fp, "\n Arrays have unequal lengths\n");
      return; 
   }

   if(lcounts_a < len) 
   { 
      () = fprintf(fp, "\n Lightcurve too short for desired segment size\n");
      return; 
   }

   dt = bin_time - shift(bin_time,-1);
   wdt = where( dt > 1.5*tbin );

   if(length(wdt) == 0)
   {
      ntrans = lcounts_a/len;

      _for j (1,ntrans,1)
      {
         s = (j-1)*len;
         f = j*len -1;

         transa = counts_a[[s:f]];
         stransa = sum(transa);
         mcounta = mcounta + stransa;

         transb = counts_b[[s:f]];
         stransb = sum(transb);
         mcountb = mcountb + stransb;
         
         transa = shift( fft(transa,1) , 1 );
         atrans = atrans + 2*( abs(transa[[0:alen-1]]) )^2 / stransa;

         transb = shift( fft(transb,1) , 1 );
         btrans = btrans + 2*( abs(transb[[0:alen-1]]) )^2 / stransb;

         ctrans = ctrans + 2*transa*Conj(transb)/sqrt(stransa*stransb);
 
         ntot = ntot + 1;
      }
   }
   else
   {
      nwdt = length(wdt) + 1;
      wwdt = [0,wdt,lcounts_a-1];

      _for k (1,nwdt,1)
      {
         ntrans = ( wwdt[k] - wwdt[k-1] )/len;
            
         _for l (1,ntrans,1)
         {
            s = (l-1) * len + wwdt[k-1];
            f = l * len + wwdt[k-1] - 1;

            transa = counts_a[[s:f]];
            stransa = sum(transa);
            mcounta = mcounta + stransa;

            transb = counts_b[[s:f]];
            stransb = sum(transb);
            mcountb = mcountb + stransb;
         
            transa = shift( fft(transa,1) , 1 );
            atrans = atrans + 2*( abs(transa[[0:alen-1]]) )^2 / stransa;

            transb = shift( fft(transb,1) , 1 );
            btrans = btrans + 2*( abs(transb[[0:alen-1]]) )^2 / stransb;

            if ((stransa*stransb) > 0)
            {
               ctrans = ctrans + 2*transa*Conj(transb)/sqrt(stransa*stransb);
            }
            else
            {
               ctrans = ctrans + 2*transa*Conj(transb)/(sqrt(stransa*stransb*(-1))*(0.+1.0i));
            }
            ntot = ntot + 1;
         }
      }
   }

   if( ntot == 0 )
   {
      () = fprintf(fp, "\n No transforms performed!\n");
      return;
   }

   freq = [1:alen] / (len * tbin);
   mcounta = mcounta / ntot;
   mcountb = mcountb / ntot;

   % I believe that this will give "Leahy Normalization", i.e.,
   % Poisson noise = 2, in absence of deadtime effects.

   atrans = (atrans * len) / ntot;
   btrans = (btrans * len) / ntot;
   ctrans = (ctrans * len) / ntot;

   ctrans = ctrans[[0:alen-1]];

   if( donorm != 1 )
   {
      atrans = atrans * tbin*len / mcounta;
      btrans = btrans * tbin*len / mcountb;
      ctrans = ctrans * tbin*len / sqrt(mcounta*mcountb);
   }

   if( gsl_fft == 1)
   {
      atrans = atrans/len;
      btrans = btrans/len;
      ctrans = ctrans/len;
   }

   return freq, atrans, btrans, ctrans, ntot, mcounta, mcountb;
}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

public define sitar_avg_psd()
{
%     Take an evenly spaced lightcurve (presumed counts vs. time), and
%     calculate the PSD in segments of length l, averaged over the
%     whole lightcurve.  Segments with data gaps are skipped.
%     Deadtime corrections to Poisson noise are *not* made. Poisson
%     noise is *not* subtracted from average PSD.

   variable j,s,f,l,k,alen,ntrans,trans,strans,freq,atrans=0,ntot=0,mcount=0;
   variable dt, wdt, nwdt, wwdt, counts, len;

   variable tbin = qualifier("tbin",1.); % Some earlier versions used tbin, 
            tbin = qualifier("dt",tbin); %  allow but overide with dt
   variable bin_time = qualifier("times",NULL);
   variable donorm = qualifier("norm",1);

   switch(_NARGS)
   {
    case 2:
      (counts, len) = ();
   }
   {
      () = fprintf(fp,  "%s\n", %{{{
`
   (f,psd,navg,avg_cnts) = sitar_avg_psd(cnts,l [;dt=#,times=array,norm=#]);

     Take an evenly spaced lightcurve (presumed counts vs. time), and
     calculate the PSD in segments of length l, averaged over the
     whole lightcurve.  Segments with data gaps are skipped.  Deadtime
     corrections to Poisson noise are *not* made. Poisson noise is
     *not* subtracted from average PSD.

   Inputs:
     cnts    : Array of total counts in each time bin
     l       : Length of individual PSD segments (use a power of 2!!!)

   Qualifiers:
     dt      : Length of evenly spaced bins (default == 1)
     times   : Times of measurements (otherwise presumed to have no gaps)
     norm    : Determine PSD normalization convention.
                =1, 'Leahy normalization'. Poisson noise PSD == 2, in 
                    absence of deadtime effects).
               !=1, 'rms' or 'Belloni-Hasinger' or 'Miyamoto' normalization, 
                    i.e., PSD == (rms)^2/Hz & Poisson noise== 2/Rate.

   Outputs:
     f       : PSD frequencies ( == 1/Input Time Unit)
     psd     : Average PSD.
     navg    : Number of data segments going into the average
     avg_cnts: Average number of counts per segment of length l
`);   %}}}
      return;
   }

   len = typecast(len,Long_Type);

   alen = len/2;

   % Lengths of the counts and time arrays, and the number of
   % transforms we will expect to be able to do

   variable lcounts, ltime;

   lcounts = length(counts);

   if(bin_time != NULL)
   {
      ltime = length(bin_time);
   }
   else
   {
      bin_time = [0:lcounts-1]*tbin;
      ltime = lcounts;
   }

   if(lcounts != ltime ) 
   { 
      () = fprintf(fp, "\n Arrays have unequal lengths \n");
      return; 
   }

   if(lcounts < len) 
   { 
      () = fprintf(fp, "\n Lightcurve too short for desired segment size \n");
      return; 
   }

   dt = bin_time - shift(bin_time,-1);
   wdt = where( dt > 1.5*tbin );

   if(length(wdt) == 0)
   {
      ntrans = lcounts/len;

      _for j (1,ntrans,1)
      {
         s = (j-1)*len;
         f = j*len -1;

         trans = counts[[s:f]];
         strans = sum(trans);
         mcount = mcount + strans;
         
         trans = shift( fft(trans,1) , 1 );
         atrans = atrans + 2*( abs(trans[[0:alen-1]]) )^2 / strans;
 
         ntot = ntot + 1;
      }
   }
   else
   {
      nwdt = length(wdt) + 1;
      wwdt = [0,wdt,lcounts-1];

      _for k (1,nwdt,1)
      {
         ntrans = ( wwdt[k] - wwdt[k-1] )/len;
            
         _for l (1,ntrans,1)
         {
            s = (l-1) * len + wwdt[k-1];
            f = l * len + wwdt[k-1] - 1;

            trans = counts[[s:f]];
            strans = sum(trans);
            mcount = mcount + strans;

            trans = shift( fft(trans,1) , 1 );
            atrans = atrans + 2*( abs(trans[[0:alen-1]]) )^2 / strans;

            ntot = ntot + 1;
         }
      }
   }

   if( ntot == 0 )
   {
      () = fprintf(fp, "\n No transforms performed! \n");
      return;
   }

   freq = [1:alen] / (len * tbin);
   mcount = mcount / ntot;

   % I believe that this will give "Leahy Normalization", i.e.,
   % Poisson noise = 2, *** in the absence of deadtime effects ***

   atrans = (atrans * len) / ntot;

   if( donorm != 1 )
   {
      atrans = atrans * tbin*len / mcount;
   }

   if( gsl_fft == 1)
   {
      atrans = atrans/len;
   }

   return freq, atrans, ntot, mcount;
}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

public define sitar_lags()
{
%     Given two input PSD (without noise subtracted), their CPD, and
%     their associated noise levels, all as functions of Fourier
%     frequency, calculate the time lag and coherence function (and
%     associated errors) vs. f
%
%     For background discussion, see:
%
%        Vaughan, B. A. & Nowak, M. A.  1997, ApJ, 474, L43
%
%     and references therein (especially Bendat & Piersol 1986)
%
   variable freq,psda,psdb,cpd,noisea,noiseb,navg=1;

   switch(_NARGS)
   {
    case 6:
      (freq,psda,psdb,cpd,noisea,noiseb) = ();
   }
   {
    case 7:
      (freq,psda,psdb,cpd,noisea,noiseb,navg) = ();
   }
   {
      () = fprintf(fp,  "%s\n", %{{{
`
   (lag,dlag,g,dg) = sitar_lags(freq,psda,psdb,cpd,noisea,noiseb [,navg]);

     Given two input PSD (without noise subtracted), their CPD, and
     their associated noise levels, all as functions of Fourier
     frequency, calculate the time lag and coherence function (and
     associated errors) vs. f

   Inputs:
     freq    : Fourier frequency array
     psda/b  : Power spectra array associated with two lightcurves
               (Noise subtraction *not* applied!)
     cpd     : Cross power spectra array associated with two lightcurves
     noisea/b: Power spectra noise level. Array *or* a constant.

   Optional Unput:
     navg    : Number of averages (over independent data segments and frequency
               bins) associated with each input frequency bin. Can be a constant
               vs. f (defaulted to 1).

   Outputs:
     lag     : Time lag vs. frequency. Negative values mean psdb lags psda
     dlag    : Associated error on the lag.
     g       : Coherence function
     dg      : Error on the coherence function
`);   %}}}
      return;
   }  

   variable lfreq = length(freq), lnoisea = length(noisea), lnoiseb = length(noiseb),
            lnavg = length(navg);

   if( lnavg==1 )
   {
      navg = ones(lfreq)*navg;
      lnavg = lfreq;
   }

   if( lnoisea==1 )
   {
      noisea = ones(lfreq)*noisea;
      lnoisea = lfreq;
   }

   if( lnoiseb==1 )
   {
      noiseb = ones(lfreq)*noiseb;
      lnoiseb = lfreq;
   }

   if( lnavg == 1 )
   {
      navg = ones(lfreq)*navg;
   }

   variable lpsda = length(psda), lpsdb = length(psdb), lcpd = length(cpd);

   if( lnoisea != lnoiseb or lnoiseb != lpsda or lpsda != lpsdb or
       lpsdb != lcpd or lcpd != lnavg or lnavg != lfreq )
   {
      () = fprintf(fp, "\n Inconsistent array sizes \n");
      return;
   }

   variable sa, sb, ns, gs, dgs, gerr, lag, dlag, acpd = abs(cpd);

   % Define noise subtracted PSDs

   sa = psda - noisea;
   sb = psdb - noiseb;

   ns = ( sa*noiseb + sb*noisea + noisea*noiseb ) / navg;

   % Estimate of intrinsic coherence, high signal to noise

   gs = ( acpd^2 - ns ) / sa / sb;

   % Error on intrinsic coherence, squared, sans Poisson noise

   dgs = 2. * ( 1 - gs )^2 / gs^2;

   % Overall error estimate

   gerr = sqrt( ( 2. * ns^2 * navg / ( acpd^2 - ns )^2 +
                  ((noisea+noiseb)/2)^2 * (1/sa^2 + 1/sb^2) + dgs )
                / navg );

   % Calculate the magnitude of the phase lag
   % (Note: PI is a s-lang intrinsic)

   lag = Real(cpd) / acpd;
   lag = acos( lag ) / 2 / PI / freq * Imag(cpd)/abs( Imag(cpd) );

   % Phase lag errors

   dlag = acpd^2 / psda / psdb;
   dlag = sqrt(1-dlag) / sqrt( dlag * 2 * navg ) / 2 / PI / freq;

   return lag, dlag, gs, dgs;
}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

public define sitar_lbin_cpd()
{
%     Logarithmically rebin two PSDs (and, optionally, their
%     associated Poisson, noise levels) and the associated CPD

   variable freq, psda, psdb, cpd, dff, noisea = NULL, noiseb = NULL, rrev = NULL;

   switch(_NARGS)
   {
    case 5:
      (freq, psda, psdb, cpd, dff) = ();
   }
   {
    case 6:
      (freq, psda, psdb, cpd, dff, noisea) = ();
   }
   {
    case 7:
      (freq, psda, psdb, cpd, dff, noisea, noiseb) = ();
   }
   {
    case 8:
      (freq, psda, psdb, cpd, dff, noisea, noiseb, rrev) = ();
   }
   {
      () = fprintf(fp,  "%s\n", %{{{
`
   (aflo,afhi,apsda,apsdb,acpd,nf,[anoisea,anoiseb]) = 
     sitar_lbin_cpd(f,psda,psdb,cpd,dff,[noisea,noiseb,&rev]);

     Logarithmically rebin two PSDs (and, optionally, their associated
     Poisson noise levels) and the associated Cross Power Spectral
     Density

   Inputs:
     f        : Array of Fourier frequencies
     psda/b   : Arrays of PSD values
     dff      : \Delta f/f value for logarithmic bin spacing

   Optional Inputs:
     noisea/b : Arrays (*or single values*) of Poisson noise levels
     rev      : If declared (i.e., 'isis> variable rev;') and input,
                rev returns the reverse inidices for the binning (i.e.,
                (apsda[i] = mean( psda[ rev[i] ] ), etc.)

   Outputs:
     aflo     : Lower boundary of Fourier frequency bin
     afhi     : Upper boundary of Fourier frequency bin
     apsda/b  : Rebinned Power Spectral Density values
     acpd     : Rebinned Cross Power Spectral Density values
     nf       : Number of frequencies going into the bin
     anoisea/b: Rebinned noise level (array, even for single value input)
`);  %}}}
      return;
   }

   variable lf, lpsda, lpsdb, lcpd, lnoisea, lnoiseb;
   variable fmin, fmax, nmax, flo, fhi, iw, rev;
   variable nfreq, afreq, afreqlo, afreqhi, apsda, apsdb, acpd, 
            anoisea, anoiseb;

   lf = length(freq);
   lpsda = length(psda);
   lpsdb = length(psdb);
   lcpd = length(cpd);

   if( noisea != NULL )
   {
      lnoisea = length(noisea);
      if( lnoisea == 1 )
      {
         noisea = ones(lf) * noisea;
      }
      else if( lnoisea != lpsda )
      {
         () = fprintf(fp, "\n Noise A and PSD arrays have unequal lengths \n");
         return; 
      }
   }

   if( noiseb != NULL )
   {
      lnoiseb = length(noisea);
      if( lnoisea == 1 )
      {
         noiseb = ones(lf) * noiseb;
      }
      else if( lnoiseb != lpsdb )
      {
         () = fprintf(fp, "\n Noise B and PSD arrays have unequal lengths \n");
         return; 
      }
   }

   if( lf != lpsda or lpsda != lpsdb or lpsdb != lcpd )
   {
      () = fprintf(fp, "\n Inconsistent Array Sizes \n");
      return;
   }

   fmin = min(freq);
   fmax = max(freq);

   nmax = typecast((log10(fmax)-log10(fmin))/log10(1+dff)+0.5,Integer_Type);

   flo = fmin * 10.^([0:nmax-1]*log10(1+dff));
   fhi = make_hi_grid(flo);

   nfreq = histogram(freq, flo, fhi, &rev);

   iw = where(nfreq !=0);

   variable liw = length(iw);

   nfreq = nfreq[iw];
   afreqlo = Double_Type[liw];
   apsda = @afreqlo;
   apsdb = @afreqlo;
   acpd  = Complex_Type[liw];

   rev = rev[iw];

   variable i;

   if( noisea == NULL and noiseb == NULL)
   {
      _for i (0, liw-1, 1)
      {
         afreqlo[i] = min(freq[rev[i]]);
         apsda[i] = sum( psda[rev[i]] );
         apsdb[i] = sum( psdb[rev[i]] );
         acpd[i] = sum( cpd[rev[i]] );
      }
   } 
   else if( noiseb == NULL )
   {
      anoisea = @afreqlo;
      anoiseb = NULL;

      _for i (0, liw-1, 1)
      {
         afreqlo[i] = min(freq[rev[i]]);
         apsda[i] = sum( psda[rev[i]] );
         apsdb[i] = sum( psdb[rev[i]] );
         acpd[i] = sum( cpd[rev[i]] );
         anoisea[i] = sum( noisea[rev[i]] );
      }

      anoisea = anoisea / nfreq;
   }
   else if( noisea == NULL )
   {
      anoiseb = @afreqlo;
      anoisea = NULL;

      _for i (0, liw-1, 1)
      {
         afreqlo[i] = min(freq[rev[i]]);
         apsda[i] = sum( psda[rev[i]] );
         apsdb[i] = sum( psdb[rev[i]] );
         acpd[i] = sum( cpd[rev[i]] );
         anoiseb[i] = sum( noiseb[rev[i]] );
      }

      anoiseb = anoiseb / nfreq;
   }   
   else
   {
      anoisea = @afreqlo;
      anoiseb = @afreqlo;

      _for i (0, liw-1, 1)
      {
         afreqlo[i] = min(freq[rev[i]]);
         apsda[i] = sum( psda[rev[i]] );
         apsdb[i] = sum( psdb[rev[i]] );
         acpd[i] = sum( cpd[rev[i]] );
         anoisea[i] = sum( noisea[rev[i]] );
         anoiseb[i] = sum( noiseb[rev[i]] );
      }

      anoisea = anoisea / nfreq;
      anoiseb = anoiseb / nfreq;
   }   


   afreqhi = make_hi_grid(afreqlo);
   afreqhi[liw-1] = max(freq[rev[liw-1]]);

   apsda = apsda / nfreq;
   apsdb = apsdb / nfreq;
   acpd = acpd / nfreq;

   if(rrev != NULL)
   {
     @rrev = rev;
   }

   if( noisea == NULL and noiseb == NULL )
   {
      return afreqlo, afreqhi, apsda, apsdb, acpd, nfreq;
   } 
   else
   {
      return afreqlo, afreqhi, apsda, apsdb, acpd, nfreq, 
             anoisea, anoiseb;
   }
}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

public define sitar_lbin_psd()
{
%     Logarithmically rebin a PSD (and, optionally, its associated
%     Poisson noise)

   variable freq, psd, dff, noise = NULL, rrev = NULL;

   switch(_NARGS)
   {
    case 3:
      (freq, psd, dff) = ();
   }
   {
    case 4:
      (freq, psd, dff, noise) = ();
   }
   {
     case 5:
      (freq, psd, dff, noise, rrev) = ();
   }
   {
      () = fprintf(fp, "%s\n", %{{{
`
   (aflo,afhi,apsd,nf,[anoise]) = sitar_lbin_psd(f,psd,dff [,noise,&rev]);

     Logarithmically rebin a PSD (and, optionally, its associated Poisson noise)

   Inputs:
     f      : Array of Fourier frequencies
     psd    : Array of PSD values
     dff    : \Delta f/f value for logarithmic bin spacing 

   Optional Inputs:
     noise  : Array (*or single value*) of Poisson noise levels
     rev    : If declared (i.e., 'isis> variable rev;') and input,
              rev returns the reverse inidices for the binning (i.e.,
              (af[i] = mean( freq[ rev[i] ] ), etc.)

   Outputs:
     aflo   : Lower bounds of Fourier frequency bins
     afhi   : Upper bounds of Fourier frequency bins
     apsd   : Rebinned PSD values
     nf     : Number of frequencies going into the bin
     anoise : Rebinned noise level (array, even for single value input)
`);   %}}}
      return;
   }

   variable lf, lpsd, lnoise, fmin, fmax, nmax, flo, fhi, iw, rev;
   variable nfreq, lnf, afreq, afreqlo, afreqhi, apsd, anoise;

   lf = length(freq);
   lpsd = length(psd);

   if(lf != lpsd) 
   { 
      () = fprintf(fp, "\n Frequency and PSD arrays have unequal lengths \n");
      return; 
   }

   if( noise != NULL )
   {
      lnoise = length(noise);
      if( lnoise == 1 )
      {
         noise = ones(lpsd)*noise;
      }
      else if( lnoise != lpsd )
      {
         () = fprintf(fp, "\n Noise and PSD arrays have unequal lengths \n");
         return; 
      }
   }

   fmin = min(freq);
   fmax = max(freq);

   nmax = typecast((log10(fmax)-log10(fmin))/log10(1+dff)+0.5,Integer_Type);

   flo = fmin * 10.^([0:nmax-1]*log10(1+dff));
   fhi = make_hi_grid(flo);

   nfreq = histogram(freq, flo, fhi, &rev);

   iw = where(nfreq !=0);

   nfreq = nfreq[iw];
   lnf = length(nfreq);

   apsd = Double_Type[lnf];
   afreqlo = @apsd;
   afreqhi = @apsd;

   rev = rev[iw];

   variable i;

   if( noise == NULL )
   {
      _for i (0, lnf-1, 1)
      {
         afreqlo[i] = min( freq[rev[i]]);
         apsd[i] = sum( psd[rev[i]] );
      }
   } 
   else
   {
      anoise = @apsd;

      _for i (0, lnf-1, 1)
      {
         afreqlo[i] = min( freq[rev[i]]);
         apsd[i] = sum( psd[rev[i]] );
         anoise[i] = sum( noise[rev[i]] );
      }

      anoise = anoise / nfreq;
   }

   afreqhi = make_hi_grid(afreqlo);
   afreqhi[lnf-1] = max(freq[rev[lnf-1]]);

   apsd = apsd / nfreq;

    if(rrev != NULL)
    {
      @rrev = rev;
    }

   if( noise == NULL )
   {
      return afreqlo, afreqhi, apsd, nfreq;
   } 
   else
   {
      return afreqlo, afreqhi, apsd, nfreq, anoise;
   }
}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

public define sitar_define_psd()
{
%     Define a PSD to be a fittable dataset in ISIS, with frequency
%     values set to keV.

   variable f_lo, f_hi, psd, err, noise=NULL;

   switch(_NARGS)
   {
    case 4:
      (f_lo, f_hi, psd, err) = ();
   }
   {
    case 5:
      (f_lo, f_hi, psd, err, noise) = ();
   }   
   {
      () = fprintf(fp,  "%s\n", %{{{
`
   id = sitar_define_psd(flo,fhi,psd,err [,noise]; noise=# or [Array]);

     Define psd to be an ISIS counts data set.

   Inputs:
     f_lo   : Low value of PSD frequency bin (Hz -> keV)
     f_hi   : High value of PSD frequency bin (Hz -> keV)
     psd    : Array of PSD values (Power -> Counts/bin)
     err    : Array of PSD Errors

   Optional Input:
     noise  : Array (*or single value*) of Poisson noise levels.
              If undefined, background is undefined.

   Qualifiers:
     noise  : Same as above.  Supercedes optional input.

   Outputs:
     id     : Dataset Index
`);
      return;
   }

   variable id;

   noise = qualifier("noise",noise);

   if(noise !=NULL)
   {
      variable lnoise = length(noise); 
      if( lnoise == 1 )
      {
         noise = ones(length(psd))*noise;
         id = define_counts(_A(f_hi),_A(f_lo),reverse(psd-noise),
                            reverse(err));
      }
      else if(lnoise != length(psd))
      {
         print("Noise array length differs from psd length");
         return;
      }
      else
      {
         id = define_counts(_A(f_hi),_A(f_lo),reverse(psd-noise),
                            reverse(err));
      }
   }
   else
   {
      id = define_counts(_A(f_hi),_A(f_lo),reverse(psd),reverse(err));
   }
   return id;
}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% The GSL library has to be imported in order to run this routine, so,
% do a check, and if it exists, load the Gregory-Loredo function.

#ifexists interp_akima

% t = event times; tmin & tmax =  max & min time of lightcurve,
% mmax = maximum partition number in trial lightcurves
% ta, adt = effective area/deadtime times and values
% dodither = Toggle for applying dither correction to lightcurve
% thresh = helps determine truncation of partitioned lightcurve

public define sitar_glvary()
{
   variable t, dodither=0; 

   switch(_NARGS)
   {
    case 1:
      t = ();
   }
   {
      () = fprintf(fp, "%s\n", %{{{
`
   gl = sitar_glvary(t [;tmin=#, tmax=#, mmax=#, thresh=#,
                         nbins=#, texp=array, frac_exp=array]);

     Use the Gregory-Loredo algorithm to find the odds ratios that
     even divisions of a lightcurve are better descriptions of the
     data than a constant lightcurve.  A 'best estimate' lightcurve
     can also be output for a lightcurve with nbins.

   Inputs:
     t      : Array of event times

   Qualifiers:
     tmin     : The minimum time to consider (Default=min(t)-1)
     tmax     : The maximum time to consider  (Default=max(t)+1)
     mmax     : Consider lightcurve divisions from 2 to mmax evenly spaced 
                bins. (Default=300)
     thresh   : Truncate the maximum number of partitions considered by ignoring
                partitionings for which--
                   \sum_m(odds ratio) < max(\sum_n(odds ratio))/exp(thresh) , 
                where n = 2 -> mmax.  I.e., the more partitionings you have,
                the less significant the results tend to become.  thresh
                sets a minimum probability: (1-p_min) ~ exp(thresh)*(1-p_peak).
                (Default=2.)
     nbins    : Create an output lightcurve with nbins. (Default=mmax)
     texp,    : A pair of arrays that give values for the fractional exposure
     frac_exp   as a function of time.  *** Only non-zero values of exposure 
                will be retained ***, which will then be interpolated and used
                to correct the lightcurve rates. These arrays must have a
                minimum of five entries. (Default is for no correction.)

   Outputs:
     gl.p         : Total probability that some evenly partitioned lightcurve, with up
                    to gl.mmax bins, is a better description than a constant lightcurve
     gl.ppart     : The probability for an individual evenly partitioned lightcurve
                    that it is a better description than a constant lightcurve
     gl.lodds_sum : The natural logarithm of the sum of the odds ratios comparing 
                    lightcurves with two or more partitions to a constant lightcurve.
                    gl.p == exp(gl.lodds_sum)/[1+exp(gl.lodds_sum)]
     gl.mpeak     : The number of partitions for the evenly partitioned lightcurve 
                    with the maximum probability.
     gl.mmax      : The maximum number of partitions actually used (influenced by
                    the setting of the thresh parameter)
     gl.m         : The number of partitions corresponding to each evenly partitioned
                    lightcurve considered (=[2:mmax])
     gl.pm        : Total probability that some evenly partitioned lightcurve
                    is a better description than a constant lightcurve for each
                    maximum number of partitions considered ([2:mmax])
     gl.nj        : The counts histogram corresponding to each partitioning above
     gl.aj        : The integrated fractional exposure for each partitioning above
     gl.a_avg     : The averaged fractional exposure.
     gl.tmin      : The value of tmin actually used (maximum of [tmin,min(texp)]);
     gl.tmax      : The value of tmax actually used (minimum of [tmax,max(texp)]);
     gl.tlc       : The output lightcurve times (an array with input nbins bins)
     gl.rate      : Best estimate of the lightcuve rates at the above times
     gl.erate     : Best estimate of the lightcuve rate errors at the above times
`);   %}}}
      return;
}

   variable tmin     = qualifier("tmin",min(t)-1.);
   variable tmax     = qualifier("tmax",max(t)+1.);
   variable mmax     = qualifier("mmax", 300);
   variable thresh   = qualifier("thresh",2);
   variable ta       = qualifier("texp",NULL);
   variable adt      = qualifier("frac_exp",NULL);
   variable nbins    = qualifier("nbins",mmax);

   if(length(ta)==length(adt) && length(adt) >=5) dodither=1;

   if(dodither)
   {
      variable inz = where(adt>0);
      variable adt_use = [inz[0]:inz[length(inz)-1]];

      ta = ta[adt_use];
      adt = adt[adt_use];

      variable na = length(adt);
      if(na < 5)
      {
         () = fprintf(fp,"\n Fractional exposure array is too short. \n");
         return;
      }

      if(tmin < ta[0]){ tmin = ta[0]; }
      if(tmax > ta[na-1]){ tmax = ta[na-1]; }
   }

   % Only look at times within tmin & tmax

   t = t[where(t>=tmin and t<tmax)];

   variable j;
   variable n=length(t);
   variable dt_int = tmax - tmin;
   variable a_avg=0.;

   if(dodither)
   {
      % Define the spline of the effective area curve

      variable spline_adt = interp_akima_init(ta,adt);
      a_avg = log(interp_eval_integ(spline_adt,tmin,tmax)/dt_int);

      % Do integration of spline in neighboring points, then sum

      variable iadt = Double_Type[na];
      _for j (1,na-1,1)
      {
         iadt[j] = interp_eval_integ(spline_adt,ta[j-1],ta[j]);
      }
      iadt=cumsum(iadt);

      % Define the spline of the integrated effective area curve

      variable spline_iadt = interp_akima_init(ta,iadt);
   }

   variable nj = Array_Type[int(mmax)-1];
   variable aj = @nj;
   variable m=[2:int(mmax):1];
   variable lods=Double_Type[int(mmax)-1];

   variable lo,hi,im;

   _for im (2,int(mmax),1)
   {
      % Create grid for partitioning lightcurve into im bins

       lo = [0:im-1];
       hi = [1:im];
       lo = __tmp(lo)*(dt_int/im);
       hi = __tmp(hi)*(dt_int/im);

      if(dodither)
      {
         % For each partitioning, create average deadtime/effective
         % area per bin.  Dead bins are set equal to the average 
         % effective area, so as not to contribute to the sum 

         aj[im-2] = (interp_eval(spline_iadt,hi)
                    -interp_eval(spline_iadt,lo))*(im/dt_int);

         aj[im-2][where(aj[im-2]==0.)] = a_avg;
      }
      else
      {
         aj[im-2] = Double_Type[im]+1;
      }

      % For each lightcurve partition, the arrays of counts per bin

      nj[im-2] = histogram(t,lo,hi);

      % The odds ratio for each partitioning.  The nj[im-2]*()
      % term is removed if effective area variation is unimportant

      if(dodither)
      {
         lods[im-2] = sum( nj[im-2]*(a_avg-log(aj[im-2])) + 
                           lngamma(nj[im-2]+1) ) +   
                      n*log(im) + lngamma(im) - lngamma(n+im);
      }
      else 
      {
         lods[im-2] = sum( lngamma(nj[im-2]+1) ) +
                      n*log(im) + lngamma(im) - lngamma(n+im);
      }
   }

   % This bit uses the return pieces from above to find what the
   % maximum odds ratio is, and temporarily takes that out (lomax),
   % and then truncates the number of lightcurve binnings kept

   variable iw = where(lods==max(lods));
   variable lomax = lods[min(iw)];

   lods = exp(lods-lomax);
   variable csum = cumsum(lods)/(m-1);

   % Truncate the number of partitionings of the lightcurve

   if(thresh < 0) thresh==2.;
   variable imax = max(where(csum >= max(csum)/exp(thresh)));

   % Return the probability (p), the odds ratio for each partioning
   % (ppart), the log of the summed odds (lodds_sum), the # of
   % partitions for the peak odds ratio (mpeak), the # of partitions
   % for each (m), the histogrammed counts for each partitioning (nj),
   % the integrated and averaged effective area for each partioning 
   %(aj, a_avg), and the used min/max times (tmin/tmax)

   variable gl_struct=   
     struct{p, ppart, lodds_sum, mpeak, mmax, pm, m, nj,
            aj, a_avg, tmin, tmax, tlc, rate, erate};

   % The # of partitions of the lightcurve with the highest odds ratio.

   gl_struct.mpeak = min(iw)+2;

   % Calculate the total probability of variability (p) and the
   % probability for each partitioning of the lightcurve (ppart).

   variable msum = csum[imax];
   variable lsum = log(msum) + lomax;
   gl_struct.p = 1/(exp(-lsum)+1);
   gl_struct.ppart = __tmp(lods)[[0:imax]]/
                            ((imax+1)*(msum+exp(-lomax)));
   gl_struct.pm = 1/(exp(-lomax)/csum+1);

   % The rest of the return values

   gl_struct.lodds_sum = lsum;
   gl_struct.mmax = imax+2;
   gl_struct.m = __tmp(m);
   gl_struct.nj = __tmp(nj);
   gl_struct.aj = __tmp(aj);
   gl_struct.a_avg = exp(a_avg);
   gl_struct.tmin = tmin;
   gl_struct.tmax = tmax;

   % Return a lightcurve with nbin bins

   variable tfrac = [0:nbins-1]/(nbins*1.);
   variable rate = Double_Type[nbins];
   variable erate=@rate;

   % We might not have used all the data ...

   variable ntot=sum(gl_struct.nj[0]); 

   variable nj_ii_k, oratio_ii, aj_ii_k, drate;

   % Best estimate of the lightcurve is calculated here.
   % Fractional area/deadtime correction is included.

   variable ii;
   _for ii (0,imax,1)
   {
      variable mloop=ii+2;
      variable k = int( tfrac*mloop );

      nj_ii_k=gl_struct.nj[ii][k]; 
      oratio_ii=gl_struct.ppart[ii];
      aj_ii_k = gl_struct.aj[ii][k];
      drate = (mloop/(ntot+mloop))*
              (__tmp(oratio_ii)*(nj_ii_k+1)/aj_ii_k);
      rate = __tmp(rate) + drate; 
      erate = __tmp(erate) + (mloop/(ntot+mloop+1))*
              (__tmp(drate)*(__tmp(nj_ii_k)+2)/__tmp(aj_ii_k));
   }

   % Here I differ from Gregory & Loredo, who only include the
   % variable part of the lightcurve (i.e., partitionings with >=2
   % bins), reasonable for p~1. But p can be lower, therefore the
   % estimated constant lightcurve part should be included.

   rate = __tmp(rate) + (1.-gl_struct.p)/gl_struct.a_avg;
   erate = __tmp(erate) + (1.-gl_struct.p)/gl_struct.a_avg^2;
   erate = sqrt(__tmp(erate)-rate^2);
 
   gl_struct.tlc = __tmp(tfrac)*(tmax-tmin)+tmin;
   gl_struct.rate = __tmp(rate)*(ntot/(tmax-tmin));
   gl_struct.erate = __tmp(erate)*(ntot/(tmax-tmin));

   return gl_struct; 
}
#endif

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


private define reb_evt_use()
{
   () = fprintf(fp,  "%s\n", %{{{
`
   lc = sitar_bin_events(t,dt,gti_lo,gti_hi [;tstart=#,tstop=#,obin=#,
                                              user_bin=array,delgap]);
    
   Inputs:
     t          : Discrete times at which rates are measured
     dt         : Width of evenly spaced time bins in rebinned lightcurve
                  (Ignored if a user_bin is input.)

   Qualifiers:
     tstart     : Beginning of first output time bin (default: min(t))
     tstop      : End of last output time bin (default: max(t))
     obin       : For a tstop, last bin width might be < dt.  Include this
                  "overflow" bin if width is > obin*dt (default: obin=0.1)
     user_bin.lo: Lower time bounds for user defined grid
     user_bin.hi: Upper time bounds for a user defined grid
                  (Overrides dt, tstart, and tstop inputs)
     minexp     : Do not include bins with exposure < minexp (default: 0)
     delgap     : If set, delete 0 exposure bins from the lightcurve,
                  otherwise set their rate & error to 0 (default)
    
   Outputs:
     lc.rate    : Mean rate of resulting binning
     lc.err     : Error on bin rate (minimum = -log(0.32)/(bin exposure))
     lc.cts     : Counts in a bin
     lc.bin_lo  : Lower time bounds of binned lightcurve
     lc.bin_hi  : Upper time bounds of binned lightcurve
     lc.expos   : Exposure of a time bin (i.e., intersection with
                  good time intervals)
`);  %}}}    
   return;
}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

public define sitar_bin_events( )
{
   % Bin a rate lightcurve to a new scale

   variable t, dt, gti_lo, gti_hi; 
   variable custf=0, istck;
   variable lc = struct{ rate, err, cts, bin_lo, bin_hi, expos };

   variable tstart = qualifier("tstart",NULL);
   variable tstop  = qualifier("tstop",NULL);
   variable bintol = qualifier("obin",0.1);
   variable minexp = qualifier("minexp",0);
   variable cust   = qualifier("user_bin",NULL);
   variable gap    = qualifier_exists("delgap");

   switch (_NARGS)
   {
    case 4:
      ( t, dt, gti_lo, gti_hi ) = ( );
   }
   {
      reb_evt_use();
      _pop_n(_NARGS);
      return;
   }

   if ( cust != NULL )
   {
      if ( is_struct_type(cust) == 1 && struct_field_exists(cust,"lo")
                                     && struct_field_exists(cust,"hi") )
      { 
         custf = 1;
         dt = 0;
      }
      else 
      {
         reb_evt_use();
         () = fprintf(fp,  "%s\n", %{{{
`
  Custom user_bin is incorrect.  Input structure requires: 
     variable user_bin = struct{ lo, hi };
`); %}}}
         return;
      }
   }

   if ( dt == NULL ) dt = 0;

   variable mint = min(t);
   variable maxt = max(t);

   if ( tstart == NULL ) tstart = mint;
   if ( tstop == NULL )  tstop = maxt;

   if( maxt < tstart || mint > tstop)
   {
      reb_evt_use( );
      () = fprintf(fp, "\n Data outside of bounds of start & stop times:\n");
      () = fprintf(fp, "      [tstart,tstop] = [%e,%e]\n",tstart,tstop);
      () = fprintf(fp, "      [min_t, max_t] = [%e,%e]\n\n",mint,maxt);
      return;
   }

   % Remove tstart from all the times to begin with, add back in at the end
   t -= tstart;
   gti_lo -= tstart;
   gti_hi -= tstart;

   if( custf == 0 && dt > 0)
   {
      variable nfrac = ( tstop - tstart )/dt;
      variable nbins = int(nfrac);
      lc.bin_lo = [0:nbins-1]*dt;
      lc.bin_hi = [1:nbins]*dt;
      if( nfrac-nbins > bintol )
      {
         lc.bin_lo = [lc.bin_lo,lc.bin_hi[nbins-1]];
         lc.bin_hi = [lc.bin_hi,tstop-tstart]; 
      }
   }
   else if ( custf )
   {
      lc.bin_lo = @cust.lo-tstart;
      lc.bin_hi = @cust.hi-tstart;
   }
   else
   {
      reb_evt_use();
      () = fprintf(fp, "\n Time bin or custom bin array must be input.\n\n");
      return;
   }

   % Create the histogram of the times to creat counts per bin
   lc.cts = histogram(t, lc.bin_lo, lc.bin_hi);
   lc.rate = 0.*lc.cts;
   lc.err = @lc.rate;
   lc.expos = @lc.rate;

   % Figure out the bins within the Good Time Intervals

   variable i,j,ilo,ihi;
   _for i (0,length(gti_lo)-1,1)
   {
      ilo=where(lc.bin_hi>gti_lo[i]);
      if(length(ilo))
      {
         ihi=where(lc.bin_lo<gti_hi[i]);
         if(length(ihi))
         {
            ilo=min(ilo);
            ihi=max(ihi);

            _for j (ilo+1,ihi-1,1)
            {
               lc.expos[j] += lc.bin_hi[j]-lc.bin_lo[j];
            }
            if(ilo==ihi)
            {
               lc.expos[ilo] += min([gti_hi[i],lc.bin_hi[ihi]]) -
                                max([gti_lo[i],lc.bin_lo[ilo]]);
            }
            else
            {
               lc.expos[ilo] += lc.bin_hi[ilo]-max([gti_lo[i],lc.bin_lo[ilo]]);
               lc.expos[ihi] += min([gti_hi[i],lc.bin_hi[ihi]])-lc.bin_lo[ihi];
            }
	 }
      }
   }

   % Figure out which exposure bins aren't empty, and make rates & errors

   variable idx = where( lc.expos > minexp);
   if ( length(idx) )
   {
      lc.rate[idx] = lc.cts[idx]/lc.expos[idx];
      lc.err[idx] =  sqrt(abs(lc.cts[idx]+log(0.32)/2)-log(0.32)/2)/lc.expos[idx];
   }

   lc.bin_lo += tstart;
   lc.bin_hi +=tstart;

   if (gap)
   {
      lc.cts = lc.cts[idx];
      lc.rate = lc.rate[idx];
      lc.err = lc.err[idx];
      lc.expos = lc.expos[idx];
      lc.bin_lo = lc.bin_lo[idx];
      lc.bin_hi = lc.bin_hi[idx];
   }

   return lc;
}
