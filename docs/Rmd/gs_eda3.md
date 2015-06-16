


##
##
## Growing season EDA: Space-Time

The plots below explore patterns in TDD through space and time.



### TDD spatial distributions by ecoregion and year | threshold

This is similar to the first plot of TDD pdfs by ecoregion and threshold
except that there is not integration across time.
Separate, overlapping annual density curves are shown.
This plot shows annual pdfs specifically for the 10% TDD variable.
Overlapping lines represent pdfs for individual years, 1982 - 2010.

![](gs_eda3_files/figure-html/plot_tdd_02a-1.png) 

The above plot shows clear inter-annual differences in TDD distribution shape, scale, and location,
but it is difficult to parse visually.
The following plot shows annual spatial distributions as boxplots with superimposed point observations,
ordered by year.
This plot reveals clearly not only the differences between ecoregions suggested by the previous plot,
but also annual changes in TDD distributions within ecoregions.

![](gs_eda3_files/figure-html/plot_tdd_02b-1.png) 

### TDD spatial distributions by year | threshold and ecoregion

In the next plot, two ecoregion cases are highlighted for direct comparison of annual spatial TDD distributions for each TDD threshold.
Note the substantial difference in variability between ecoregions, independent of time,
in addition to differences in means and how these may appear to correlate during some sub-periods and not others.

![](gs_eda3_files/figure-html/plot_tdd_02c-1.png) 

### Time series by sampled locations | threshold and ecoregion

Continuing with the comparison of these two regions,
rather than showing spatial distributional summaries as before,
the time series of TDD below is now broken out by lines representing the continuity of TDD through time at sampled pixels within each ecoregion.

![](gs_eda3_files/figure-html/plot_tdd_02d-1.png) 

### Spatial standard deviation by ecoregion and year

Focusing again specifically on the 10% TDD threshold,
this heat map displays the standard deviation in 10% TDD by ecoregion and year.

![](gs_eda3_files/figure-html/plot_tdd_02e-1.png) 
