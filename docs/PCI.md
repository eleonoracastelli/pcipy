# PCI
(from John's AAS poster)

Principal Component Interferometry (PCI) is an analysis technique whose theory we outline in two papers:
- [first](https://journals.aps.org/prd/abstract/10.1103/PhysRevD.103.042006) 
- [second](https://journals.aps.org/prd/abstract/10.1103/PhysRevD.104.122001)

PCI is based on Principal Component Analysis (PCA):

```{eval-rst}
.. note::

    PCA is a data science workhorse: it recasts a data matrix into independent components ranked by their variance. Typical applications include data compression by restricting to the high-variance components.
```  

## Alternatives to TDI
PCI is presented here as an alternative approach to TDI (Time Delay Interferometry) to suppress laser noise in the LISA measurement chain.

TDI has been demonstrated to work, but other alternative approaches might prove to be useful to deal with some of the various scenarios for LISA. In particular, other methods of working around the loud laser noise in the phasemeter channels may have advantages. In particular, TDI required precise auxiliary information about the constellation to work correctly.

Possibly a mix of method may eventually be part of the analysis and diagnostic toolkit for LISA operations.

A handful of alternative approaches have been explored at some level:

| Technique  |  Description |
|---|---|      
| TDI  | Standard approach to LISA initial noise reduction. current baseline model for LISA <br/>ground segment design and specifically for the level 1 data products. There is a <br/>series of prescriptions which can accomodate the effects of constellation <br/>dynamics to different derivative orders.|
| aPCI  | Our method uses no explicit information about the constellation structure or <br/>motion, identifying the laser noise correlations from their dominant imprint <br/>on the data. |
| TDI-Ranging | this variant following TDI, expect that the precise time delays ar  |
| TDI-infinity  | this approach encodes the information about the laser noise propagation delays <br/>in a discretely time-dependent design matrix, then computationally identifies <br/>channels free of laser noise with any approximations for the constellation's motion.|

## How does PCI work?

Since it's an alternative to TDI, we expect to obtain the same results as for TDI.

### What does TDI do?

Each TDI channel $A_\alpha$ is produced by some linear combination of various delayed measurements

$$
 A_\alpha (t) = \sum_{k=1}^6 \sum_{n=1}^{n_\text{max}} c_{kn\alpha}(t)\mathcal{D}[d_{kn\alpha}(t)] y_k(t),
$$

where $c_{kn\alpha}(t)$ are the coefficients, $\mathcal{D}$ is the delay operator, and $y_k(t)$ are the interferometer measurements (maybe phasemeter?)

In practice, we have discretely sampled data. Thus, a delay is a fractional delay filter of some finite half-width $n_h$

$$
A_{i\alpha} = \sum_{l=-n_h}^{n_h} \sum_{k=1}^6 \sum_{n=1}^{n_\text{max}} c_{kn\alpha}(t_i)f_{lkn}\alpha(t_i)\left[\mathcal{D}y_k\right](t_i)
$$

Gathering all the integer-time shifted data into a matrix $X$, we can write the previous equation as
 
$$
A_{i\alpha} = X_i g_\alpha (t_i),
$$
where $X_i$ is the matrix of shifted observations, containing all the time-shifted versions of $y$, and  $g_\alpha (t_i)$ is the vectors of coefficients.

data driven principal component interferometri (PCI) means that we're deriving alternatives of $A_alpha$ by exploring general combinations of $X$.

But how you do that when the coefficients are time-dependend?

approximate them linearly in time as

$$
g_{l\alpha}(t_i)\simeq g_{l\alpha}^{(0)} + (t_i - t_0) g_{l\alpha}^(1)
$$

Then we get the same form as before:

$$
A_{i\alpha} = Z_i G_\alpha,
$$

where 

$$
Z_i = (X_i, (t_i - t_0)X_i)
$$


This is how TDI actually applies to the phasemeter data channels y. For LISA, we assume six phasemeter data channels, two per each arm $y_{er}$, with $e,r = 1,2,3$.

TDi makes a linear combination of the ys with some delays. This is realized on discretely sampled data with fracitonal delay filtering. 

if we look over the course of a day or so, the TDI coefficients and delays will vary slowly, which we can approximate as say a linear drift.
the sresult can be abtractly summarized  in the last eq.

here Z is a data matrix including the ys with rows that are copies ar various time-shifts up to the required stencil width, and some rows that are also multiplied by the relative sample time. For TDI, G encodes the linear combination of those columns needed to get the TDI channels.


### What about PCI?
For aPCI now forget about the TDI version of $G$ and instead use principal component analysis (PCA) to find a new one which provides linear combinations of the data ordered by their variance. We will be interested in the lowest variance channels, which must lack the laser noise.


To find the coefficient matrix $G$ in practice, we decompose $Z$ by using PCA:

$$
Z = USV\dagger
$$

The components of $S$ are going to be ordered by variance. 

In our case, when applying it to LISA data:
- the high variance components are going to be dominated by laser noise
- the low variance components are what is left-over, including the GW signal

with PCA the data matrix $Z$ is organized into principale components sorted by variances as they appear in the diagonalize matrix $S$, with the help of unitary matrices $U$ and $V$. Our aPCI analog of 

$$
A = Z G
$$

is 

$$
E = Z V_{(q)}
$$

where we truncate the $V$ matrix to keep only a few of the lowest variance columns (marked by $q$), which are our new channel projectors.

We need $q$ to be large enough to get the gravitational wave content, but small enough to avoid significant laser noise.


### Open questions
- how many $q$ components to keep?
- what is the optimal length of a data stretch to use in the channel identification? (right now, hours to days)
- how many time orders to include in the data matrix? (more for a longer stretch, in principle)
- while PCA orders the components by variance, LISA's red noise data can't be naturally treated as zero-mean. An alternative is to select the components with the smalles RMS values. 