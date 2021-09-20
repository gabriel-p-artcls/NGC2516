
# NGC2516

Analysis of the total mass, binary fraction, binary mass ratio, and IMF of the NGC2516 cluster.

<!-- MarkdownTOC levels="1,2,3" autolink="true" style="ordered" -->

1. [TODO](#todo)
    1. [Relevant articles](#relevant-articles)
1. [Analyzed region](#analyzed-region)
1. [Structural analysis](#structural-analysis)
1. [pyUPMASK](#pyupmask)
1. [Membership](#membership)
1. [Fundamental parameters with ASteCA](#fundamental-parameters-with-asteca)
1. [Coarse identification of single/binary systems](#coarse-identification-of-singlebinary-systems)
1. [Coarse per-star mass estimation](#coarse-per-star-mass-estimation)
    1. [Single systems](#single-systems)
    1. [Binary systems](#binary-systems)
1. [Total mass estimation](#total-mass-estimation)

<!-- /MarkdownTOC -->



## TODO

* Cross-match with [Cantat-Gaudin (2018)](http://dx.doi.org/10.1051/0004-6361/201833476) members. According to Bouma et al. (2021) there are 1106 "candidate cluster members" up to G=18 mag.


Read:

4. About the 'systemic PM': https://www.redalyc.org/pdf/571/57115758041.pdf
5. *Mass-loss rates and the mass evolution of star clusters*, Lamers et al. (2010)
6. *An analytical description of the disruption of star clusters in tidal fields with an application to Galactic open clusters*, Lamers (2005)
7. *Clusters in the solar neighbourhood: how are they destroyed?*, Lamers & Gieles (2006)


### Relevant articles

* [What a local sample of spectroscopic binaries can tell us about the field binary population, Fisher et al. (2005)](https://ui.adsabs.harvard.edu/abs/2005MNRAS.361..495F/abstract)
* [Study of the mass-ratio distribution of spectroscopic binaries - I. A novel algorithm, Shahaf et al. (2017)](https://ui.adsabs.harvard.edu/abs/2017MNRAS.472.4497S/abstract)
* [Study of the mass-ratio distribution of spectroscopic binaries - II, Shahaf & Mazeh (2019)](https://ui.adsabs.harvard.edu/abs/2019MNRAS.487.3356S/abstract)
* [Stellar, brown dwarf and multiple star properties from a radiation hydrodynamical simulation of star cluster formation, Bate (2012)](https://ui.adsabs.harvard.edu/abs/2012MNRAS.419.3115B/abstract)
* [A Survey of Stellar Families: Multiplicity of Solar-type Stars, Raghavan et al. (2010)](https://ui.adsabs.harvard.edu/abs/2010ApJS..190....1R/abstract)
* [Modeling Unresolved Binaries of Open Clusters in the Color-Magnitude Diagram. I. Method and Application of NGC 3532, Li et al. (2020)](https://ui.adsabs.harvard.edu/abs/2020ApJ...901...49L/abstract)

A lot more can be found in [NGC 2516 , the SIMBAD biblio](http://simbad.u-strasbg.fr/simbad/sim-id-refs?Ident=NGC%20%202516&Name=NGC%20%202516)



## Analyzed region

```
Source                              RA           DE
---------------------------------------------------
Aladin                        119.51667   -60.75333
Cantat-Gaudin et al. (2020)   119.527       -60.800
```

We downloaded a 4x4 deg^2 region from the Gaia EDR3 survey using Vizier, with a magnitude cut on `G=18` mag. This region contains 114512 stars.



## Structural analysis

The field contamination is most severe below G=15 mag. I performed a structural analysis with a maximum magnitude cut on G=15 mag and with no magnitude cut.

For G<15 mag (and automatic `Nmemb` estimation) the rc and rt are estimated to be ~(26.6, 85) arcmin, respectively, with approximately 330 members within the region.
For the analysis with no magnitude cut the field density is much more difficult to estimate, which affects the estimated number of members. These two parameters affect in turn the King profile fit since they are fixed in the Bayesian inference process.

The automatic field density estimation for the frame with no magnitude cut results in 7308 st/deg^2. Fixing `Nmemb` to different values results in:

```
Nmemb   rc     rt    cmmts on fit
---------------------------------
500     10     65    great in central RDP, not so good in 15'-25'
800     13     96    excellent fit up to 40'
1000    14    112    slightly above in 30'-40' but very reasonable
1200    16    123    above the RDP 10'-50'
1500    20    129    semi-reasonable central part, bad from 10'-50' 
1700    23    130    bad fit
2000    29    133    bad fit
```

The auto value seems a bit high. Lowering the field density to 7000 st/deg^2 and the `Nmemb` to different values results in:

```
Nmemb   rc     rt    cmmts on fit
---------------------------------
500     10     71    great in central RDP, bad in 10'-30'
800     13    100    good in central RDP, better in 10'-30'
1000    15    114    very reasonable overall
1200    17    123    not so good in central RDP, much better in 10'-30'
1500    21    129    worse central RDP than before, still reasonable
1700    25    131    very similar to the run with 1500
2000    30    133    the fit looks clearly inferior
```


## pyUPMASK

We use the `pyUPMASK` code with the `GMM` clustering method using `Plx pmRA pmDE` as input data.


## Membership

A first step in determining the number of members is a manual analysis of the cluster's data. Isolating the overdensity in proper motions and parallax is very easy due to the proximity of this cluster. The number of (member) stars selected with this coarse approach inside a radius of 1.5 deg is approximately 1300.

A second step involves a structural analysis that includes the field density and radial density profile. This is hindered due to the already mentioned heavy and non-uniform contamination in the field.








## Fundamental parameters with ASteCA

After removing 15 stars with incomplete photometry or large errors, the final number of members within the defined cluster region is 1553 as seen in **Fig 3**.

![Alt text](figs/NGC2516_A2.png?raw=true)

<p align="center"><b>Fig 3</b>: selected region for members of NGC2516.</p>

The above selection of members is processed with `ASteCA` v0.3.1 using Gaia EDR3 isochrones. The ranges for the fundamental parameters used were:

```
RZ   0.0001  0.0291
RA     7.5      9.5
RE       0        1
RV     3.1
RD       8      8.5
RM     100    10000
RB       0        1
```

with a minimum  binary mass fraction of `0.7`. The result is shown in **Fig 4** and summarized in the table below.

![Alt text](figs/NGC2516_D3.png?raw=true)

<p align="center"><b>Fig 4</b>: ASteCA results for NGC2516.</p>

```
--------------------------------------------------------------------------
z_mean      z_MAP   z_median     z_mode     z_16th     z_84th      z_std 
0182004  0.0176466  0.0181344  0.0178594  0.0174533  0.0189221 0.000689662

a_mean      a_MAP   a_median     a_mode     a_16th     a_84th      a_std
7.92354     7.9185    7.91928    7.91973    7.91847    7.92642  0.0114528

E_mean      E_MAP   E_median     E_mode     E_16th     E_84th      E_std
0.189961    0.18436   0.188103   0.185361   0.181808   0.202378  0.0115996

d_mean      d_MAP   d_median     d_mode     d_16th     d_84th      d_std
8.11194    8.11616    8.11477    8.11817    8.09084    8.13043   0.020638

M_mean      M_MAP   M_median     M_mode     M_16th     M_84th      M_std
4035.92    4087.11    4033.78    4029.57    3842.57     4232.4    206.798

b_mean      b_MAP   b_median     b_mode     b_16th     b_84th      b_std
0.213004   0.259339    0.22155   0.254616    0.16449   0.261205  0.0435521
--------------------------------------------------------------------------
```




## Coarse identification of single/binary systems





Using the isochrone from the fit found with ASteCA we separate between single and binary systems using the system-isochrone distance. We manually set a maximum distance threshold from the star to the isochrone in the CMD. Stars closer than this value to the isochrone are classified as single systems, and the remaining stars are classified as binaries. An example is shown in **Fig 5**. To the right, we show the distribution of masses for the stars classified as single systems.

![Alt text](figs/NGC2516_CMD.png?raw=true)

<p align="center"><b>Fig 5</b>: example of the splitting between single/binary systems.</p>

This distance threshold is the most important parameter in this analysis. The larger this distance threshold the smaller the estimated binary fraction, and vice-versa.




## Coarse per-star mass estimation

Once all the stars are classified as single/binaries, we use the theoretical isochrone to estimate the systems' masses.

### Single systems

For single systems the process is simply to find the closest point in the isochrone and assign its mass to the observed star.


### Binary systems

For binary systems we need to identify the combination of primary-secondary stars that are most likely to shift the position of the primary star from the isochrone to the current observed location in the CMD. The process is as follows:

1. Given the `(mag, col)` position in the CMD for the observed binary system
2. Propose a range for the magnitude of the primary star. This range is `[mag, mag + 1]`, where `mag` is the observed magnitude.
3. Propose a range for the secondary mass fraction `q=[0, 1]`
4. Randomly generate a magnitude for the primary star from the given range (`mag1`), and a secondary mass fraction (`q`)
5. Find the primary mass `m1` associated to `mag1`
6. Find the secondary mass as `m2=q*m1`
7. Find all the magnitudes for both these masses
8. Combine the magnitudes to generate a synthetic binary system
9. Estimate the photometric distance between this synthetic binary and the observed binary under analysis
10. Minimize this distance repeating 4-9

In **Fig 6** we show this process for an observed binary system (green circle), and the estimated two stars (black circles) that produce the closest synthetic binary system (red square.)

![Alt text](figs/binar_masses.png?raw=true)

<p align="center"><b>Fig 6</b>: example of the binary masses estimation for and observed binary system.</p>


## Total mass estimation

![Alt text](figs/hurley98.png?raw=true)

From Hurley & Tout (1998)


![Alt text](figs/fisher05.png?raw=true)

From Sollima et al. (2010)
