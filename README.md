# Heliolinc3D C++ Implementation #

## Quick-start for developers

To bootstrap your environment and build the Python bindings, do the following:

```
mamba create -c conda-forge -n heliolinx-dev $(cat conda-requirements.txt)
mamba activate heliolinx-dev
pip install -e .
```

We assume you have a C++14 (or higher) compliant C++ compiler installed
(that can be invoked with `c++`)

To build the binary packages for pypi, do the following on a machine where
you have access to docker:

```
docker run --rm -v "$PWD":/io quay.io/pypa/manylinux2014_x86_64 /io/tools/build_wheels.sh
twine upload dist/*
```

## Why download heliolinc? ##

The motivation for this implementation of the heliolinc algorithm is to enable asteroid discovery using the LSST survey strategy of taking just two images of each field per night, rather than the usual practice of taking four images per field per night. This simple change in survey strategy requires a paradigm shift in asteroid detection. Surveys taking the usual four images per night can identify candidate asteroid discoveries based on just a single night's data. With only two images per night, the attempt to do this would result in overwhelming numbers of false positives. Identifying a reasonable discovery candidate requires **linking** multiple detection pairs across multiple nights of data. The new C++ implementation of heliolinc has proven capability to link simulated **and real** asteroid detections across multiple nights, fast enough that ingesting detection catalogs spanning two weeks of output for major surveys -- including detections from more than one observing site -- is computationally tractable. 

This does not mean that heliolinc has rendered obsolete the typical survey strategy of taking four images per field per night. This strategy continues to have significant advantages over a two image per night strategy, even with heliolinc. Relative to LSST, the currently-operating four-image surveys have much better sensitivity to small asteroids that have brief, close encounters with Earth. These include small asteroids passing within the 0.01 AU of Earth, and very small asteroids on their 'final plunge' toward impact with the Earth -- several of which have been discovered by currently operating surveys. The high-performance implementation of the heliolinc algorithm that we present here should **not** be interpreted as a reason for ongoing asteroid surveys to abandon their highly successful four-image strategy.

Instead, heliolinc has something to offer to surveys using the conventional four-image strategy -- even though it was implemented to enable a completely different strategy. What heliolinc can offer to a four-image survey is the likelihood of successful linking and discovery of objects detected fewer than four times on a given night. This can happen for many reasons: the survey might not acquire four images of a particular field because of weather or other contingencies, or the asteroid might not be successfully detected on all four images because of varying image quality, superposition on a star, rotational brightness variations, or a host of other reasons. Such 'missing' detections are most likely for the faintest objects, hence heliolinc has the potential to extend a survey's sensitivity toward fainter magnitudes.

Additionally, heliolinc has the potential to enable remarkable asteroid discoveries from data sets not acquired with the objective of finding asteroids at all.

## To try it out, obtain the following files from GitHub: ##

### Source code directory (heliolinx/src): ###

##### Tracklet creation: #####
```
make_tracklets.cpp
```

##### Tracklet linking (heliolinc proper): #####
```
heliolinc.cpp
heliolinc_omp.cpp
```

##### Post-processing: #####
```
link_purify.cpp
link_planarity.cpp
```

##### Library for solar system dynamics and geometry, etc: #####
```
solarsyst_dyn_geo01.cpp
solarsyst_dyn_geo01.h
```

##### Useful auxiliary codes: #####
```
parse_clust2det_new.cpp
parse_clust2det_MPC80.cpp
modsplit_hlfile.cpp
```

#### Makefile: ####
```
Makefile
```

### General-use input files (heliolinx-aux/tests): ###

```
Earth1day2020s_02a.csv
ObsCodesNew.txt
colformat_LSST_02.txt
heliohyp_rmb00a.txt
```

### Test data files (heliolinx-aux/tests): ###

```
test_TenObjects01a.csv
```

## Compile the heliolinc suite ##

#### In the directory where you have checked out heliolinx: ####

```
cd src
make
make install
```

## What the programs do: ##

**make_tracklets:** Perform image-based pairing to create input pair/tracklet files for heliolinc

**heliolinc:** Link together pairs/tracklets produced by make_tracklets into candidate asteroid discoveries, using (by default) a k-d tree range query for linking. The k-d tree linking is significant because previous versions used the DBSCAN algorithm, which is more sophisticated but actually inappropriate for the specific use case of heliolinc. Significant performance gains were realized by switching from DBSCAN to the simpler k-d tree range query.

**heliolinc_omp:** Multi-threaded version of heliolinc.

**link_purify:** Post-process linkages produced by heliolinc. This includes Keplerian orbit-fitting of every linkage; iterative (one at a time) rejection of astrometric outliers until the RMS astrometric residual falls below a threshold; identification of sets of mutually overlapping linkages; selection of the best linkage within each mutually overlapping set; and rejection of all the other overlapping linkages. The output from `link_purify` is a final set of linkages guaranteed to be non-overlapping (i.e., composed of detections not shared by any other linkage) and that have been successfully fit by orbits with RMS astrometric residual below the specified threshold.

**link_planarity:** A version of `link_purify` that uses pre-screening based on lack of coplanarity of the inferred 3-D positions to reject some outliers prior to the (much more computationally intensive) Keplerian orbit fitting. With a well-chosen coplanarity criterion, `link_planarity` achieves results almost as good as those of `link_purify` with runtimes a factor of a few shorter.

## Testing your installation: ##

### Testing make_tracklets ###

If you've downloaded the recommended test data files, you can test your installation immediately. In a test directory of your choice, accumulate the following files from heliolinx-aux/tests:

```
Earth1day2020s_02a.csv
ObsCodesNew.txt
colformat_LSST_02.txt
heliohyp_rmb00a.txt
test_TenObjects01a.csv
```

Make sure that the heliolinx/bin directory (where the heliolinc executables produced by make install have been copied) is in your path: e.g. put it in your .bashrc file or type:

```
export PATH=$PATH:/home/user/git_heliolinc/heliolinx/bin

```

You can then run make_tracklets, with the following minimal invocation:

```
make_tracklets -dets test_TenObjects01a.csv -earth Earth1day2020s_02a.csv \
 -obscode ObsCodesNew.txt -colformat colformat_LSST_02.txt
```

The final output written to the terminal should be something like:

```
Output image catalog outimfile01.txt, with 123 lines, has been written
Writing paired detection file pairdetfile01.csv with 109 lines
Writing tracklet file trackletfile01.csv with 42 lines
Writing trk2det file trk2detfile01.csv with 109 lines
```

If the minimalist invocation failed, check that you have the input files used in the example:

```
Earth1day2020s_02a.csv
ObsCodesNew.txt
colformat_LSST_02.txt
test_TenObjects01a.csv
```

##### Brief description of the output files: #####

`outimfile01.txt` is the **output image catalog**. It exists to feed information to heliolinc about the time (MJD), boresight right ascension and declination (RA and Dec), and observer position and velocity (in heliocentric Cartesian coordinates) at the midpoint of each image. 

`pairdetfile01.csv` is the **paired detection file**. It is a simplified, reformatted version of the input catalog containing only the detections that make_tracklets linked into a pair or longer tracklet. The columns are limited to information required by heliolinc or likely to be useful in post-processing -- regardless of what additional columns were present in the input file.

`trackletfile01.csv` is the **tracklet file**. A tracklet is defined as a set of observations on a single night that (could) correspond to a single object. Heliolinc (and therefore the make_tracklets output files that feed it) represents each tracklet as a pair of detections. A pair of detections can be converted into a celestial position and an angular velocity -- which, in heliolinc, are combined with a hypothetical heliocentric distance and radial velocity to specify a full 3-D position and velocity. The tracklet file encodes the time by providing an index to the image catalog, but gives the RA, Dec coordinates in full. If the tracklet has more than two points (i.e., is not just a pair), the image indices are for representative images near the start and end of the tracklet, and the RA, Dec positions are from a Great Circle fit evaluated at the times corresponding to the representative images. 

`trk2detfile01.csv` is the **tracklet-to-detection file**. This is a simple, two-column file, with integers in each column. The first column is an index to the tracklet file, starting from 0. The second is an index to the paired detection file. Hence, the trk2det file enables heliolinc and the downstream, post-processing codes to connect the tracklets specified in the tracklet file to specific individual detections listed in the paired detection file.


##### Specifying output file names: #####

Rather than letting ```make_tracklets``` generate default output names, you can specify them as follows:

```
make_tracklets -dets test_TenObjects01a.csv -outimgs outim_TenObjects01a_01.txt \
-pairdets pairdets_TenObjects01a_01.csv  -tracklets tracklets_TenObjects01a_01.csv \
-trk2det trk2det_TenObjects01a_01.csv -earth Earth1day2020s_02a.csv -obscode ObsCodesNew.txt \
-colformat colformat_LSST_02.txt
```


### Testing heliolinc ###

As a first test of heliolinc, try this invocation:

```
heliolinc -imgs outim_TenObjects01a_01.txt -pairdets pairdets_TenObjects01a_01.csv \
-tracklets tracklets_TenObjects01a_01.csv -trk2det trk2det_TenObjects01a_01.csv -obspos Earth1day2020s_02a.csv \
-heliodist heliohyp_rmb00a.txt
```

After printing a bunch of status output, it should end by printing the following:

```
De-duplicating output set of 26 candidate linkages
Final de-duplicated set contains 15 linkages
Writing 15 lines to output cluster-summary file sumfile_test.csv
Writing 152 lines to output clust2det file clust2detfile_test.csv
Automatically calculated reference MJD was 61107.73
```

As with make_tracklets, you may want the specify the output file names. It turns out there is another parameter, the clustering radius (described in more detail below), for which the default value (1e5 km) is smaller than optimal for this data set. Hence, you can try an additional test, explicitly specifying not only the output file names, but also selecting a clustering radius of 2e5 km:


```
heliolinc -imgs outim_TenObjects01a_01.txt -pairdets pairdets_TenObjects01a_01.csv \
-tracklets tracklets_TenObjects01a_01.csv -trk2det trk2det_TenObjects01a_01.csv -obspos Earth1day2020s_02a.csv \
-heliodist heliohyp_rmb00a.txt -clustrad 2e5 -outsum sum_TenObjects01a_01.csv \
-clust2det clust2det_TenObjects01a_01.csv 
```
The final lines of output should be as follows:

```
De-duplicating output set of 51 candidate linkages
Final de-duplicated set contains 27 linkages
Writing 27 lines to output cluster-summary file sum_TenObjects01a_01.csv
Writing 288 lines to output clust2det file clust2det_TenObjects01a_01.csv
Automatically calculated reference MJD was 61107.73
```

The increase in clustering radius has resulted in more linkages being identified (27 vs. only 15 before), and the output files now have customized names, rather than the default sumfile_test.csv and clust2detfile_test.csv.

##### Brief description of the output files: #####

```sum_TenObjects01a_01.csv``` is the **cluster summary file**, which presents a summary of the candidate linkages ('clusters') identified by ```heliolinc```, one cluster per line. It includes information such as the total number of unique points (that is, individual detections) in a linkage; the spread (RMS deviation from the mean) of the corresponding cluster in position and velocity space; the total time spanned, the number of distinct observing nights represented, and other characteristics to be described later. Though linkages/clusters are much more complex than tracklets, the cluster summary file output by heliolinc is closely analogous to the tracklet file output by make_tracklets -- albeit with many more columns, which are required to adequately represent the essential characteristics of the more complex objects.

```clust2det_TenObjects01a_01.csv``` is the **cluster-to-detection file**, which is exactly analogous to the tracklet-to-detection file output by ```make_tracklets```.

##### Clusters vs. Linkages #####

Herein, we use the terms 'cluster' and 'linkage' interchangeably because they represent two different ways of looking at the same object. A 'linkage' is the end product of the heliolinc analysis: a set of detections, spanning multiple nights, that likely all correspond to the same asteroid observed at different times.

A 'cluster' describes how such a linkage appears when identified by the heliolinc algorithm. The algorithm receives detections already organized into tracklets by make_tracklets, and then, by means of a hypothesis about the object's heliocentric distance (and its first two time derivatives), converts each tracklet into a 3-D position and velocity (referred to as 'state vectors'), which are accurate if the hypothesis is. State vectors fully specify an object's orbit, and can be integrated forward or backward in time. The heliolinc algorithm proceeds by integrating the state vectors corresponding to each input tracklet to a constant reference time. Tracklets corresponding to the same object, when transformed into state vectors using an accurate hypothesis for the heliocentric distance and then integrated to the reference time, will all have nearly the same 3-D position and velocity. Hence, after integration to the reference time, a correct linkage appears as a set state vectors that are tightly clustered in both position and velocity space. It is by identifying such clusters that heliolinc forms its output linkages.

To summarize, the language of 'clusters' is relevant when discussing the 'clustering radius' and related parameters in heliolinc -- but to the end user, the desired entity is a 'linkage' of disparate detections that have been identified as likely corresponding to the same object. Hence, both terms are usefully descriptive and we employ either depending on context.


### Testing link_purify ###

In the test immediately above, heliolinc reported writing 27 lines to the summary file, meaning it identified 27 distinct candidate linkages. However, our test data set actually includes simulated observations of only ten distinct objects. There are several reasons why there are more candidates than real objects:

* Multiple non-identical subsets of the detections for a given asteroid might lead to multiple distinct linkages.
* Spurious linkages consisting of detections from more than one asteroid can exist.

In the heliolinc suite, the programs `link_purify` and `link_planarity` exist to deal with these cases and cull the linkages down to an optimized, non-redundant set. They do this by identifying sets of mutually overlapping candidate linkages, choosing the best one, and rejecting all the inferior linkages that overlap (share detections) with it. The criteria used to define the *best* linkages will be described in more detail below. We will describe `link_purify` first, because `link_plarity` is designed to achieve exactly the same result more quickly, at the cost of greater complexity and (very slightly) increased chance of rejecting a good linkage.

To test your installation of `link_purify`, you must first create a text file listing the output from `heliolinc` that will be analyzed. The reason `link_purify` requires a list of files is that it can be used to analyze output from multiple executions of `heliolinc` at one time. Each line of this text file gives the names of both the output files from a single execution of heliolinc: first the summary file (e.g., `sum_TenObjects01a_01.csv`) and then the clust2det file (e.g., `clust2det_TenObjects01a_01.csv`). For this test, your link file list will have only one line since you are analyzing results from just one execution of `heliolinc`. You can construct the list file by hand in a text editor such as emacs, or automatically with a command such as:

```
echo "sum_TenObjects01a_01.csv clust2det_TenObjects01a_01.csv" > TenObjects01a_lflist
```

Its contents should be simply the single line:

> `sum_TenObjects01a_01.csv clust2det_TenObjects01a_01.csv`

You can then test `link_purify` with the invocation:

```
link_purify -imgs outim_TenObjects01a_01.txt -pairdets pairdets_TenObjects01a_01.csv -lflist TenObjects01a_lflist
```

The final status output should read:

```
Writing 10 lines to output cluster-summary file LPsumfile_test.csv
Writing 103 lines to output clust2det file LPclust2detfile_test.csv
```

Here, the writing of 10 lines to the output cluster summary file indicates that exactly ten distinct linkages have been found: the exact number of simulated asteroids included in the test data. The output files `LPsumfile_test.csv` and `LPclust2detfile_test.csv` are formatted identically to the summary and clust2det files, respectively, that were originally outputted by `heliolinc`. This allows `link_purify` (and `link_planarity`) to be run recursively -- that is, its output files can be cycled back as inputs to a new execution. Reasons why you might want to do this will be discussed below.

### Testing link_planarity ###

You can test link_planarity with the same form of invocation as link_purify, on the same input files:

```
link_planarity -imgs outim_TenObjects01a_01.txt -pairdets pairdets_TenObjects01a_01.csv \
-lflist TenObjects01a_lflist
```

The final status output should read:

```
Writing 8 lines to output cluster-summary file LPLsumfile_test.csv
Writing 80 lines to output clust2det file LPLclust2detfile_test.csv
```

This is not as good a result as with `link_purify`: just 8 out of the 10 asteroids have apparently been identified. As with `heliolinc`, the reason for this is that the default value of one of the search parameters turns out not to be optimal for this data set. This parameter is the the maximum RMS distance out-of-plane, or `oop`, which is the key parameter that distinguishes `link_planarity` from `link_purify`.

When `link_purify` reads in a candidate linkage from `heliolinc`, it immediately performs a Keplerian orbit fit. It then rejects bad points from the linkage based on their large astrometric residuals from the best-fit Keplerian orbit.

By contrast, `link_planarity` first reproduces the `heliolinc`-style calculation of initial state vectors: that is, from input RA,Dec coordinates and a hypothesis about the heliocentric distance, it infers the object's 3-D position at the moment of observation. If all the observations in a linkage really do correspond to the same object, and the heliocentric hypothesis is accurate, all the inferred state vectors should be coplanar: they should all lie in the plane of the object's orbit around the sun. Hence, before doing any orbit fitting, `link_planarity` calculates the mean orbit plane of the inferred 3-D positions, and the RMS out-of-plane distance. If this RMS `oop` distance is larger than the threshold, `link_planarity` rejects the point farthest from the mean, and re-evaluates the mean plane and the RMS `oop` distance, iteratively rejecting points until the RMS `oop` has fallen below the threshold. Only then does it proceed to orbit-fitting.

The motivation for this is that the coplanarity analysis is computationally much faster than a Keplerian orbit fit, and therefore by pre-screening bad points by means of their lack of coplanarity, `link_planarity` does fewer computationally expensive orbit fits and, for large data sets, typically runs several times faster than `link_purify`. The disadvantage is that imperfect planarity is not as definitive an indication of a bad point as a large residual from an orbit fit is. For example, even in the case of a correct linkage, a substantially imperfect hypothesis about the heliocentric distance can lead to large deviations from coplanarity of the inferred 3-D positions -- but a Keplerian orbit fit will usually still converge with small astrometric residuals.

Hence, `link_planarity` is at greater risk of discarding good points than is `link_purify`. Rejecting points can lead to a whole linkage being rejected if it falls below the minimum criteria for a valid discovery (at least two detections per night on three distinct nights). This is why the example run above found only 8 linkages rather than the total 10 that were in the data. To prevent these bad outcomes, it is necessary to set the maximum permitted RMS `oop` distance to a value permissive enough (that is, large enough) that the coplanarity test will almost never reject a point that the Keplerian orbit fit would pass. For the simulated LSST test data used here, we have found that 10000 km is a good choice, while the default used in `link_planarity` is 1000 km. Setting this value explicitly, and specifying the output file names as well, we can re-run `link_planarity` as follows:


```
link_planarity -imgs outim_TenObjects01a_01.txt -pairdets pairdets_TenObjects01a_01.csv \
-lflist TenObjects01a_lflist -oop 10000.0 -outsum LPLsum_TenObjects01a_01.csv \
-clust2det LPLclust2det_TenObjects01a_01.csv
```

Now, the final status output should read:

```
Writing 10 lines to output cluster-summary file LPLsum_TenObjects01a_01.csv
Writing 103 lines to output clust2det file LPLclust2det_TenObjects01a_01.csv
```

With the relaxed out-of-planarity threshold `oop`, `link_planarity` has found all of the valid linkages. In this example with a tiny data set, the faster runtime enabled by `link_planarity` is not evident, but if millions (or even just tens of thousands) of candidate linkages are being analyzed, the speedup can be very valuable. For a well-chosen value of the `oop` parameter, the losses from using `link_planarity` rather than `link_purify` for the inital screening are negligible.

##### Brief description of the output files: #####

The cluster summary files output by the post processing codes `link_purify` and `link_planarity` have the same format as those output by `heliolinc`. However, some columns are populated by the post processing codes that are left empty (that is, filled with dummy placeholder values) by `heliolinc`. These are columns holding orbit-fitting information, such as Keplerian semimajor axis and eccentricity -- but more importantly, the full state vectors (X, Y, Z) and (VX, VY, VZ) for the best-fit orbit at a customized epoch MJD (also provided) near the center of the time spanned by the observations in the linkage. Since the post processing codes typically reject (as they should) a large number of spurious and overlapping/redundant linkages output by `heliolinc`, the cluster summary files from the post processing are generally much shorter than those output directly by heliolinc.

The cluster-to-detection files output by `link_purify` and `link_planarity` have exactly the same meaning and format as those output by `heliolinc`. Like the cluster summary files to which they provide indices, the cluster-to-detection files ouput by the post processing codes are generally much shorter than those originally written out by `heliolinc`.

Since the files output by the post processing codes have exactly the same format as the files they take as input, it is possible to run `link_purify` or `link_planarity` recursively, using as input the files that either program wrote out on a previous execution. The main reason you might want to do this is to combine linkages from multiple separate executions of post processing codes (and potentially of heliolinc also) into a single optimized and de-duplicated master set. This will be discussed more below.

This completes the test run of the regular heliolinc suite of programs. If all the tests have been successful, you probably want to run the heliolinc suite on your own data. In order to accomplish that, all you need to do is convert your data into CSV (comma separated values) format and supply a column formatting file, which, as described below, tells `make_tracklets` which columns in your CSV file hold the information it requires.

## Understanding the File Formats ##

### Original Input Data ###

To input new astronomical data into the heliolinc suite, you just need to create a single file: the input data file for `make_tracklets`. The format required is a CSV file with a one-line header followed by one line for each detection, with columns containing at least the following quantities:

* A string identifier for each detection (max length 19 characters).
* The time the detection was made (e.g. mid-exposure) in Modified Julian Days (MJD)
* The right ascension (RA) of the detection, in decimal degrees.
* The declination (Dec) of the detection, in decimal degrees.
* The magnitude of the detection.
* The photometric band in which the observation was made (a single character such as u, g, r, i, etc).
* The 3-character observatory code of the observing site (e.g. 695 for Kitt Peak or I11 for Gemini South)

Any number of additional columns is allowed, and the required columns can be in any order. The columns that hold the required data are communicated to `make_tracklets` through a column format file, introduced on the command line by the `-colformat` keyword, which simply tells it which column has each of the seven required quantities. For example, the full contents of the column format file `colformat_LSST_02.txt` used in our example invocation above are:

```
IDCOL 2
MJDCOL 3
RACOL 6
DECCOL 7
MAGCOL 11
BANDCOL 9
OBSCODECOL 16
```
That's all there is to it. To run `make_tracklets` on your own data, create a CSV file with any columns you want (as long as it includes the required seven) and then create a column format file that tells `make_tracklets` where to find the seven things it needs. Note that it starts counting columns from 1, not from 0.

**Running with missing columns:** in principle, `make_tracklets` and the rest of the heliolinc suite can run with only three inputs per observation: MJD, RA, and Dec. There are many reasons **not** to do this, but there may be some cases where it's helpful. In order to make it possible, `make_tracklets` has a `-forcerun` keyword that will cause it to tolerate input files missing the string ID (which then defaults to "null"), the magnitude (defaults to 99.999), the photometric band (defaults to V), and/or the observatory code (defaults to 500 -- i.e., the geocenter). Including `-forcerun` (or `-f`) as an argument to `make_tracklets` will cause it to run (but print warnings) if any or all of these are missing. In that case, they should also be omitted from the column formatting file -- hence, this file could have as few as three lines, e.g.:

```
MJDCOL 3
RACOL 6
DECCOL 7
```
Importantly, the four omitted quantities are not created equal. The magnitude, photometric band, and string ID are not explicitly used in calculations by programs in the `heliolinc` suite, so omitting them has no bad effects other than making them unavailable for downstream analysis. By contrast, the observatory code **is** used in calculation. In particular, `make_tracklets` uses it to generate Cartesian coordinates for the observer, which are passed on to `heliolinc` itself. Bad observer positions result not only in bad inferred positions for the objects (which may be close enough to work OK), but also in bad inferred velocities that are probably **not** close enough to work OK. These bad velocities can keep `heliolinc` from identifying perfectly valid linkages of real asteroids. Hence, running `make_tracklets` without an observatory code column is **generally inadvisable**. Cases where it might work include those where the data are relatively sparse, so a very large (~1e6 km) clustering radius can be used in `heliolinc` without much confusion; and those where the objects being sought are very distant (>20 AU), so the difference between geocentric and topocentric coordinates is less significant in angular terms.

**Supplying additional columns:** Besides containing an unlimited number of columns not used or retained by the `heliolinc` suite program, your input data file may optionally supply the following, which are retailed by `heliolinc` and included in the paired detection file:

##### Optional columns in the input file #####

Name                    | Description
:---                    | :---
`trail_len`             | trail length in arcsec, if source is trailed
`trail_PA`              | trail position angle, degrees east from celestial north. If the source is not trailed, defaults to 90.0 degrees.
`sigmag`                | magnitude uncertainty
`sig_across`            | cross-track astrometric uncertainty (Dec uncertainty if not trailed)
`sig_along`             | along-track astrometric uncertainty (RA uncertainty if not trailed)
`known_obj`             | Integer quantifying the likelihood that this is a known object, 0 to 999
`det_qual`              | Integer quantifying the likelihood that this detection is real, 0 to 999

While the `origindex` in the paired detection file (see below) enables mapping of final linkages output by the heliolinc suite all the way back to line numbers in the original input file, explicitly including the above optional data makes mapping those specific quantities much easier. This is because `make_tracklets` will write the data directlyt to the paired detection file, which is read and indexed by both `heliolinc` and its post-processing programs. Here are example lines you would add to your column formatting file to read all the optional data from the table above, assuming it to be present in your input file.

```
TRAILPACOL 19
TRAILLENCOL 20
SIGMAGCOL 6
SIGACROSSCOL 12
SIGALONGCOL 15
KNOWNOBJCOL 17
DETQUALCOL 16
```

These quantities will be revisited in the discussion of the output paired detection file below. We briefly note here that, if accurately supplied by the input file, they enable in-depth analysis after the `heliolinc` suite processing is complete. This analysis could include: (1) Evaluating whether all the detections in a linkage likely corresponded to a known object. (2) Identifying a probably-spurious linkage by the fact that all detections had a low detection quality. (3) Bolstering confidence in a linkage corresponding to a fast-moving NEO, by noting that all its individual detections were trailed with mutually consistent trail length and orientation. 

### Input Image File ###

If you have easy access to information about the time (MJD) of mid-exposure and the bore-sight pointing for each image in your input data, you can (and should) supply this information to `make_tracklets` in the form of a 4-column file that has **no header** and is **not a CSV** but has columns separated by spaces. The name of this file is supplied to `make_tracklets` using the command line keyword `-imgs`. The four columns in the file are MJD, RA, Dec, and observatory code. The RA and Dec must be in decimal degrees. The observatory code must be a three-character string, and it and the MJD **must match** the corresponding entries in the input detection data file. Here is an example of the first few lines in an acceptable input image file.

```
61100.29125 180.85436 0.74503 X05
61100.29169 180.85428 0.74508 X05
61103.01003 110.54100 12.43000 X05
61103.01048 110.54100 12.43000 X05
```
If you do not have easy access to this information, don't worry: `make_tracklets` reconstructs it internally from the input detection data. If you want `make_tracklets` to write out its internally-reconstructed image log in the required format, supply the name of an output image file to be written using the command line keyword `-outimgs`. The output file has eight additional columns. The first six record the observer state vectors (X,Y,Z) and (VX,VY,VZ) at the time of the exposure. The last two record the first and last entries (actually, one after the last entry) in the output paired detection catalog (see below) that come from that image. For example:

```
61100.29125 180.85436 0.74503 X05 -139581193.4 49872573.9 -4829.9 -10.4089 -28.5146 0.1586 0 1
61100.29169 180.85428 0.74508 X05 -139581589.0 49871489.9 -4823.9 -10.4076 -28.5145 0.1585 1 2
61103.01003 110.54100 12.43000 X05 -141885389.1 43212107.0 -6841.1 -9.5603 -28.6235 0.0159 2 3
61103.01048 110.54100 12.43000 X05 -141885760.8 43210994.1 -6840.4 -9.5600 -28.6246 0.0163 3 4
```
As regards the last two columns, this is an unrealistic example because our test data set is a tiny subset of a full-scale run, and therefore there is only one paired detection on each of the input images. This is why the indices in the last two columns advance by only one per line. Here is the beginning of a file actually run on real data, with approximately one thousand paired detections per image:

```
59858.000134   8.219640 -42.054112 M22 146089786.4 32164326.5 -6402.1 -7.1231 29.2645 -0.1290 0    994
59858.000232 289.089807 -40.000181 W68 146088372.4 32157591.3 -3174.8 -6.5624 29.1733 -0.0908 994  2239
59858.000683 277.645892 -24.269771 W68 146088117.0 32158726.6 -3178.3 -6.5633 29.1741 -0.0912 2239 3297
59858.000806   9.759413 -52.596351 M22 146089372.4 32166027.1 -6409.6 -7.1248 29.2636 -0.1286 3297 4075
```

Because `make_tracklets` simply ignores input image columns beyond the first four, the output image file written by one invocation of `make_tracklets` can be used as an input image file for future invocations on detections from the same images. A case where you might want to do this is if the output image file was produced by `make_tracklets` running on all available data, and you were re-running `make_tracklets` after aggressively culling the input data down to a much smaller number of detections. Because `make_tracklets` estimates the bore-sight RA and Dec of a given image by simply averaging over all the detections taken at the same time and observatory code, the old values based on a larger number of detections will be more accurate than what can be reconstructed from the smaller, culled set of detections.

### Observer Location Files ###

##### Executive summary: just use the ones we gave you. #####

The remaining files required by `make_tracklets` are the ones used to determine the observer's position in the solar system at the time of each observation. In our example invocation, these are `Earth1day2020s_02a.csv` and `ObsCodesNew.txt`, which are introduced with the command line keywords `-earth` and `-obscode`, respectively. The files used in the example are broadly applicable, and it's unlikely you'll need different files. The `-earth` file specifies the Earth's position relative to the sun as a function of time, while the `-obscode` (observatory code) file allows `make_tracklets` to look up the location of the observatory on the Earth's surface based on the three-character observatory codes specified in the input detection file. 

Cases where you might need different files include observations before December 2019 or after December 2030; or made from a new observatory that doesn't have an existing observing code. You can get a new `-earth` file using the JPL Horizons ephemeris generator (`https://ssd.jpl.nasa.gov/horizons/app.html`). Set Ephemeris Type to `Vector Table`; Target Body to `Earth` (**not** the Earth-Moon barycenter); and Coordinate Center to `Sun (body center)` (try entering '@sol' or '@10' in the search box). Set Time Specification to whatever range of dates you want, at 1-day sampling. Under 'Table Settings', click the box for CSV format and leave everything else set to the defaults. Click `Generate Ephemeris` and then `Download Results`; save the resulting file; and feed it directly into `make_tracklets` using the command line keyword `-earth`.

For a new observatory code, your best option is probably to go to `https://minorplanetcenter.net/iau/lists/ObsCodesF.html`, find the lines for the new observatories you need, and just add them to the existing file `ObsCodes.txt`. Alternatively, you can copy and paste the whole file. Alternatively, `make_tracklets` is able to directly ingest the unformatted html file at `https://minorplanetcenter.net/iau/lists/ObsCodes.html`.

### Output files from make_tracklets ###

The format of the output files is CSV (comma separated values) -- except for the image catalog, which is space-separated and only used internally. Most of the files have single-line headers that should give you some idea of what the columns mean. Common formats have been used when possible. For example, both the tracklet-to-detection file output by `make_tracklets` and the cluster-to-detection files output by `heliolinc` and by the post-processing codes `link_purify` and `link_planarity` have the same format: they are two-column files connecting long-integer indices to two diffent catalogs. Similarly, the cluster summary files output by `heliolinc` and by the post-processing codes all have the same format. The files with unique formats (by necessity) are the paired detection catalog and the tracklet file, both output by `make_tracklets`.

#### Paired detection file ####

The paired detection file ouput by `make_tracklets`, called `pairdets_TenObjects01a_01.csv` in our test example, is essentially a culled and reformatted copy of the input detection catalog (`test_TenObjects01a.csv` in our example). It contains only the detections that were included in a tracklet by `make_tracklets`. By default, the shortest (and typically most numerous) class of tracklets is those consisting of just two points -- i.e., pairs -- hence our description of this as the **paired** detection file.

Besides culling the input down to only those detections that formed tracklets, the paired detection file also includes only those columns from the input data that are relevant for `heliolinc`. Many of these columns, though potentially valuable, are not essential, so it is acceptable to provide an input data catalog that lacks them. In this case, the missing columns must also be omitted from the column formatting file, `colformat_LSST_02.txt` in the example here. The essential columns are a string identifier for each detection (which can have from 1-19 characters); the time (given as a Modfied Julian Day, or MJD); the RA and Dec; the magnitude; the photometric band; and the observatory code. 


The optional output image file introduced with the `-outimgs` command line keyword has already been described above. The primary output files produced by `make_tracklets` are the catalog of paired detections ('paired detection file') and the file specifying the pairs and tracklets ('pair file').

The default name of the output **paired detection file** is `pairdetfile01.csv`. A different name can be specified using the command line keyword `-pairdets`. This file echoes the data in the input detection catalog, but in a standardized form expected by `heliolinc`, and with two important changes. First, the observer's position in Cartesian ecliptic coordinates relative to the center of the sun at the time of each detection is calculated and recorded. Second, only detections that were included in some viable pair or tracklet are included -- hence, the paired detection file does not necessarily have as many entries as the input detection catalog. Here is a description of the columns:

##### Columns in the paired detection file #####

Column &nbsp;&nbsp;| Name                    | Description
:---               | :---                    | :---
1                  | `MJD`                   | Modified Julian Day
2 	           | `RA`                    | right ascension in decimal degrees
3 	           | `Dec`                   | declination in decimal degrees
4                  | `mag`		     | magnitude
5                  | `trail_len`             | trail length in arcsec, if source is trailed
6                  | `trail_PA`              | trail position angle, degrees east from celestial north. If the source is not trailed, defaults to 90.0 degrees.
7 	           | `sigmag`                | magnitude uncertainty
8                  | `sig_across`            | cross-track astrometric uncertainty (Dec uncertainty if not trailed)
9                  | `sig_along`             | along-track astrometric uncertainty (RA uncertainty if not trailed)
10                 | `image`                 | image index to the output image catalog (internally calculated)
11 	           | `idstring`              | string identifier copied from input detection catalog (max 19 characters).
12 	           | `band`                  | photometric band, single character: u, g, r, i, etc
13	           | `obscode`               | Observatory code, three-character string
14                 | `known_obj`             | Integer quantifying the likelihood that this is a known object, 0 to 999
15                 | `det_qual`              | Integer quantifying the likelihood that this detection is real, 0 to 999
16 	           | `origindex`             | Line number of corresponding entry in the original input detection
                   |                         | catalog (internally calculated)


These quantities are not restricted to those required, or even used, by the programs in the baseline suite. We have already seen that `make_tracklets` and heliolinc can run with only MJD, RA, and Dec -- and run without degredation with only these three columns plus the Observatory Code. Instead of a list of strictly required quantities, we have tried to maintain the paired detection file as a carefully optimized set of quantities to which users will likely want convenient access in some important contexts. The emphasis here is on **convenient**: with a little customized index-tracing, access to any column in the freely formatted input file can be obtain through the `origindex` column. In other words, `origindex`, though not explicitly used by any program in `heliolinc` suite, is retained to enable mapping of the final linkages back to lines in the original detection catalog. This is likely to be useful in the (quite probable) case that you had (and wrote into your input detection catalog) a bunch of interesting data besides the columns preserved in the paired detection file. You can then use the preserved line numbers to recover this original, detailed information for every linked detection. Note that `origindex` starts counting from zero on the first data line of your input detection catalog. Since this file is expected to have a one-line header, origindex=0 actually corresponds to the second line of the file: the first line is a header, and the second line is the first **data** line and is assigned origindex=0. Note also that the user does not include `origindex` as a column in the original input file to `make_tracklets`: rather, `make_tracklets` calculates and outputs it regardless of what other columns are presentin the original input file.

Besides the internally calculated `origindex`, the quantities in the paired detection file can be divided into three catagories: Used internally (if available) by all versions of the `heliolinc` suite; used only for processing of trailed sources; not currently used by `heliolinc`, but may be used in future versions and/or very useful for externally evaluating `heliolinc` results. The columns in each category are as follows:

##### Paired detection file: columns always used by `heliolinc` when available #####

Column &nbsp;&nbsp;| Name                    | Description
:---               | :---                    | :---
1                  | `MJD`                   | Modified Julian Day
2 	           | `RA`                    | right ascension in decimal degrees
3 	           | `Dec`                   | declination in decimal degrees
4                  | `mag`		     | magnitude
10                 | `image`                 | image index to the output image catalog (internally calculated)
11 	           | `idstring`              | string identifier copied from input detection catalog (max 19 characters).
12 	           | `band`                  | photometric band, single character: u, g, r, i, etc
13	           | `obscode`               | Observatory code, three-character string

##### Paired detection file: columns used for trail analysis  #####

Column &nbsp;&nbsp;| Name                    | Description
:---               | :---                    | :---
5                  | `trail_len`             | trail length in arcsec, if source is trailed
6                  | `trail_PA`              | trail position angle, degrees east from celestial north. If the

##### Paired detection file: columns not currently used, retained for future improvements or external evaluation #####

Column &nbsp;&nbsp;| Name                    | Description
:---               | :---                    | :---
7 	           | `sigmag`                | magnitude uncertainty
8                  | `sig_across`            | cross-track astrometric uncertainty (Dec uncertainty if not trailed)
9                  | `sig_along`             | along-track astrometric uncertainty (RA uncertainty if not trailed)
14                 | `known_obj`             | Integer quantifying the likelihood that this is a known object, 0 to 999
15                 | `det_qual`              | Integer quantifying the likelihood that this detection is real, 0 to 999
16 	           | `origindex`             | Line number of corresponding entry in the original input detection
                   |                         | catalog (internally calculated)

#### Tracklet file ####

The tracklet file has eight columns. The first six comprise two sets of three, specifying the start and end of the tracklet, respectively. The final two columns, number of points and tracklet ID number, describe the tracklet as a whole.

##### Columns in the tracklet file #####

Column &nbsp;&nbsp;| Name        | Description
:---               | :---        | :---
1                  | `Image1`    | Index in the `make_tracklets` output image catalog of the image
                   |             | with the first representative detection
2 	           | `RA1`       | First representative RA (decimal degrees)
3 	           | `Dec1`      | First representative Dec (decimal degrees)
4                  | `Image2`    | Index in the `make_tracklets` output image catalog of the image
                   |             | with the second representative detection
5 	           | `RA2`       | Second representative RA (decimal degrees)
6 	           | `Dec2`      | Second representative Dec (decimal degrees)
7 	           | `npts`      | Number of unique detections included in the tracklet
8                  | `trk_ID`    | Tracklet index (counts sequentially from zero)

In the case of a two-point tracklet, `heliolinc` requires the time (MJD), RA, and Dec of the first observation and of the second observation -- and additionally the observer's location (in 3D heliocentric Cartesian coordinates) at the instant of each observation. This sounds like six quantities (MJD, RA, Dec, plus observer X, Y, Z) for **each** observation in the tracklet. We get by with only three (saving memory space, read time, etc.) by taking advantage of the fact that the time and the observer's location is the same for all detections on a given image. Hence, we supply just the image index in the tracklet file. Then `heliolinc`, having read the image file output by `make_tracklets` into a relatively short internal array, uses the image index to extract the MJD and observer X, Y, Z corresponding to a given detection almost instantly.

If the tracklet is just a pair (i.e., it contains only two points) the "first representative" RA and Dec are just the RA and Dec of the first point, and similarly for the second point. If the tracklet has more than two points, `make_tracklets` still preserves a simple API with respect to `heliolinc` by representing it as just two "representative" points. These are, for example, the first and last points in a tracklet with three or four points, but the second and fifth point in a six-point tracklet; and the third and eighth point in a ten-point tracklet. The representative RA and Dec for these multi-point tracklets are not simply the RA and Dec of the representative points: using those values would needlessly discard the information provided by the unused points. Instead, `make_tracklets` performs a constant-velocity Great Circle fit to the tracklet, and evaluates this fit at the times corresponding to the two representative points to get the respective RA and Dec values. This fit averages down the noise in the astrometry and hence provides values more precise than those of individual detections would be. For tracklets with five or more points, the representative detections are not the absolute first and last detections in the tracklet, in order to ensure that the Great Circle fit is evaluated in a regime where it is as constrained as possible by the data. The MJD and observer X, Y, Z coordinates for the reference points, being accurately known, are not determined by a fit but are read directly from the image arrays as usual.

Readers very familiar with tracklet-based asteroid discovery may worry that the "representative point" scheme outlined above will produce bad results in the case of an NEO passing very close to the Earth, because such an object's on-sky motion can deviate significantly (arcseconds) from a constant velocity Great Circle even over a short time such as 30-60 minutes. This danger is avoided by means of a maximum Great Circle residual (maxGCR) parameter in `make_tracklets`. Multi-point tracklets with GCR greater than the maximum threshold are broken up into pairs, which will then be individually linked by `heliolinc` using an actual physical orbit rather than a geometrical approximation such as a Great Circle.

#### Tracklet-to-detection file ####

Abbreviated as the trk2det file, this is a simple, two-column file, with integers in each column. The first column is an index to the tracklet file, starting from 0. The second is an index to the paired detection file. Hence, the trk2det file enables heliolinc and the downstream, post-processing codes to connect the tracklets specified in the tracklet file to specific individual detections listed in the paired detection file. The trk2det file preserves the information, otherwise lost in the tracklet file itself, of which detections correspond to multi-point tracklets. For example, the `heliolinc` algorithm reads a ten-point tracklet from the tracklet file exactly as it would a two-point tracklet (i.e., a pair), and performs the identical processing on it. After linking all the tracklets for a given object, however, `heliolinc` and the post-processing codes use the trk2det file to map the tracklet ID numbers in a linkage back to individual detections. This ability is particularly important since multiple tracklets in a linkage may contain some of the same detections, so only using the mapping provided by the trk2det file can `heliolinc` or post-processing codes determine how many **unique** detections are actually included in a given linkage.

### Input Files for heliolinc ###

The most important `heliolinc` input files are the ones generated by `make_tracklets` and described above. It also requires an Earth position file, but we have already discussed this since it is also a required input for `make_tracklets`. The only `heliolinc` input file that has yet to be described is the one containing the hypotheses about heliocentric radial motion, which is introduced by the command line keyword `-heliodist`, and was called `heliohyp_rmb00a.txt` in our example invocation. This file has three columns. It has a one-line header and is **not a csv** but rather a space-delimited file like the image files read and written by `make_tracklets`. The four required columns are as follows:

Column &nbsp;&nbsp;| Name            | Description
:---               | :---            | :---
1                  | `r(AU)`         | Distance from the sun at the reference time
2 	           | `rdot(AU/day)`  | Radial velocity at the reference time
3 	           | `mean_accel`    | Time-derivative of radial velocity, divided by (-GMsun/r^2)

The solar gravity always imposes, on any object in the solar system, a vector acceleration of GMsun/r^2 toward the center of the sun. The `mean_accel` column in the hypothesis file is not the vector acceleration but rather the time-derivative of the heliocentric radial velocity (equivalently, the second time-derivative of the heliocentric distance). We are inclined to call this the **radial acceleration**, and it is generally not equal to the vector acceleration. For example, in a circular orbit the radial velocity and radial acceleration are both always zero, even though the vector acceleration is still GMsun/r^2 toward the heliocenter. This is why we use units of -GMsun/r^2 for the radial acceleration term `mean_accel`: the value is 0 for a circular orbit and 1.0 for an object with zero angular momentum (Keplerian eccentricity exactly 1.0) that is falling directly toward the center of the sun.

Realistic values for actual bound orbits in the solar system are typically between -1.0 and +1.0, with a distribution that is highly dependent on the heliocentric distance and radial velocity. This is why we used integration of actual asteroid orbits to generate heliocentric hypothesis files.

Besides the very short example file `heliohyp_rmb00a.txt` provided for test runs, we have included the following additional heliocentric hypothesis files:

```
radhypo_mb02a.txt
heliohyp_TNO_02a.txt
radhyp_NEOF_01b.txt
radhypo_ie01a.txt
```

Here, `radhypo_mb02a.txt` targets main-belt asteroids with 213 different hypotheses, radhyp_NEOF_01b.txt targets near-Earth asteroids (NEOs) outside the Earth's orbit with 172095 hypotheses, radhypo_ie01a.txt targets NEOs interior to Earth's orbit with 19933 hypotheses, and `heliohyp_TNO_02a.txt` targets distant, Trans-Neptunian objects (TNOs) with 14 hypotheses. Each is the best hypothesis set we have tested in their respective regions.

The short example file provided for test runs, `heliohyp_rmb00a.txt`, probes only 6 hypotheses, but actually finds most of the possible main belt asteroids (though not as many as radhypo_mb02a.txt) in realistic simulations of LSST data. This is illustrative of a common feature in our tests: a relatively small number of hypotheses find a majority of target objects, but finding the last few percent often requires increasing the number of hypotheses by a large factor. The reasons for this likely include the fact that only a small proportion of objects will have their detections distributed in a maximally difficult way (e.g., smallest possible number of detections spread over the longest possible time interval), or will be in unusual orbits not effectively probed by smaller hypothesis sets.

#### Splitting `heliolinc` hypothesis files, and recombining output ####

For analyses involving large number of hypotheses, such as those targeting NEOs, you can break the hypothesis files into smaller pieces for embarrassingly parallel runs. For example, the file radhyp_NEOF_01b.txt could be broken into 35 pieces with 5000 hypotheses each (except the last, which would have only 2095. **Don't forget that each individual hypothesis file needs its own one-line header**. There is also a multithreaded version of `heliolinc`, called `heliolinc_omp`, which will be discussed below.

The NEO files with tens of thousands of hypotheses have been run successfully only on carefully chosen subsets of input data: trailed sources in the case of radhyp_NEOF_01b.txt, and detections at less than 90 degrees solar elongation in the case of radhypo_ie01a.txt. Processing unculled input data, such as is expected for runs targeting the main belt with radhypo_mb02a.txt, might not be computationally tractable with the NEO files. If you want to find NEOs without culling the input data (e.g., if you are targeting NEOs outside Earth's orbit that might not be moving fast enough to make trails on your images), a way to make this computationally tractable is to cut out the the hypotheses closest to Earth, where the sampling is finest. For example, radhyp_NEOF_01b.txt has 172095 hypotheses, but only 53135 correspond to asteroids more than 1.1 AU from the sun, and only 10761 correspond to those more distant than 1.2 AU. Using only the more distant hypotheses can greatly reduce the runtimes.

If you divide a hypothesis file, such as radhyp_NEOF_01b.txt, into parts with extensions _p00, _p01, _p02, etc, and you want to run it on a set of data called `mydata` that you have already processed with `make_tracklets`, you can run `heliolinc` on it as follows:

```
heliolinc -imgs outim_mydata.txt -pairdets pairdets_mydata.csv -tracklets tracklets_mydata.csv \
-trk2det trk2det_mydata.csv -mjd 59816.00 -obspos Earth1day2020s_02a.csv -heliodist radhyp_NEOF_01b_p00.txt \
-clustrad 500000.0 -mingeodist 0.01 -mingeoobs 0.005 -minimpactpar 50000.0 -outsum sum_mydata_NEOF_01b_p00.csv \
-clust2det clust2det_mydata_NEOF_01b_p00.csv

heliolinc -imgs outim_mydata.txt -pairdets pairdets_mydata.csv -tracklets tracklets_mydata.csv \
-trk2det trk2det_mydata.csv -mjd 59816.00 -obspos Earth1day2020s_02a.csv -heliodist radhyp_NEOF_01b_p01.txt \
-clustrad 500000.0 -mingeodist 0.01 -mingeoobs 0.005 -minimpactpar 50000.0 -outsum sum_mydata_NEOF_01b_p01.csv \
-clust2det clust2det_mydata_NEOF_01b_p01.csv

heliolinc -imgs outim_mydata.txt -pairdets pairdets_mydata.csv -tracklets tracklets_mydata.csv \
-trk2det trk2det_mydata.csv -mjd 59816.00 -obspos Earth1day2020s_02a.csv -heliodist radhyp_NEOF_01b_p02.txt \
-clustrad 500000.0 -mingeodist 0.01 -mingeoobs 0.005 -minimpactpar 50000.0 -outsum sum_mydata_NEOF_01b_p02.csv \
-clust2det clust2det_mydata_NEOF_01b_p02.csv

etc.
```
The arguments -mingeodist, -mingeoobs, and -minimpactpar used here can be important for good results on NEOs, and will be explained in detail below.

The output files can all be fed to the post processing code (`link_purify` or `link_planarity`) at one time. For example, if you create a file called clusterlist_mydata_NEOF_01b, with contents as follows:

```
sum_mydata_NEOF_01b_p00.csv clust2det_mydata_NEOF_01b_p00.csv
sum_mydata_NEOF_01b_p01.csv clust2det_mydata_NEOF_01b_p01.csv
sum_mydata_NEOF_01b_p02.csv clust2det_mydata_NEOF_01b_p02.csv

etc.
```
You can then run `link_planarity` as follows:

```
link_planarity -imgs outim_mydata.txt -pairdets pairdets_mydata.csv -lflist clusterlist_mydata_NEOF_01b \
-oop 10000.0 -outsum LPLsum_mydata_NEOF_01b.csv -clust2det LPLclust2det_mydata_NEOF_01b.csv
```
It is because they need to have this capability of combining multiple files that the post processing codes require an "lflist": that is, a file listing the heliolinc output files on which they are to be run. Otherwise we would have made them accept the heliolinc output filenames directly, e.g., with command line keywords -insum and -inclust2det.

If the individual output files from a "split" run of heliolinc are very large, feeding all of them into the post processing code at the same time might be undesirable either because of excessive runtime or memory usage. In this case, a post-processing code such as `link_planarity` can be run individually on each file, and then re-run on the whole list of smaller output files from each of the individual runs. This is made possible by the fact that the post processing codes use input and output files of identical format: therefore they can re-ingest their own output whenever needed. In this case, the split run of heliolinc could be followed by a split run of `link_planarity`:

```
link_planarity -imgs outim_mydata.txt -pairdets pairdets_mydata.csv -lflist clusterlist_mydata_NEOF_01b_p00 \
-oop 10000.0 -outsum LPLsum_mydata_NEOF_01b_p00.csv -clust2det LPLclust2det_mydata_NEOF_01b_p00.csv

link_planarity -imgs outim_mydata.txt -pairdets pairdets_mydata.csv -lflist clusterlist_mydata_NEOF_01b_p01 \
-oop 10000.0 -outsum LPLsum_mydata_NEOF_01b_p01.csv -clust2det LPLclust2det_mydata_NEOF_01b_p01.csv

link_planarity -imgs outim_mydata.txt -pairdets pairdets_mydata.csv -lflist clusterlist_mydata_NEOF_01b_p02 \
-oop 10000.0 -outsum LPLsum_mydata_NEOF_01b_p02.csv -clust2det LPLclust2det_mydata_NEOF_01b_p02.csv

etc.
```

Where the lflist files clusterlist_mydata_NEOF_01b_p00, clusterlist_mydata_NEOF_01b_p01, etc., each have only a single line, giving the output file names from a single execution of `heliolinc`.

The split runs of `link_planarity` in this example could be run in parallel if they were split to reduce runtime; or in series if they had to be split to reduce memory usage. The output files (`LPLsum_mydata_NEOF_01b_p00.csv`, etc.), will be much smaller than the original heliolinc files from which they came. Each individual output file will have only non-overlapping linkages that could be fit with a Keplerian orbit that did not have enormous astrometric residuals. However, there will still be overlap, duplications, etc. across the output files from the different split runs of `link_planarity`. To resolve this duplication and choose the overall best linkage from each mutually overlapping set, you can run an additional round of post processing. To begin with, you would construct a single, additional lflist file called clusterlist_mydata_NEOF_01b, with contents:

```
LPLsum_mydata_NEOF_01b_p00.csv LPLclust2det_mydata_NEOF_01b_p00.csv
LPLsum_mydata_NEOF_01b_p01.csv LPLclust2det_mydata_NEOF_01b_p01.csv
LPLsum_mydata_NEOF_01b_p02.csv LPLclust2det_mydata_NEOF_01b_p02.csv

etc.
```
For combining files that have already had one round of post-processing, we recommend `link_purify` rather than `link_planarity`, because the significant speed increase that is the only reason to use `link_planarity` applies only in the first round of post-processing. Hence, to combine the files into a single master output, you would use:

```
link_purify -imgs outim_mydata.txt -pairdets pairdets_mydata.csv -lflist clusterlist_mydata_NEOF_01b \
-useorbMJD 1 -outsum LPsum_mydata_NEOF_01b.csv -clust2det LPclust2det_mydata_NEOF_01b.csv
```

Here, -useorbMJD is a specific option that improves results for second-round post processing, as we explain in more detail below.

#### Customized `heliolinc` hypothesis files ####

You should feel free to make your own heliocentric hypothesis files rather than using the ones we have provided. Remember that the velocity units are AU/day, and the "radial acceleration" units are -GMsun/r^2. Unless you want to probe interstellar trajectories, be sure to calculate the solar escape velocity at each distance you probe, and do not include radial velocities with absolute value greater than this escape velocity.

### Files output by heliolinc ###

Two output files are produced for each execution of `heliolinc`: the "summary" file that provides a one-line summary of each cluster/linkage; and the clust2det file that maps clusters to individual detections in the input paired detection file. These are both csv-formatted files. They are complementary, and the post-processing programs such as `link_purify` require both files as their input.

#### The cluster summary file ####

The file providing summary lines for each cluster/linkage has a one-line header and 31 columns, as follows:

Column &nbsp;&nbsp;| Name               | Description
:---               | :---               | :---
1                  | clusternum         | Integer specifying which candidate linkage ('cluster') this is 
2                  | posRMS             | RMS scatter of position state vectors for this linkage (km)
3                  | velRMS             | RMS scatter of velocity state vectors for this linkage, scaled by characteristic time (km)
4                  | totRMS             | Total RMS scatter for this linkage in 6-D parameter space (km)
5                  | astromRMS          | Astrometric residual RMS from best-fit orbit (arcsec, written only in post-processing)
6                  | pairnum            | Number of (potentially overlapping) pairs or tracklets included in this linkage 
7                  | timespan           | Total time spanned by detections in this linkage (days )
8                  | uniquepoints       | Total number of unique detections in this linkage
9                  | obsnights          | Number of distinct observing nights represented in this linkage
10                 | metric             | Value of a quality metric for this linkage
11                 | rating             | Is this linkage pure or mixed, based on input string IDs of constituent detections
12                 | reference_MJD      | Main reference MJD used by initial run of heliolinc
13                 | heliohyp0          | Heliocentric distance at the reference time for the hypothesis that led to this linkage (AU)
14                 | heliohyp1          | Heliocentric radial velocity for the hypothesis that led to this linkage (km/sec) 
15                 | heliohyp2          | Second time-derivative of heliocentric distance for the hypothesis that led to this linkage, (m/sec^2)
16                 | posX               | Mean X-component of position state vectors at the reference time (km)
17                 | posY               | Mean Y-component of position state vectors at the reference time (km)
18                 | posZ               | Mean Z-component of position state vectors at the reference time (km)
19                 | velX               | Mean X-component of velocity state vectors at the reference time (km/sec)
20                 | velY               | Mean Y-component of velocity state vectors at the reference time (km/sec)
21                 | velZ               | Mean Z-component of velocity state vectors at the reference time (km/sec)
22                 | orbit_a            | Semimajor axis of best-fit orbit, in AU (written only in post-processing)
23                 | orbit_e            | eccentricity of best-fit orbit (written only in post-processing)
24                 | orbit_MJD          | MJD at epoch of best-fit orbit (written only in post-processing)
25                 | orbitX             | X-component of position state vector for best-fit orbit, in km (written only in post-processing)
26                 | orbitY             | Y-component of position state vector for best-fit orbit, in km (written only in post-processing)
27                 | orbitZ             | Z-component of position state vector for best-fit orbit, in km (written only in post-processing)
28                 | orbitVX            | X-component of velocity state vector for best-fit orbit, in km/sec (written only in post-processing)
29                 | orbitVY            | Y-component of velocity state vector for best-fit orbit, in km/sec (written only in post-processing)
30                 | orbitVZ            | Z-component of velocity state vector for best-fit orbit, in km/sec (written only in post-processing)
31                 | orbit_eval_count   | Iterations used to converge to best-fit orbit (written only in post-processing)

##### Notes on specific, potentially confusing columns #####

The quality metric (column 10) has one definition when written be `heliolinc` and another definition, based on more detailed and definitive information, when re-written by the post processing codes. In both cases it is used to identify the best candidate among a set of mutually redundant or overlapping clusters. For `heliolinc`, the metric is equal to the total time span (column 7), multiplied by the number of unique points (column 8), multiplied by the number of distinct observing nights (column 9), divided by the total RMS of the cluster (column 4). For the post-processing codes, the total RMS of the cluster (in km) is replaced by the RMS astrometric residual from the best-fit orbit, which is much more diagnostic but is not available to `heliolinc`. The post-processing codes also provide options for the user to determine what quantities should be included in the metric, and to what powers they should be raised. These will be discussed further below, when we describe the command line arguments for the post-processing codes.

The heliocentric hypotheses (columns 13-15) have the same meaning as those in the input hypothesis file supplied to heliolinc using the command line keyword `-heliodist`, but they do **not** all have the same units. The distance hypothesis does have the same units (AU), but the velocity, which was in AU/day in the input hypothesis file, is converted into km/sec for the output file, and the deriviative of the velocity, which was given in units of the local GMsun/r^2 in the hypothesis file, is converted into m/sec^2. These physical units are likely more intuitive, and in the case of the velocity they match the state vector units. We use m/sec^2 rather than km/sec^2 for the derivative of the velocity in order to keep the numerical values from being inconveniently small.

The MJD at the epoch of the best fit orbit (colum 24), like all the columns from 22 to 31, is written only by the post-processing code, since `heliolinc` does not carry out an orbit-fit. This MJD at the epoch is chosen to lie at the midpoint of the time spanned by the observations corresponding to a particular linkage -- i.e., the time when the orbit fit is most strongly constrained by the data and the corresponding state vectors are likely to have the least error. This can differ from the overall reference MJD (column 12) because that is referenced to the time spanned by the data as a whole, while the observations corresponding to one particular object will in general span a smaller time interval and have a different central time. 

#### The cluster-to-detection (clust2det) file ####

The clust2det file is a simple, two-column CSV, giving the cluster index (i.e., clusternum from the summary file) in the first column, and the detection number (line number from the paired detection file) in the second column. Hence, it maps the final linkages produced by heliolinc and the post processing codes to the detections in the paired detection file -- from which, if necessary, they can be mapped back to the lines in the original input file via the origindex column in the paired detection file. Note that both clusternum and detnum start counting from 0.

This concludes the description of file formats used by programs in the heliolinc C++ suite.

## Taking Control: understanding the optional arguments ##

All the programs in the heliolinc C++ suite allow many optional arguments beyond the minimal invocations used above. In many cases these optional arguments are the key to getting the specific behavior you want from the programs. 

### Arguments for make_tracklets ###

The table below describes all of the arguments 

keyword     | type         | Status      | Units      | Default           | Description
:---        | :---:        | :---:       | :---       | :---              | :---
-dets       | file name    | Required    | NA         | none            | CSV-formatted input file containing the detection catalog
-imgs       | file name    | Optional    | NA         | none              | Input image catalog
-outimgs    | file name    | &nbsp;Recommended&nbsp; | NA      | outimfile01.txt   | Output image catalog
-pairdets   | file name    | Recommended | NA         | pairdetfile01.csv | Output paired detection catalog
-tracklets  | file name    | Recommended | NA         | trackletfile01.csv&nbsp; | Output file recording tracklets
-trk2det    | file name    | Recommended | NA         | trk2detfile01.csv | Output file mapping tracklets to paired detection catalog
-colformat  | file name    | Critical    | NA         | See below         | Tells what columns in the input detection catalog hold each of the required values
-imrad      | float        | Optional    | degrees    | 2.0               | Center-to-corner size of images
-maxtime    | float        | Optional    | hours      | 1.5               | Maximum time between detections in a pair
-mintime    | float        | Optional    | hours      | 0.000278 (1 sec)  | Minimum time between detections in a pair
-maxGCR     | float        | Optional    | arcsec     | 0.5               | Maximum Great Circle residual for a tracklet with more than two points
-mintrkpts  | integer      | Optional    | none       | 2                 | Minimum number of points for a valid tracklet
-max_netl   | integer      | Optional    | none       | 2                 | Maximum number of points for a non-exclusive tracklet
-minvel     | float        | Optional    | deg/day    | 0.0               | Minimum angular velocity for a valid tracklet
-maxvel     | float        | Optional    | deg/day    | 1.5               | Maximum angular velocity for a valid tracklet
-minarc     | float        | Optional    | arcsec     | 0.0               | Minimum angular arc for a valid tracklet
-earth      | file name    | Required    | km, km/sec | none              | Heliocentric ephemeris for Earth, in CSV format, from JPL Horizons 
-obscode    | file name    | Required    | NA         | none              | Maps MPC observatory codes to Earth latitude, longitude, and elevation
-forcerun   | on/off flag  | Optional    | NA         | Off               | If set, pushes through all non-fatal errors
 
#### A few more details about previously discussed arguments ####

The column formatting file is technically not required, because `make_tracklets` defaults to the following column specifications:

```
IDCOL 1
MJDCOL 2
RACOL 3
DECCOL 4
MAGCOL 5
BANDCOL 6
OBSCODECOL 7
```
However, in the very likely event that you have constructed your input detection catalog with additional useful data beyond the minimum required columns and/or with a different column order than the above, you will need to create a new column formatting file (using the above as a template) and supply it using the `-colformat` command line keyword. This is why we characterize the column formatting file as "Critical" in the table above -- meaning you can run without it, but that is not recommended.

As previously discussed, the input image file is optional. When no input image file is supplied `make_tracklets` will construct one internally, which will be written out for use by `heliolinc` and the post-processing programs.

The names of the output image file, and the output paired detection file, though technically optional because of the default names describe above, should always be specified.

The only keyword that does not require an associated argument is the on/off flag "-forcerun". If the string "-forcerun" is in the argument list, `make_tracklets` will run in a mode that ignores all non-fatal errors. The errors to be ignored in this context have to do with missing information from the input detection catalog and the column formatting file. For example, `make_tracklets` will run with an input file that does not include the observatory code corresponding to each detection (it will assume they are all geocentric), but only if the -forcerun flag is set. Without -forcerun, `make_tracklets` will exit with an error message when it determines that the observatory codes are missing.

#### Purpose and use of the other arguments ####

Here, we will attempt to describe how you can use the arguments that haven't yet been discussed to fine-tune the behavior of `make_tracklets` and obtain the desired results for your particular science objective.

**Image radius** `-imrad`: This is nominally the center-to-corner distance, in degrees, for each of your images. For example, with a square field 1 degree on a side, it would be 0.707 degrees. The image radius is used when `make_tracklets` loops over all images, considering each in turn as 'Image A' and constructing a list of possible 'Image B' images whose detections might form pairs or tracklet with detections on Image A. A viable Image B must have a central RA and Dec within twice the image radius (plus an additional margin accounting for the movement of the objects) of the central RA and Dec of image A. Hence, setting the image radius too small can prevent the formation of valid pairs, while setting it too large wastes compute time probing 'candidate' pairs on images that are really too far apart on the sky. Note that if no input image file is provided, `make_tracklets` estimates the central RA and Dec of each image through an average of detections made at the same time from the same observing site. If the average number of detections per image is small, the estimated image center may be quite inaccurate. The possible bad effects of this can be mitigated by using a nominal image radius larger than the true value.

**Maximum time interval** `-maxtime`: This parameter again relates to identifying 'Image B' candidates whose detections might pair with those of a given Image A. If the images are separated by more than the specified maximum time, in hours, no pairs will be sought. This helps avoid excessive numbers of pairs, which can greatly increase the computational load both for `make_tracklets` itself, and for `heliolinc` when it is run on the paired detection files produced by `make_tracklets`.

**Minimum time interval** `-mintime`: This parameter defaults to 1 second (0.000278 hours). It can be raised to a larger value to avoid making tracklets spanning so little time that, for many objects, negligible motion will be detected. Typically, `heliolinc` cannot link objects with negligible motion because the velocity is essentially undetermined -- so eliminating such objects to begin with could reduce computation time. We have typically not used this option, but it might be valuable when targeting very slow-moving objects. In many cases, the minimum arc length (see below) is likely to be a better way to achieve the same effect.

**Maximum Great Circle Residual** `-maxGCR`: This parameter takes effect only in cases where `make_tracklets` detects overlapping pairs that might indicate more than two detections of the same object, and hence attempts to construct a 'tracklet' with more than two detections. It attempts to fit such a candidate tracklet using a Great Circle trajectory on the sky, at constant angular velocity. If the RMS residual from this Great Circle fit is less than the value specified (in units of arcseconds), the tracklet is considered valid and the pairs are merged. If the residual is too high, the tracklet is considered invalid and its constituent pairs are written to the output file individually. Note that unlimited sharing of detections between different pairs is allowed, but tracklets are entirely exclusive. In other words, if a set of detections are found to make up a valid tracklet, none of them will be included in any other pair (or tracklet). The maximum Great Circle residual defaults to 0.5 arcsec. For surveys with 1-2 arcsecond pixels, larger values up to 2 arcseconds can be reasonable. If for some reason it is desired that **only pairs** should be generated (no tracklets with more than two points), one way of ensuring this is to set the maximum Great Circle residual to an unreasonably small value such as 0.0001 arcseconds.

**Maximum non-exclusive tracklet length** `-max_netl`: For tracklets with only two detections (i.e., pairs), `make_tracklets` always allows unlimited overlap. In other words, a given detection can be paired with an unlimited number of other detections, and each such pairing will be passed to `heliolinc` as a distinct tracklet to be considered. However, by default `make_tracklets` does not extend this allowance of overlap to tracklets with more than two points. If a tracklet is constructed with three or more points lying along a consistent Great Circle trajectory, `make_tracklets` does not allow any of these points to be included in a different pair or tracklet. This default behavior was chosen because some of our test data included tracklets with very large numbers (~60) distinct detections. If such a long tracklet is non-exclusive, there are 1770 possible distinct two-point subtracklets, 34220 three-point subtracklets, and 487635 four-point subtracklets -- a time-consuming and utterly useless computational task for `heliolinc`, since the single 60-point tracklet encodes all the relevant information with higher precision.

Clearly, a 60-point tracklet has to be exclusive. But making three-point or four-point tracklets non-exclusive is not necessarily such a bad idea. For example, a non-exclusive three-point tracklet has only three possible sub-tracklets, each with two points -- and a non-exclusive four-point tracklet has only ten possible sub-tracklets (6 with two points and 4 with three points). The `-max_netl` parameter enables the user to instruct `make_tracklets` to accept such overlapping tracklets (using, respectively, -max_netl 3 or -max_netl 4). The point is that there really can exist crossing or overlapping tracklets, and it can be worthwhile to incur some amount of duplication in the simpler cases in order to avoid missing them. In tests on simulated data, going from the default (`-max_netl` 2) to `-max_netl` 4 increased the completeness of the `heliolinc` suite as a whole from 98.52% to 98.70%. This was for LSST data, where the aim is only two images per field per night -- meaning that 4-point tracklets are relatively rare. The gain might be considerably larger for a survey taking more image per field per night.

**Minimum number of points in a tracklet** `-mintrkpts`: This defaults to 2, and values smaller than 2 are of course meaningless. Any value larger than 2 disallows pairs, and will in general greatly reduce the size of the output 'pair file' and hence the runtime of `heliolinc` when this file is used as input. In the context of LSST, the `heliolinc` suite of programs is explicitly designed to enable multi-night linking of pairs, and setting the `-mintrkpts` to 3 or 4 is self-defeating. For surveys or data sets that have more than two images per field per night, however, such settings can greatly reduce the number of candidate objects, the resulting `heliolinc` runtimes, and the purity of the output.

**Minimum angular velocity** `-minvel`: This is the minimum angular velocity, in degrees per day, for a valid pair or tracklet. It defaults to zero, but can be used to exclude extremely slow-moving candidate detections. One reason you might want to do this is if you find that vast numbers of spurious slow-moving candidates are being produced due to false-positive detections caused by stationary stars.  As with the minimum time interval, however, the minimum arc length (see below) may provide a better way to achieve the same effect.

**Maximum angular velocity** `-maxvel`: This defaults to 1.5 degrees per day, and should be set with care because it has an extremely strong effect on the number of pairs and tracklets produced. For a given detection, the area of sky containing other detections that could be pair-partners after a given interval of time scales with the square of the maximum angular velocity. Hence, setting large values of the `-maxvel` can result in a huge increase in spurious pairs. `heliolinc` can handle this, up to a point, but runtimes and file sizes may be unnecessarily increased. An important consideration is at what angular velocity the detections would be noticeably trailed. If trailed-source information is available for the input data, it might make sense to run `make_tracklets` two or more times with different target regimes. For example, run it once with `-maxvel` set to the lowest angular velocity where detections are expected to be clearly trailed -- and then create a new input catalog culled down to **only** the trailed detections, and run `make_tracklets` again with a much faster maximum angular velocity. Note that the margin for object motion, mentioned in the discussion of image radius above, is equal to the maximum angular velocity times the time separation of two images being considered.

**Minimum angular arc** `-minarc`: This defaults to zero, and has units of arcseconds. It refers to the minimum angular separation for the endpoints of a valid pair or tracklet. It is very useful for rejecting spurious pairs/tracklets caused by stationary stars. A value similar to that used for `-maxGCR` -- e.g., two or three times the astrometric precision on the faintest detectable objects -- will often be a sensible setting.

**Force run** `-forcerun`: As already mentioned, including the `-forcerun` in the argument list for `make_tracklets` will cause `make_tracklets` to push through all non-fatal errors, such as the failure of the input file/column formatting file to supply an observatory code or other important information for each detection. Because `-forcerun` is a simple on/off flag, it is not necessary to supply a value or filename following the keyword, as it is for the other arguments.

### Command-line arguments for heliolinc ###

The table below describes all of the arguments 


keyword          | type         | Status      | Units       | Default           | Description
:---             | :---:        | :---:       | :---        | :---              | :---
-imgs            | &nbsp;file name&nbsp;    | &nbsp;Required &nbsp;   | NA          | none              | space-delimited image catalog produced by `make_tracklets`
-pairdets        | file name    | Required    | NA          | none              | CSV-formatted paired detection catalog produced by `make_tracklets`
-tracklets       | file name    | Required    | NA          | none              | CSV-formatted tracklet catalog produced by `make_tracklets`
-trk2det         | file name    | Required    | NA          | none              | CSV-formatted, two-column tracklet-to-detection mapping produced by `make_tracklets`
-mjd             | float        | Optional    | days (MJD)  | midpoint of timespan | Reference time, in Modified Julian Days (MJD), to which the state vectors from all tracklets will be integrated
-autorun         | 0 or 1       | Optional    | none        | 1                 | -autorun 0 turns off the default setting of the reference MJD to the midpoint of the time spanned by the input observations, and causes -mjd to become a required argument.
-obspos          | file name    | Required    | km, km/sec  | none              | Heliocentric ephemeris for Earth, in CSV format, from JPL Horizons 
-heliodist       | file name    | Required    | AU, AU/day, GMsun/r^2 | none    | Space-delimited list of heliocentric hypotheses to be probed, where each line gives a hypothetical heliocentric distance, radial velocity, and "acceleration" (time derivative of radial velocity).
-clustrad        | float        | Optional    | km          | 1.0e5		| clustering radius used to find clusters in the space of state vectors that have been integrated to the reference MJD. Note: the clustering radius automatically scales with the geocentric distance in AU.
-clustchangerad  | float        | Optional    | AU          | 0.5               | A geocentric distance within which the clustering radius no longer decreases linearly with decreasing geocentric distance, but remains fixed at the value it had at clustchangerad.
-npt             | integer      | Optional    | none        | 3                 | Minimum number of points (i.e., tracklets, or equivalently state vectors) that are required for a valid cluster.
-minobsnights    | integer      | Optional    | none        | 3                 | Minimum number of distinct observing nights required for a valid cluster
-mintimespan     | float        | Optional    | days        | 1.0               | Minimum total time that must be spanned for a valid linkage (cluster).
-mingeodist      | float        | Optional    | AU          | 0.1               | Minimum geocentric distance that will be probed. Tracklets whose inferred geocentric distances at the reference time (under a given hypothesis) are more than a factor of geologstep less than this will not be considered (under that hypothesis).
-maxgeodist      | float        | Optional    | AU          | 100.0             | Maximum geocentric distance that will be probed. 
-geologstep      | float        | Optional    | none        | 1.5               | Logarithmic scaling for bins of geocentric distance. The default value of 1.5 means that a bin centered on x AU will range from x/1.5 to x\*1.5, and the next (overlapping) bin will be centered on x\*1.5.
-mingeoobs       | float        | Optional    | AU          | 0.0               | Minimum inferred geocentric distance at time of observation (in contrast to -mingeodist, which applies at the reference MJD.
-minimpactpar    | float        | Optional    | AU          | 0.0               | Minimum inferred impact parameter w.r.t. Earth, if trajectory is extrapolated along a straight line.
-use_univar      | 0 or 1       | Optional    | none        | 0                 | For integrating Keplerian orbits, use the universal variable formulation, which is slightly slower, but can handle unbound (i.e., interstellar) trajectories.
-vinf            | float        | Optional    | km/sec      | 0.0               | Maximum v_infinity w.r.t. the sun. To probe interstellar trajectories, you must choose -vinf to be greater than zero **and** set -use_univar to 1.
-outsum          | file name    | Optional    | NA          | sumfile_test.csv  | Output cluster summary file, CSV formatted, one line per cluster
-clust2det       | file name    | Optional    | NA          | clust2detfile_test.csv | CSV-formatted, two-column cluster-to-detection mapping.
-verbose         | integer      | Optional    | none        | 0                 | Level of verbosity for status output



#### Purpose and use of arguments not previously described ####

**Reference time** `-mjd`: This is the reference time to which the inferred state vectors corresponding to each tracklet are propagated. By default, it is set to the midpoint of the time spanned by the data -- for example, if the earliest detection in your input detection catalog was made at MJD 60002.1 and the last detection was made at MJD 60012.2, the reference time will default to MJD 60007.15. This auto-generated reference time is rounded to two decimal places, to make it easier to remember and type if necessary. If you enter the keyword -mjd, whatever number you type after it overrides the default reference time, unless it is invalid. Any reference time that does not lie within the range spanned by the data is invalid: that is, the reference time cannot be before your earliest input detection, nor after your latest detection.

**autorun** `-autorun`: This is a simple switch (set to either 0 or 1) that controls whether a default reference time is in fact generated if the user does not supply one. `-autorun` defaults to on (`-autorun 1`). If you turn it off by entering `-autorun 0`, the program will exit with an error unless you explicitly supply a reference time. This option is provided in case you want to avoid ever using an automatically generated reference time.

**Clustering radius** `-clustrad`: This is the clustering radius in units of kilometers for identifying candidate linkages after propagating the orbits to the reference time. Recall that for each hypothesis about the heliocentric distance as a function of time, the heliolinc algorithm converts each pair or tracklet into a 3-dimensional position and velocity 'state vectors' which completely define an orbit. This orbit is propagated to the reference time using the Keplerian two-body approximation, producing 3-D position and velocity state vectors at the reference time. Hence, each input pair or tracklet, under each heliocentric hypothesis, becomes a distinct point in the 6-dimensional parameter space of position and velocity state vectors at the reference time. Linkages are identified in this 6-D parameter space using a simple range query on a KD-tree. Previous versions of `heliolinc` have used the DBSCAN algorithm rather than a KD range query, but DBSCAN, though more sophisticated, actually turns out to be demonstrably inappropriate for the specific use case of `heliolinc`. Changing to a KD range query has produced significant gains in completeness, and simultaneously a large reduction in the number of heliocentric hypotheses that must be used.

Prior to clustering, the velocity state vectors are converted to position units through multiplication by a characteristic time equal to one fourth the total time spanned by the input data.  For example, if the input data span 14.0 days -- the specified linking range for LSST -- the characteristic timecale is 14.0/4.0 = 3.5 days, and the velocity state vectors are multiplied by 3.5 days to convert them into distance units. The motivation for the one-fourth multiplier is the fact that, for observations uniformly distributed in time, the average distance from the central time is one-fourth of the total time span.

An alternative version of `heliolinc`, called `heliolincRR`, does not use the velocity state vectors explicitly. Instead, it integrates every tracklet to two different auxiliary reference times, one characteristic timecale before and after the main reference time. Clustering is performed in the 6-D parameter space of X,Y,Z position at the first auxiliary reference time and X,Y,Z position at the second auxiliary reference time. This idea was originally implemented by Ben Engebreth. It has the elegant property of producing a parameter space with naturally homogeneous units (while still effectively encoding the velocity via the separation between the two position vectors). In most tests, it also produces slightly better performance: hence, it is likely to be adopted for production with LSST.

Internally, `heliolinc` scales the clustering radius linearly with the (inferred, hypothesis-dependent) distance from Earth. To enable this, it divides the state vectors under each heliocentric hypothesis into overlapping, logarithmically spaced and sized bins in terms of geocentric distance. Clustering is performed independently in each bin, with a clustering radius equal to the specified values (default 100,000 km) times the bin-center geocentric distance in AU. Hence, a smaller clustering radius is used for objects close to the Earth. This distance-dependent clustering radius was added to resolve a problem encountered in testing where spurious clusters comprising unreasonable numbers (e.g. thousands) of tracklets were constructed at small geocentric distances. The linear scaling with distance was suggested by the fact that the original input data effectively correspond to angles on the sky as observed from Earth, and physical length at a given angular scale has the same linear dependence on distance.

Note that the post-processing programs `link_purify` and `link_planarity` also have a parameter `-maxrms`, which is analogous to the clustering radius and also defaults to 100,000 km.

**Radius at which the cluster-scaling changes** `-clustchangerad`: This is a minimum geocentric distance at which the linear scaling of the clustering radius, described above, still applies. At geocentric distances smaller than `clustchangerad`, the clustering radius remains at a fixed minimum, equal to the scaled value it had at `clusterchangerad`. This is to prevent the clustering radius from becoming excessively small (e.g., only 1000 km) at very small geocentric distances. Clustering radii this small require unreasonably fine sampling of the heliocentric hypotheses, leading to an excessive number of hypotheses and unacceptably long runtimes. In extreme cases, the clustering radius could even become smaller than the error on the Keplerian two-body approximation used to propagate the orbits, in principle making it impossible to find any asteroids regardless of the fineness of the hypothesis sampling.

Detection of objects very close to Earth with `heliolinc` is a difficult optimization problem, and `-clustchangerad` as well as the clustering radius, and three other parameters discussed below (`-mingeodist`, `-mingeoobs`, and `-minimpactpar`), may have to be carefully tuned. If `-clustchangerad` is set too small, the clustering radius for objects very near Earth will be too small and either none will be found or the hypothesis sampling will have to be so fine that many millions of hypotheses are probed and the computational cost is enormous. If `-clustchangerad` is too large, the clustering radius near Earth will be too large and the clusters will be large, impure mixtures of intractably confused detections of different asteroids and/or fully spurious detections. More on how to avoid this below.

**Minimum number of cluster points** `-npt`: The minimum number of points 'npt' in a cluster. Note that a 'point' here means a point in 6-D parameter space that originated as a pair or tracklet. Hence, the default value npt=3 means a valid linkage must consist of at least 3 *pairs*, which means at least 6 observations in all. The LSST specifications require npt=3. It is possible to set npt=2, but for orbital/geometric reasons, this tends to produce unreasonable numbers of false positives. Using npt>3 can greatly increase the purity of output linkages. It also reduced sensitivity and completeness relative to npt=3, but for densely sampled input data (e.g., with coverage of overlapping regions of the sky every night) the loss may not be severe.

**Minimum number of observing nights** `-minobsnights`: This is the minimum number of distinct nights on which detections must exist for a valid linkage. It is not redundant with `-npt` because some objects could have multiple pairs/tracklets on the same night -- hence, the default of minobsnights=3 imposes further constraints not already implied by npt=3.

**Minimum time span for a valid linkage** `-mintimespan`: This is the minimum total time spanned from the first to the last detection in a candidate linkage. It has units of days, and defaults to 1.0 -- a value small enough that it generally does not produce a separate constraint not already covered by the default minobsnights=3. In general, however, the two constraints are different: for example, you could set mintimespan=10 days with minobsnights=3, and then you would only get objects seen on at least three distinct nights for which the first night and the last night were separated by at least 10 days. Reasons why you might want to use the constraint include the fact that longer time spans result in more accurate orbits.

**Minimum geocentric distance** `-mingeodist`: This is the central distance for the first of the overlapping bins in geocentric distance. It has units of AU, and defaults to 0.1 AU. If no linkages are found in a given geocentric bin, `heliolinc` skips to the next bin with no problems or error messages.

**Maximum geocentric distance** `-maxgeodist`: The central distance for the outermost of the geocentric distance bins is guaranteed to be equal to or greater than this value. Each geocentric bin is a factor of `-geologstep` (see below) farther out than the one before, until on the final bin the central distance exceeds `-maxgeodist`, when no more bins are generated. The units of `maxgeodist` are AU, and the default is 100 AU. Again, the existence of geocentric bins that may not contain any points is benign.

**Logarithmic step size (and half-width) for geocentric bins** `-geologstep`: This unitless parameter defaults to 1.5. It means that the bin with central geocentric distance 'x' extends from x/1.5 to x\*1.5, and that the next bin outward will be centered on x\*1.5. Redundant linkages produced by the overlapping bins are harmless, and are naturally weeded out by the post-processing codes. The reason it's necessary for the bins to overlap is that otherwise a cluster split across a bin boundary could be incorrectly missed.

**Min. geocentric distance at observation** `-mingeoobs`: While the geocentric bins defined by `-mingeodist`, `-maxgeodist`, and `-geologstep` refer to geocentric distance after a tracklet has been integrated to the reference time, `-mingeoobs` refers to the inferred distance at the instant of observation. The parameter exists so that tracklets with extremely small inferred geocentric distances (under a specific hypothesis) can be rejected from further consideration. While obviously dangerous -- highly interesting discoveries could be discarded -- there is a compelling reason why rejecting such tracklets (under specific hypotheses) actually increase our ability to find asteroids. To see why, remember that each heliocentric hypothesis probed by `heliolinc` is simply a model of how the heliocentric distance changes with time. At a given instant in time, a given hypothesis proposes that any object observed just then is somewhere on a spherical shell, centered on the sun, with radius equal to the hypothesized heliocentric distance at that time. Heliolinc infers tracklet positions from the hypothesis by projecting the tracklet's RA, Dec coordinates outward from the observatory until they intersect this heliocentric sphere. Suppose that on a given night (or even shorter period of time), a given heliocentric hypothesis proposes a distance just slightly (say, 5000 km) more distant from the sun than the observatory itself. Then (since most observations are made in directions away from the sun) nearly all the tracklets from that night will have inferred (and very incorrect!) positions within a few thousand kilometers of the observatory -- and of one another. Thus, regardless of their actual distance, the inferred 3-D positions of tracklets from that night under that hypothesis will form a tight cluster very near the Earth. Such a spurious cluster can easily contain tens of thousands of points -- and (see below) it is likely to remain tight even when the tracklets have all been integrated to the reference time.

**Min. inferred impact parameter w.r.t. Earth** `-minimpactpar`: This is the minimum inferred impact parameter with respect to the Earth, when a tracklet's inferred position (at the time of observations) is extrapolated in a straight line at the (constant) inferred velocity. In the case of a large, spurious cluster such as described immediately above, the inferred impact parameters are likely to be even smaller than the inferred geocentric distances at time of observation. This is because, due to their small inferred distances from the observer, the tracklets' real angular velocities will correspond to unreasonably small inferred tangential velocities relative to the observer, and hence their inferred velocity will be almost entirely due to the heliocentric hypothesis and the observer's own motion with Earth's orbit. Since most detections are made in the direction approximately opposite the sun, a radial velocity away from the sun will correspond to a very small impact parameter relative to Earth -- considerably smaller than the tracklets inferred geocentric distance at time of observation. Meanwhile, a **real** asteroid observed at a very small geocentric distance would have an appreciable tangential velocity and (most likely) a larger impact parameter.

Because of this, `-mingeoobs` and `-minimpactpar` can be used together to eliminate spurious mega-clusters with very little risk of rejecting real asteroids. This works because an tracklet is rejected only if its inferred geocentric distance at observation is less than `-mingeoobs` **and** its inferred impact parameter is less than `-minimpactpar`. Tracklets with geocentric distances less than `-mingeoobs` but larger impact parameters would be retained.

The reason it is important to reject spurious mega-clusters close to the Earth is the fact that, as already mentioned, the tracklets in such a cluster all have very similar velocities. This means that they remain clustered -- in all six dimensions -- even when integrated to the reference time, and are likely to be identified as a bona fide candidate linkage. As a result, `heliolinc` and the post-processing codes spend useless computations analysing huge and completely implausible clusters, potentially missing real discoveries in the process.

**Use universal variables** `-use_univar`: Heliolinc's default Keplerian integrator uses the f and g functions, as derived in J. M. A. Danby's *Fundamentals of Celestial Mechanics*. This method, though very efficient, cannot handle unbound (i.e., interstellar, or hyperbolic) orbits. An alternative formulation, also from Danby's book, uses universal variables and handles bound and unbound orbits equally well. Since this formulation is slightly slower, we did not make it the default -- however, it works fine and should be invoked with `-use_univar 1` whenever `heliolinc` is used to search for interstellar objects. In this case it is also necessary to set `-vinf` (see below) to a value greater than the default, 0.0 km/sec.

Since its initial definition, `-use_univar` has been given additional possible values (from 2 to 15, besides the values of 0 and 1 already discussed). These enable the user to choose between various different "flavors" of `heliolinc` -- meaning, permutations of different algorithms for various stages of the analysis. The discussion of all these flavors is too involved to be included here: see the separate section on `heliolinc` flavors below.

**Maximum v_infinity w.r.t. the sun** `-vinf`: In the case of an interstellar object, v_infinity is the velocity it has in the limit of infinite distance from the sun -- in other words, the velocity it had (relative to the sun) long before the encounter when it was essentially unaffected by the sun's gravity, and the velocity it will have again (neglecting planetary encounters and other complications) when it has departed the sun's gravitational influence after passing through the solar system. If E is the object's total energy (kinetic plus gravitational potential), then v_infinity is the velocity it will have when all the energy is kinetic -- i.e., v_infinity = sqrt(2E/m), where m is the object's mass. This value was about 26 km/sec for 'Oumuamua and 32 km/sec for Comet Borisov. It could potentially be larger, but the odds of an interstellar object with v_infinity greater than, e.g., 200 km/sec (which would have to come from a star in the Galactic halo population) are rather low. Hence, `heliolinc` provides the `-vinf` parameter to enable the user to avoid wasting computation on highly improbable trajectories. If you are searching for interstellar objects, you **must** explicitly supply a value for `-vinf` greater than the default 0.0 km/sec, as well as setting `-use_univar` to 1.

We have also provided the somewhat strange option of specifying a negative value for v_infinity. For bound orbits with (by definition) negative total energy E, `heliolinc` interprets v_infinity as -sqrt(-2E/m). Hence, setting a negative maximum value for v_infinity will reject objects in "barely bound" solar orbits, such as long period comets in the inner solar system. This option might possibly be useful to reduce computation time in the case of a large, confusion-prone data set being searched for a specific population such as main-belt asteroids or Jupiter Trojans, which are never in "barely bound" orbits.

**Verbosity level** `-verbose`: The default of 0 writes a modest number of progress updates to stdout (i.e., to the terminal). For more detailed information, try `-verbose 1`. On the other hand, `-verbose -1` will cause the program to run in eerie silence.

### Flavors of heliolinc ###

This table gives the flavors of heliolinc: i.e., various permutations of algorithms that are described in detail below it.

`-use_univar` | Name        |Parameter space            | Clustering algorithm | Keplerian propagator | Hypothesis interpretation
:---:         | :---:       | :---                      | :---:                | :---                 | :---     
0             | heliolinc3D | pos & vel at 1 ref. time  | kd-tree range query  | f and g functions    | Keplerian
1             | heliolinc3D | pos & vel at 1 ref. time  | kd-tree range query  | universal variables  | Keplerian
2             | heliolincRR | position at 2 ref. times  | kd-tree range query  | f and g functions    | Keplerian
3             | heliolincRR | position at 2 ref. times  | kd-tree range query  | universal variables  | Keplerian
4             | heliolincR  | position at 1 ref. time   | kd-tree range query  | f and g functions    | Keplerian
5             | heliolincR  | position at 1 ref. time   | kd-tree range query  | universal variables  | Keplerian
6             | heliolinc3D | pos & vel at 1 ref. time  | DBSCAN               | f and g functions    | Keplerian
7             | heliolinc3D | pos & vel at 1 ref. time  | DBSCAN               | universal variables  | Keplerian
8             | heliolinc3D | pos & vel at 1 ref. time  | kd-tree range query  | f and g functions    | Taylor Series
9             | heliolinc3D | pos & vel at 1 ref. time  | kd-tree range query  | universal variables  | Taylor Series
10            | heliolincRR | position at 2 ref. times  | kd-tree range query  | f and g functions    | Taylor Series
11            | heliolincRR | position at 2 ref. times  | kd-tree range query  | universal variables  | Taylor Series
12            | heliolincR  | position at 1 ref. time   | kd-tree range query  | f and g functions    | Taylor Series
13            | heliolincR  | position at 1 ref. time   | kd-tree range query  | universal variables  | Taylor Series
14            | heliolinc3D | pos & vel at 1 ref. time  | DBSCAN               | f and g functions    | Taylor Series
15            | heliolinc3D | pos & vel at 1 ref. time  | DBSCAN               | universal variables  | Taylor Series

#### Comparing the different clustering parameter spaces ####

Three different options are available for the parameter space in which the clustering is performed.

The first, default, and most thoroughly explored option is the heliolinc3D parameter space pioneered by Siegfried Eggl. Here, the state vectors derived from each tracklet are integrated (using the Keplerian 2-body approximation) to a single reference time (which can be specified by the user as `-mjd` or allowed to default to the center of the time spanned by the input data). Clustering is performed in the 6-D parameter space of position and velocity at the reference time. The velocity vectors are converted to length units through multiplication by a characteristic timescale equal to 1/4 the total time spanned by the input data.

The second option is the heliolincRR parameter space, invented by Ben Engebreth. Here, the state vectors derived from each tracklet are integrated to two different reference times. The two reference times are before and after the central time by an interval equal to 1/4 the total time spanned by the data (note: Ben Engebreth has typically used times more closely spaced than this). Clustering is performed in the 6-D parameter space of 3-D position at the first reference time and 3-D position at the second reference time. This parameter space has the advantage of naturally homogeneous length units, with no conversion required.

The third option is the heliolincR parameter space, where just a single reference time is used and only the position vector is used for clustering (velocity is ignored). This is the only 3-D parameter space: the others are all 6-D.

Tests of heliolinc3D vs. heliolincRR suggest very similar performance, with different metrics favoring one or the other. There is evidence that heliolincRR works best at larger heliocentric distances, but it is not conclusive. HeliolincRR has delivered the highest completeness in tests on simulated main-belt asteroids and TNOs. However, experiments aimed at quantitatively optimizing hypothesis sampling suggest that heliolincRR requires finer hypothesis sampling for the same performance -- and there is some evidence even in the aforementioned simulations that heliolinc3D provides longer linkages (more total detections on average) than heliolincRR, even though the later provides a slightly larger number of linkages.

#### Comparing the different clustering algorithms ####

DBSCAN is much more sophisticated than the simple kd-tree range query that we prefer as an alternative. The advantage of DBSCAN is its ability to find clusters with odd shapes such as filaments, sheets, or shells, when they are extended well beyond the clustering radius. DBSCAN can only do this if the cluster contains many more than three points, however, and the ability to accumulate clusters with odd shapes and with sizes larger than the clustering radius also makes DBSCAN more likely to create impure, 'mixed' clusters containing tracklets from more than one object.

The asteroid-linking specification for LSST requires the discovery of any object that produces three (or more) unique tracklets in a two-week interval. This specification can be robustly met only if the hypothesis sampling, clustering radius, and other parameters are set such that any set of three tracklets corresponding to the same object and falling within the two-week time span will, for some hypothesis, be clustered tightly enough to fall within a single clustering radius of each other. If **any** set of **three** tracklets within a two-week time span will be clustered this tightly, regardless of how they are distributed in time, it follows that **any number** of tracklets that fall within the two-week time span (and all correspond to the same object) will fall within a single clustering radius. Hence, a simple kd-tree range query will find **the entire cluster**, regardless of how many points it has. This means that DBSCAN's elegant ability to trace an oddly-shaped yet contiguous cluster over a distance much greater than the clustering radius becomes unnecessary -- and indeed a liability because it increases the likelihood of mixed clusters. For these reasons, we prefer the kd-tree range query for linking asteroids in LSST.

This is not to say, however, that DBSCAN has no value in the `heliolinc` context. For one thing, the kd-tree range query finds a much larger number of mutually overlapping clusters (since a cluster is sought centered on every single point individually). This bloats up the output files and increase post-processing times. For LSST, the fact that it delivers greater completeness and purity more than compensates -- but in other contexts, the tradeoffs might be different. In particular, if every cluster is expected to include many more than three tracklets, then DBSCAN's ability to trace extended, oddly-shaped clusters might be valuable. There may be other advantages to DBSCAN in particular situations that we have not explored. It's a good idea to test multiple flavors of `heliolinc` when using it in a new context for the first time, since one of them might offer unexpected advantages.

#### Comparing different Keplerian propagators ####

As has already been discussed, Keplerian integration using the universal variables formulation of the Kepler problem is able to handle unbound (i.e., hyperbolic, interstellar) orbits, while the formulation using the f and g functions is not. We default to using the f and g functions because this makes the code run somewhat faster. However, there are contexts in which you might want to use the universal variables formulation instead. Even when not searching for strictly interstellar orbits, astrometric errors in the input data (especially for slow-moving objects such as TNOs) can cause a specific tracklet that really corresponds to a bound object to have a formally unbound trajectory, even for a well-chosen hypothesis. Such a tracklet will be rejected by the integrator using the f and g functions, but retained and likely successfully linked if using the universal variables (combined with a modest positive value for `-vinf`). The flip side, of course, is that admitting formally unbound tracklets will increase confusion, resulting in longer runtimes and a greater probability of false linkages.



### Command-line arguments for the post-processing codes ###

We will discuss `link_purify` first, as `link_planarity` has almost the same set of arguments.

The table below describes all of the arguments for `link_purify`.


keyword          | type         | Status      | Units       | Default           | Description
:---             | :---:        | :---:       | :---        | :---              | :---
-imgs            | &nbsp;file name&nbsp;    | &nbsp;Required &nbsp;   | NA          | none              | space-delimited image catalog produced by `make_tracklets`
-pairdets        | file name    | Required    | NA          | none              | CSV-formatted paired detection catalog produced by `make_tracklets`
-lflist          | file name    | Required    | NA          | none              | Each line gives the name of a summary file and of the corresponding clust2det file, separated by a space. The files can be directly output by `heliolinc`, or else produced by a previous run of post processing. Any number of lines is allowed.
-simptype        | integer      | Optional    | none        | 0                 | Type of simplex used to intialize orbit-fitting (0, 1, or 2). See below for more.
-rejfrac         | float        | Optional    | none        | 0.5               | Maximum fraction of points that can be rejected from a linkage.
-rejnum          | integer      | Optional    | none        | 50                | Maximum number of points that can be rejected from a linkage. Of rejfrac and maxrej, whichever results in the smallest number of rejected points is operative.
-max_astrom_rms  | float        | Optional    | arcsec      | 1.0               | Maximum RMS astrometric residual from the best-fit orbit. Points will be iteratively rejected until the RMS residual drops below this value, or a rejection limit is reached, or the linkage becomes invalid.
-minobsnights    | integer      | Optional    | none        | 3                 | Minimum number of distinct observing nights for a valid linkage.
-minpointnum     | integer      | Optional    | none        | 6                 | Minimum number of distinct points for a valid linkage.
-useorbMJD       | 0 or 1       | Optional    | none        | 0                 | For `-useorbMJD 1`, use the orbit_MJD and associated state vectors (orbitX, orbitY, etc.) in the summary file, if available, as the starting point for orbit fitting. Default is to use `reference_MJD` and its associated state vectors (posX, posY, etc.).
-ptpow           | integer      | Optional    | none        | 1                 | Power to which the number of unique points should be raised, in calculating the cluster quality metric
-nightpow        | integer      | Optional    | none        | 1                 | Power to which the number of distinct observing nights should be raised, in calculating the cluster quality metric
-timepow         | integer      | Optional    | none        | 0                 | Power to which the total time spanned by a linkage should be raised, in calculating the cluster quality metric
-rmspow          | integer      | Optional    | none        | 2                 | Power to which the RMS astrometric residual from the best-fit orbit should be raised, in calculating the cluster quality metric. Note that the metric is **divided by** the RMS raised to the selected power, rather than being multiplied by it as for the other quantities involved.
-maxrms          | float        | Optional    | km          | 1.0e5             | Maximum total RMS (column 4 from the cluster summary file) for a cluster to be considered. Clusters with larger total RMS will be rejected without consideration.
-outsum          | file name    | Optional    | none        | LPsumfile_test.csv | Name for the output cluster summary file
-clust2det       | file name    | Optional    | none        | LPclust2detfile_test.csv | Name for the output cluster-to-detection file
-verbose         | integer      | Optinal     | none        | 0                 | Level of verbosity for status output



#### Purpose and use of arguments not previously described ####

**Simplex Type** `-simptype`: The orbit-fitting employed by `link_purify` uses the Method of Herget to shrink the dimensionality of the problem from 6 to 2 dimensions, which become the geocentric distances at the times of the first and last observation being fit -- i.e., *geodist1* and *geodist2*. Within this 2-D parameter space, the downhill simplex method (see, e.g., Numerical Recipes in C) is used to find the best fit orbit by chi-square minimization. The -simptype keyword determines the characteristics of the 2-D simplex with which the orbit-fitting is initialized.

In 2-D, a simplex is just a triangle, and therefore it is defined by three distinct (*geodist1,geodist2*) points. The allowed values for `-simptype` are 0 (the default), 1, 2, or 3. In every case, one point of the simplex is simply the best value of *geodist1,geodist2* obtained either from heliolinc or optionally (if this is not the first round of post-processing) from a previous round of orbit fitting. This leaves only two points to be defined by -simptype. The scale of the simplex is set by a declared parameter, SIMPLEX_SCALEFAC = 0.2 (or *simpscale* for short), which the user cannot adjust -- however, simptype affects how this parameter is interpreted.

For -simptype 0, the size of the simplex scales with starting distances, so it will be larger for objects more distant from the Earth. The two additional points of this simplex are (*geodist1\*(1-simpscale),geodist2*) and (*geodist1,geodist2\*(1-simpscale)*). This will be an isoceles right triangle if *geodist1* = *geodist2*.

For -simptype 1, the simplex does **not** approximate an isoceles right triangle, but instead is deliberately elongated in the sense that geodist1 and geodist2 have relatively large fractional changes (+/- *simpscale*), while the ratio geodist2/geodist1 changes by smaller amounts related to *simpscale\^2* (it must change somewhat or the simplex would have zero area and not work). The motivation is that the direction in which this elongates the simplex -- that defined by a constant ratio of geodist2/geodist1 -- is often the least-constrained direction in the minimization problem.

For -simptype 2, the simplex does not scale with the initial distances but is instead an exact equilateral right triangle, with points (*geodist1-simpscale,geodist2*) and (*geodist2-simpscale,geodist1*). For -simptype 3, the simplex is a botched variant of -simptype 1 that produces extraordinarily bad results: hence, `-simptype 3` should never be used except as a demonstration that the definition of the initial simplex type matters. We have found `-simptype 1` to be somewhat better than the default in most case; `-simptype 2` has not been extensively tested.

**Maximum RMS astrometric residual from best-fit orbit** `-max_astrom_rms`: After each iteration of orbit fitting, `link_purify` calculates the RMS astrometric residual from the best-fit orbit, in arcseconds. If this exceeds `-max_astrom_rms`, the single largest astrometric outlier will be rejected, and `link_purify` will perform a new iteration of orbit fitting -- unless the rejection of this final outlier has made the linkage invalid (see below).

**Rejection fraction** `-rejfrac`: This determines the maximum fraction of the detections that can be iteratively rejected from the orbit fit. On each iteration, the most extreme astrometric outlier from the best-fit orbit is rejected, until the astrometric RMS residual drops below `max_astrom_rms`. If the specified fraction of points has already been rejected and the RMS astrometric residual still exceeds `max_astrom_rms`, the linkage is rejected without further analysis.

**Rejection number** `-rejnum`: This determines the maximum number of detections that can be iteratively rejected from the orbit fit. If `-rejfrac` and `-rejnum` are both set, the more conservative limit (fewer detections rejected) is operative. We provide both parameters to provide a way of handling very large spurious linkages. For example, if a linkage has 9000 points, it is almost certainly spurious in most contexts -- and yet with only `-rejfrac`, `link_purify` will go through 4500 computationally expensive iterations of orbit fitting in a vain attempt to clean up this intractable mess (for `-rejfac 0.5`, which is the default). The additional parameter `-rejnum` reduces the wasted compute by terminating the effort after only 50 iterations (for the default value `-rejnum 50`). If you have reason to think that for your particular data, very large linkages might be valid and worth purifying, you can set `rejnum` to some much larger value.

**Minimum number of observing nights** `-minobsnights`: Means exactly the same as for `heliolinc`. If the rejection of an astrometric outlier causes the remaining observations to be spread over a smaller number of nights, the linkage is considered to have become invalid and is therefore rejected without futher analysis.

**Minimum number of distinct detections** `-minpointnum`: Unlike the parameter `-npt` in `heliolinc`, this does not refer to tracklets but to individual detections. If the number of distinct detections in a linkage drops to `minpointnum`, and the astrometric RMS residual still exceeds `-max_astrom_rms`, the linkage is rejected without further analysis, since no more points are allowed to be rejected.

**Use orbit-fit MJD** `-useorbMJD`: For `-useorbMJD 1`, take the starting point for orbit fitting (initial MJD, *geodist1*, and *geodist2*) from the results of a previous orbit fit (`orbit_MJD`, `orbitX`, `orbitY`, etc., in the input cluster summary file), rather than from the raw estimate output by heliolinc (`reference_MJD`, `posX`, `posY`, etc.). If no previous orbit-fit exists, `link_purify` will detect this via the dummy values written for the orbit-fit parameters in the input cluster summary file, and `-useorbMJD 1` will not be operative: it will use the raw estimate from heliolinc, since that is all that is available. In other words, the default behavior, which corresponding to `-useorbMJD 0`, will occur unless `-useorbMJD 1` is explicitly set **and** valid results for the necessary orbital parameters are found in the input cluster summary file. 

**Point number power (exponent)** `-ptpow`: The power to which the number of distinct observations ('points') in a given linkage will be raised, when calculating the cluster quality metric. This defaults to 1 (quality metric scales linearly with number of points), but we have frequently set it as high as 3, with good results. Note that the definition of the cluster quality metric is very important. At the final stage of post-processing, for every set of mutually overlapping linkages only the one with the best quality metric will be retained: all the rest are permanently discarded.
1
**Night number power (exponent)** `-nightpow`: The power to which the number of distinct observing *nights* represented in a given linkage will be raised, in calculating the cluster quality metric.

**Time power (exponent)** `-timepow`: The power to which the total time spanned by a linkage will be raised, in calculating the cluster quality metric. This defaults to zero, not because we do not want to prioritize linkages that span long times if possible, but because this is better captured by the number of distinct observing nights (i.e., `-nightpow`) and, to a lesser extent, the number of points (`-ptpow`). All other things (including the number of individual points and distinct observing nights) being equal, a linkage than spans a longer time is **more** likely to be spurious than one spanning a briefer interval. This is because a longer time span produces more flexibility in orbit-fitting, increasing the odds that a formally acceptable orbit-fit might be found for a spurious linkage.

**RMS power (exponent)** `-rmspow`: The power to which the astrometric RMS residual (in arcseconds) is raised in calculating the cluster quality metric. Note that although a positive integer is expected (the default is 2), the quality metric is actually **divided** by the RMS raised to this power, so that a higher astrometric RMS produces a lower quality metric, as it should. 

**Maximum RMS of the original cluster** `-maxrms`: The maximum acceptable RMS, in km, for a cluster to be analyzed at all. Unlike the maximum astrometric RMS (`-max_astrom_rms`), this parameter refers to tracklets rather than individual observations, applies in the 6-D space where heliolinc performs its clustering, and (if exceeded), prevents any analysis of a given cluster rather than merely (like `-max_astrom_rms`) guiding the parameters of such an analysis. If `heliolinc` itself was run with sensible parameters, all you need to do with `maxrms` is make it a few times larger than the clustering radius used in `heliolinc`, and it will pass all of the clusters.

All of the above make it appear that `-maxrms` is a useless parameter, which is mostly true. However, if you ran heliolinc with what you later decided was too generous a clustering radius, and you want to analyze only the cleanest among an excessive number output clusters (at the cost of rejecting many that could possibly have been purified), you may want to set a smaller value of `-maxrms`. Otherwise, you should just make sure to set it large enough to pass all the clusters output by heliolinc -- an outcome which is guaranteed by making it larger than the clustering radius used in `heliolinc`.

**Verbosity level** `-verbose`: The default of 0 writes a modest number of progress updates to stdout (i.e., to the terminal). These include specific information on every thousandth cluster that is analyzed, to aid in monitoring of long runs. For more detailed information, try `-verbose 1`, which will output information on every cluster rather than just every thousandth cluster -- something that may be useful when diagnosing unexpected behavior using small input files. For even more status output (probably much more than you want), try `-verbose 2`. All larger values are equivalent to `-verbose 2`. Unlike `heliolinc`, `link_purify` does not offer an extra-silent option: setting `-verbose` to negative values simply produces the default behavior.

### Arguments specific to `link_planarity` ###

There is only one parameter specific to `link_planarity` and not shared by `link_purify`. It has far-reaching implications, however, as we now explain:

**Max out-of-plane distance** `-oop`: This parameter, which defaults to 1000 km but should usually be larger, sets the maximum distance from the mean plane for a valid observation. This is a computationally inexpensive proxy for the astrometric residual in arcseconds as used by `link_purify`. It is calculated by projecting each observation individually onto the heliocentric sphere defined by the hypothesis, thus mapping it to a 3-D position vector in the solar system. Since individual points rather than tracklets are used, there is no velocity inference and no orbital integration to a reference time: each point is considered at its own time-of-observation. The mean plane of the 3-D position vectors corresponding to all the points is calculated, and outliers with out-of-plane distance exceeding the `-oop` parameter are interatively rejected until none remain. The surviving points, all guaranteed to have inferred 3-D positions within a distance `-oop` of the mean orbit plane, are passed to the (much more computationally expensive) full Keplerian orbit fit. After this, `link_planarity` proceeds exactly like `link_purify`, iterative rejecting astrometric outliers from the best-fit orbit until the RMS residual of the survivors falls below `-max_astrom_rms`.

The reason `link_planarity` exists is because the pre-screening and rejection of outlying points using their out-of-plane distances reduces the number of full Keplerian orbit fits (which are much more computationally expensive) and speeds up the analysis, typically by a factor of several. This major speedup comes at a cost: a small but nonzero fraction of good linkages are incorrectly rejected. Hence, you should use `link_planarity` **only** if post-processing time is a limiting resource in your analysis (which, unfortunately, it is quite likely to be).

The keyword `-useorbMJD` is available in `link_planarity`, but should **not** be used: it should always be left to the default value of 0.

Relatedly, if the post-processing codes are run multiple times in series, `link_planarity` should be used **only for the first round** of post-processing. For subsequent rounds, since astrometric fitting has already been performed and outliers eliminated, there is no further advantage to `link_planarity`, and `link_purify` should be used instead. To summarize, use `link_planarity` **only** on raw output from heliolinc. 

## Converting to MPC 80-column format ##

The output detection file from heliolinc provides linkage membership, MJD, RA, Dec, magnitude, photometric band, and observatory code for all successfully linked detections. This enables easy conversion to the Minor Planet Centers dated but still useful 80-column format for representing asteroid observations. We have provided a program for doing so. It is not as mature as the heliolinc suite proper, but we hope it may be useful. Its command-line arguments, in order, are the paired detection file produced by `make_tracklets`, the output detection file from `heliolinc`, the output summary file from `heliolinc`, a prefix for temporary string identifiers in the 80-column format, and an output file. Somewhat confusingly, the `heliolinc` detection file is supplied with the command line keyword `-clust` and the summary file is supplied with the command line keyword `rms`. It can be run on files output directly by `heliolinc`, but since those usually have large numbers of overlapping linkages, we strongly recommend that **`cluster2mpc80a` should only be run on files post-processed by `link_refine`**. Hence the following example:

```
cluster2mpc80a -pairdet LSST_pairdets_01.csv -clust LSST_linkref_all01.csv \
-rms LSST_linkref_summary01.csv -prefix hl -outfile LSST_hl_01.mpc
```

The table below describes all of the arguments.

| command line keyword | type      | Required? | Description                                            |
| -------------------- | --------- | --------- | --------------------------------------------- |
-pairdet               | file name | yes       | csv-formatted input file produced by `make_tracklets` and containing the paired detection catalog |
-clust                 | file name | yes       | output file from `heliolinc` or (**strongly preferred**) `link_refine` giving all detections for each linkage |
-rms                   | file name | yes       | output file from `heliolinc` or (**strongly preferred**) `link_refine` giving one-line summary statistics for each linkage |
-prefix                | string    | yes       | short string prefix (only 2 characters recommended) for MPC 80-column formatted output file |
-outfile               | file name | yes       | name of output file |


If you get error messages, check very carefully that the input files are all of the correct format, are mutually compatible, and are introduced with the correct keywords. Note well that `cluster2mpc80a` is not recommended for use on output files directly from `heliolinc`, but only on the identically-formatted, refined and culled-down versions produced by `link_refine` and `link_refine_mulisite`.

The output files produced by `cluster2mpc80a` should be suitable for orbit-fitting with any program that accepts MPC 80-column format as input. They can also be pasted into the Minor Planet Checker(`https://minorplanetcenter.net/cgi-bin/checkmp.cgi`) to check for matches to known objects (remember to click the 'or around these observations' button on the web form).
