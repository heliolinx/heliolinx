// heliolinc_TNV: First prototype implementation using Ben
// Engebreth's concept of SAD (Swept Angle Discrepancy) tracklets.
//
// Implements in C++ (with some modifications) the Heliolinc3D
// algorithm developed by Siegfried Eggl, which in turn was based
// on the original HelioLinC algorithm developed by Matthew Holman
// (Holman et al. 2018, Astronomical Journal, 156, 135).
//
// The asteroid linking problem consists in identifying
// subsets of detections in an input catalog of astronomical
// transient detections that are in fact repeated detections of the
// same previously unknown asteroid. For example, by repeatedly
// observing a 1000 square-degree area of sky over a two-week
// period, an astronomical survey might identify millions of
// cases where an apparent source is detected at a position
// where there is no known star or asteroid. Each detection
// is uniquely identified by its celestial coordinates (e.g.,
// right ascension and declination, or RA and Dec) and the time
// of observation (e.g., modified Julian day, or MJD). These
// detections could be spurious (cosmic ray hits, electronic
// crosstalk, fluctuations of Poisson noise); stationary
// transients (e.g., supernovae or previously unknown
// variable stars); or asteroids. Although thousands of unknown
// asteroids may each have been detected many times, none can
// be said to have been 'discovered' until the linking problem
// has been solved and it is known which sets of detections are
// actually repeated measurements of the same unknown asteroid.
// Then the asteroid has really been discovered; its orbit can
// be calculated; future followup observations can be planned; and
// serendipitous past detections ('precoveries') may be identified.
//
// Heliolinc3D solves the asteroid linking problem by probing a
// number of hypotheses for the time-dependent heliocentric distance
// of undiscovered asteroids that may be present in the input data.
// It requires five input files: one is the detection catalog,
// one is a 'pair' file (see below), and one specifies the hypotheses
// to be probed. The remaining two files are reference files that
// will not typically need to be changed from one execution to another:
// one is a JPL Horizons ephemeris for the Earth, and the other is
// an observatory code file based on the one available from the Minor
// Planet Center at https://minorplanetcenter.net/iau/lists/ObsCodesF.html
// For digestion by the current version of heliolinc, the incomplete
// lines that relate to space-based observatories in this file must
// be deleted. Useful defaults for these last two files,
// Earth1day2020s_02a.txt and ObsCodes.txt, are available in the
// same github repository as this source code.
//
// The input detection catalog and preliminary 'pair' file required
// by heliolinc can be produced by using make_tracklets, available in the
// same github repository as the current code. The input to make_tracklets
// is a .csv detection catalog of very flexible formatting, which is
// described in more detail in the documentation for make_tracklets.
// The pair file lists pairs of entries in the detection catalog that
// might be detections of the same object within a relatively short
// span of time (e.g., 1.5 hours). These 'pairs' can also be 'tracklets'
// consisting of more than two detections, if the input catalog
// includes cases where more than two observations obtained close
// together in time lie along a consistent linear trajectory.
//
// The input file providing the hypotheses that should be probed, where
// each hypothesis proposes a particular dependence of unknown asteroids'
// heliocentric distance as a function of time, has one line per
// hypothesis, and four columns, as follows:
//
// 1. Heliocentric distance (AU)
// 2. Heliocentric radial velocity (AU/day)
// 3. Normalization (0: skip this line. Positive: don't skip).
// 4. 2nd time-derivative of heliocentric distance, scaled by -GMsun/r^2
//
// Columns 1, 2, and 4 effectively approximate the time-dependent
// heliocentric distance as a Taylor Series, providing the 0th,
// 1st, and 2nd time-derivatives of the heliocentric distance.
// The forth column can be thought of as an acceleration, but it's
// important to note that it is **not** the vector acceleration
// away from the sun, which always -GMsun/r^2 for any orbit. Instead,
// it's the second derivative of the distance from the sun, which is
// zero for a circular orbit, and reaches -GMsun/r^2 (1.0 in the
// normalized units used by heliolinc), only in the limit of an orbit
// with zero angular momentum (i.e., with eccentricty exactly 1.0)
//
//
//
//
// DEVELOPMENT NOTES IN REVERSE CHRONOLOGICAL ORDER:
//
// April 26, 2024: Uses Ben Engebreth's heliolincRR algorithm,
// which calculates positions at two different times, on either
// side of the master reference time, instead of calculating
// position and velocity at the reference time.
//
// Feb 22, 2024: Uses a simple KD range query rather than the DBSCAN
// algorithm, because the latter has problematic side effects.
//
// June 22, 2022: besides the viewing geometries
// handled by all earlier versions, also handles the previously-ignored
// geometry, possible only at solar elongation less than 90 degrees,
// where the line-of-sight from the observer intersects the
// the heliocentric sphere defined by the asteroid distance hypothesis
// on the near side of the sun -- that is, the asteroid is at the
// point where the line-of-sight vector pierces the heliocentric
// sphere from the **outside**. An alternative specification for
// this geometry is that the Sun-object-observer phase angle is
// **greater** than 90 degrees. It corresponds to a second solution
// of the quadratic equation for the intersection of a line with
// a sphere. This second solution was ignored in previous versions,
// but the current program handles **both** solutions correctly.
// Hence, it is not a supplement but a replacement for previous
// versions.
//
// March 15, 2022: reads the pairdets file in csv format
// with a header. Also replaces the
// daysteps parameter with obsnights = daysteps+1 (on output).
// The parameter daysteps is still used internally exactly as before.
//
// Note well that the format of the pair file (as opposed to the
// pairdets file) has NOT changed. The pair file has an unavoidably
// specialized format because of the need to handle both pairs and
// tracklets. This specialized format does not lend itself to the
// header+csv convention. Also, the pair file is not needed by any other
// program except projectpairs itself, so its specialized format
// is not expected to create compatibility issues for programs external
// to the heliolinc suite.
// 
// March 11, 2022: reads and keeps track of magnitude,
// photometric band, and observatory code for each observation.
//
// March 04, 2022: uses a clustering radius that
// depends on geocentric distance. The user-defined clustering
// radius will now be taken to apply to a standard geocentric
// distance of 1.0 AU.
//
// January 24, 2022: reads a more sophisticated pair
// file that aggregates tracklets into single effective pairs.
//
// January 07, 2022: uses an integerized form of the
// state vectors for the k-d tree and DBSCAN implementations,
// for speed.
//
// January 06, 2022: streamlined a bunch of kludgy
// aspects of the DBSCAN implementation that were good for debugging
// but not appropriate for production.
//
// January 04, 2022: aimed at performing the further
// step of clustering state vectors propagated to the reference
// time in order link asteroid discoveries. In support of doing
// this very fast in cases where we have millions of detections,
// the observers heliocentric coordinates will be supplied in
// the input paired detection file, along with string identifiers
// indicating which detections correspond to the same unique object,
// which are required for testing purposes. NOTE: I first successfully
// implemented the DBSCAN algorithm at this stage of developement
//
// December 07, 2021: Read pair files output by maketrack01a.cpp, and use a
// constant-acceleration model of the heliocentric distance.
// This is the first version that searches a grid
// of heliocentric distance and velocity. The grid is
// defined in a file whose name is supplied in the
// configuration file, which also supplies best-guess
// heliocentric accelerations for each heliocentric distance
// and radial velocity.
//
// Analysis performed:
// 1. Calculate the observer's precise barycentric
//    position at the moment of each observation.
// 2. Use an approximate distance for the target from
//    the Solar System barycenter at each observation to calculate
//    the target's barycentric coordinates from the observed RA, Dec.
// 3. Difference pairs of observations to obtain the target's
//    barycentric velocity.
// 4. Integrate the position and velocity information from each
//    pair of observations to calculate the object's location
//    at an input reference time. Hence, each pair of observations
//    is used to obtain an independent estimate of the object's
//    barycentric position and velocity at a fixed reference time
//    that is the same for all of the pairs.
// 5. The whole point of this exercise is to investigate the distribution
//    of these independently derived estimates in the 6-D space of
//    barycentric position and velocity. The tighter they cluster, the
//    easier the linking problem.

#include "solarsyst_dyn_geo01.h"
#include "cmath"

static void show_usage()
{
  cerr << "Usage: heliolinc_TNV -imgs imfile -pairdets paired detection file -tracklets tracklet file -trk2det tracklet-to-detection_file -autorun 0=do_not_auto-generate_MJDref -mjd mjdref -earth earthfile -heliodist heliocentric_dist_vel_acc_file -clustrad clustrad -clustchangerad min_distance_for_cluster_scaling -cluster_minpts cluster_minpts -minobsnights minobsnights -mintimespan mintimespan -mingeodist minimum_geocentric_distance -maxgeodist maximum_geocentric_distance -geologstep logarithmic_step_size_for_geocentric_distance_bins -mingeoobs min_geocentric_dist_at_observation(AU) -min_obstanvel min_tanvel_wrt_observer -minimpactpar min_impact_parameter(km) -tanveltol tanveltol -veltol_changerad veltol_changerad -vinf max_v_inf -nightstep nightstep -use_univar 1_for_universal_variables -use_rr 1_for_heliolincRR -use_dbscan 1_for_DBSCAN -use_taylor 1_for_TaylorSeries_hypotheses -use_smallclust 1_for_lowmem_clusters -use_smalltrk 1_for_lowmem_tracklets -use_streaming 1_for_output_streaming -outsum summary_file -clust2det clust2detfile -verbose verbosity\n";
  cerr << "\nor, at minimum:\n\n";
  cerr << "heliolinc_TNV -imgs imfile -pairdets paired detection file -tracklets tracklet file -trk2det tracklet-to-detection file -earth earthfile -heliodist heliocentric_dist_vel_acc_file\n";
  cerr << "\nNote that the minimum invocation leaves some things set to defaults\n";
  cerr << "that you may well wish to specify: in particular, the output file names\n";
  
}
    
int main(int argc, char *argv[])
{
  vector <hldet> detvec = {};
  vector <hldet> pairdets = {};
  vector <hlimage> image_log;
  vector <hlradhyp> radhyp;
  vector <EarthState> earthpos;
  HeliolincConfig config;
  vector <hlclust> outclust;
  vector <longpair> clust2det;
  string imfile,pairdetfile,trackletfile,trk2detfile,earthfile,accelfile;
  string sumfile = "sumfile_test.csv";
  string clust2detfile = "clust2detfile_test.csv";
  int default_clustrad, default_clustchangerad, default_npt, default_minobsnights;
  int default_mintimespan, default_mingeodist, default_maxgeodist;
  int default_geologstep,default_clust2detfile,default_sumfile;
  int default_mingeoobs, default_minimpactpar;
  int default_use_univar, default_max_v_inf;
  int default_tanveltol, default_veltol_changerad, default_min_obstanvel;
  default_clustrad = default_clustchangerad = default_npt = default_minobsnights = 1;
  default_mintimespan = default_min_obstanvel = 1;
  default_mingeodist = default_maxgeodist = default_geologstep = 1;
  default_clust2detfile = default_sumfile = 1;
  default_mingeoobs = default_minimpactpar = 1;
  default_use_univar = default_max_v_inf = 1;
  default_tanveltol = default_veltol_changerad = 1;
  int default_use_dbscan,default_use_taylor,default_use_rr;
  default_use_dbscan = default_use_taylor = default_use_rr = 1;
  int default_use_smallclust,default_use_smalltrk,default_use_streaming;
  default_use_smallclust = default_use_smalltrk = default_use_streaming = 1;
  int default_nightstep = 1;
  ofstream outstream1;
  long i=0;
  long clustct=0;
  int status=0;
  int use_sad=0;
  
  i=1;
  while(i<argc) {
    cout << "Checking out argv[" << i << "] = " << argv[i] << ".\n";
    if(string(argv[i]) == "-i" || string(argv[i]) == "-imgs" || string(argv[i]) == "-inimgs" || string(argv[i]) == "-img" || string(argv[i]) == "--inimgs" || string(argv[i]) == "--img" || string(argv[i]) == "--image" || string(argv[i]) == "--images") {
      if(i+1 < argc) {
	//There is still something to read;
	imfile=argv[++i];
	i++;
      }
      else {
	cerr << "Image file keyword supplied with no corresponding argument\n";
	show_usage();
	return(1);
      }
    } else if(string(argv[i]) == "-pairdet" || string(argv[i]) == "-pairdets" || string(argv[i]) == "-detpairs") {
      if(i+1 < argc) {
	//There is still something to read;
        pairdetfile=argv[++i];
	i++;
      }
      else {
	cerr << "Input paired detection file keyword supplied with no corresponding argument\n";
	show_usage();
	return(1);
      }
    }  else if(string(argv[i]) == "-tracklets" || string(argv[i]) == "-tracklet" || string(argv[i]) == "-tf" || string(argv[i]) == "-trkfile" || string(argv[i]) == "-trackletfile" || string(argv[i]) == "--trackletfile") {
      if(i+1 < argc) {
	//There is still something to read;
	trackletfile=argv[++i];
	i++;
      }
      else {
	cerr << "Input tracklet file keyword supplied with no corresponding argument\n";
	show_usage();
	return(1);
      }
    } else if(string(argv[i]) == "-trk2det" || string(argv[i]) == "-t2d" || string(argv[i]) == "-tracklet2detection" || string(argv[i]) == "-trk2detfile" || string(argv[i]) == "--tracklet2detection") {
      if(i+1 < argc) {
	//There is still something to read;
	trk2detfile=argv[++i];
	i++;
      }
      else {
	cerr << "Input trk2det file keyword supplied with no corresponding argument\n";
	show_usage();
	return(1);
      }
    } else if(string(argv[i]) == "-autorun" || string(argv[i]) == "-auto" || string(argv[i ]) == "-automjd" || string(argv[i]) == "-autoMJD" || string(argv[i]) == "--auto_mjd" || string(argv[i]) == "--autoMJD" || string(argv[i]) == "--autorun" || string(argv[i]) == "--automjd") {
      if(i+1 < argc) {
	//There is still something to read;
	config.autorun=stoi(argv[++i]);
	i++;
      }
      else {
	cerr << "Autorun keyword supplied with no corresponding argument";
	show_usage();
	return(1);
      }
    } else if(string(argv[i]) == "-m" || string(argv[i]) == "-mjd" || string(argv[i ]) == "-mjdref" || string(argv[i]) == "-MJD" || string(argv[i]) == "--mjd" || string(argv[i]) == "--MJD" || string(argv[i]) == "--modifiedjulianday" || string(argv[i]) == "--ModifiedJulianDay") {
      if(i+1 < argc) {
	//There is still something to read;
	config.MJDref=stod(argv[++i]);
	i++;
      }
      else {
	cerr << "Reference MJD keyword supplied with no corresponding argument";
	show_usage();
	return(1);
      }
    } else if(string(argv[i]) == "-earth" || string(argv[i]) == "-e" || string(argv[i]) == "-Earth" || string(argv[i]) == "--earthfile" || string(argv[i]) == "--Earthfile" || string(argv[i]) == "--earth" || string(argv[i]) == "--Earth" || string(argv[i]) == "-obspos" || string(argv[i]) == "-op" || string(argv[i]) == "-obsvec" || string(argv[i]) == "--observer" || string(argv[i]) == "--observer_position" || string(argv[i]) == "--observer_statevec") {
      if(i+1 < argc) {
	//There is still something to read;
	earthfile=argv[++i];
	i++;
      }
      else {
	cerr << "Earth file keyword supplied with no corresponding argument\n";
	show_usage();
	return(1);
      }
    } else if(string(argv[i]) == "-heliodist" || string(argv[i]) == "-hd" || string(argv[i]) == "-heliodva" || string(argv[i]) == "-hdva" || string(argv[i]) == "--heliodistvelacc" || string(argv[i]) == "--heliodva") {
      if(i+1 < argc) {
	//There is still something to read;
	accelfile=argv[++i];
	i++;
      }
      else {
	cerr << "Heliocentric distance, velocity, and acceleration\nfile keyword supplied with no corresponding argument\n";
	show_usage();
	return(1);
      }
    } else if(string(argv[i]) == "-clustrad" || string(argv[i]) == "-cr" || string(argv[i]) == "-crad" || string(argv[i]) == "-cluster" || string(argv[i]) == "--clustrad" || string(argv[i]) == "--clusterradius" || string(argv[i]) == "--clusterrad") {
      if(i+1 < argc) {
	//There is still something to read;
	config.clustrad=stod(argv[++i]);
	default_clustrad = 0;
	i++;
      }
      else {
	cerr << "Clustering radius keyword supplied with no corresponding argument\n";
	show_usage();
	return(1);
      }
    } else if(string(argv[i]) == "-clustchangerad" || string(argv[i]) == "-clustchangerad" || string(argv[i]) == "-clustchangerad") {
      if(i+1 < argc) {
	//There is still something to read;
	config.clustchangerad = stod(argv[++i]); 
	default_clustchangerad = 0;
	i++;
      }
      else {
	cerr << "Transition distance for cluster scaling keyword supplied with no corresponding argument\n";
	show_usage();
	return(1);
      }
    } else if(string(argv[i]) == "-cluster_minpts" || string(argv[i]) == "-minpts" || string(argv[i]) == "-cluster_min_points" || string(argv[i]) == "-min_cluster_pts" || string(argv[i]) == "-minclusterpts" || string(argv[i]) == "-min_clusterpts" || string(argv[i]) == "-npt" || string(argv[i]) == "-npoints" || string(argv[i]) == "-minpts" || string(argv[i]) == "-np" || string(argv[i]) == "--npt" || string(argv[i]) == "--dbscan_npt" || string(argv[i]) == "--DBSCANnpt") {
      if(i+1 < argc) {
	//There is still something to read;
	config.cluster_minpts=stoi(argv[++i]);
	default_npt = 0;
	i++;
      }
      else {
	cerr << "Cluster minimum number of points keyword supplied with no corresponding argument\n";
	show_usage();
	return(1);
      }
    } else if(string(argv[i]) == "-minobsnights" || string(argv[i]) == "-mon" || string(argv[i]) == "-minobs" || string(argv[i]) == "-minobsnight" || string(argv[i]) == "--minobsnights" || string(argv[i]) == "--minobsnight" || string(argv[i]) == "-minobsn") {
      if(i+1 < argc) {
	//There is still something to read;
	config.minobsnights=stoi(argv[++i]);
	default_minobsnights = 0;
	i++;
      }
      else {
	cerr << "Min. number of observing nights keyword supplied with no corresponding argument\n";
	show_usage();
	return(1);
      }
    } else if(string(argv[i]) == "-mintimespan" || string(argv[i]) == "-mts" || string(argv[i]) == "-minspan" || string(argv[i]) == "-mintspan" || string(argv[i]) == "--mts" || string(argv[i]) == "--mintimespan" || string(argv[i]) == "--mintspan") {
      if(i+1 < argc) {
	//There is still something to read;
	config.mintimespan=stod(argv[++i]);
	default_mintimespan = 0;
	i++;
      }
      else {
	cerr << "Minimum time span keyword supplied with no corresponding argument\n";
	show_usage();
	return(1);
      }
    }  else if(string(argv[i]) == "-mingeodist" || string(argv[i]) == "-mingd" || string(argv[i]) == "-mingeo" || string(argv[i]) == "-mingeod" || string(argv[i]) == "--mingeodist" || string(argv[i]) == "--minimum_geocentr_dist") {
      if(i+1 < argc) {
	//There is still something to read;
	config.mingeodist=stod(argv[++i]);
	default_mingeodist = 0;
	i++;
      }
      else {
	cerr << "Minimum geocentric distance keyword supplied with no corresponding argument\n";
	show_usage();
	return(1);
      }
    }  else if(string(argv[i]) == "-maxgeodist" || string(argv[i]) == "-maxgd" || string(argv[i]) == "-maxgeo" || string(argv[i]) == "-maxgeod" || string(argv[i]) == "--maxgeodist" || string(argv[i]) == "--maximum_geocentr_dist") {
      if(i+1 < argc) {
	//There is still something to read;
	config.maxgeodist=stod(argv[++i]);
	default_maxgeodist = 0;
	i++;
      }
      else {
	cerr << "Maximum geocentric distance keyword supplied with no corresponding argument\n";
	show_usage();
	return(1);
      }
    } else if(string(argv[i]) == "-geologstep" || string(argv[i]) == "-gls" || string(argv[i]) == "-geostep" || string(argv[i]) == "-glogstep" || string(argv[i]) == "--geologstep" || string(argv[i]) == "--geodistlogstep" || string(argv[i]) == "--geodiststep") {
      if(i+1 < argc) {
	//There is still something to read;
	config.geologstep=stod(argv[++i]);
	default_geologstep = 0;
	i++;
      }
      else {
	cerr << "Geocentric distance logarithmic step keyword supplied with no corresponding argument\n";
	show_usage();
	return(1);
      }
    } else if(string(argv[i]) == "-mingeoobs" || string(argv[i]) == "-mgo" || string(argv[i]) == "-minobsdist" || string(argv[i]) == "-mindistobs" || string(argv[i]) == "--min_geocentric_obsdist" || string(argv[i]) == "--min_observation_distance" || string(argv[i]) == "--mingeoobs") {
      if(i+1 < argc) {
	//There is still something to read;
	config.mingeoobs=stod(argv[++i]);
	default_mingeoobs = 0;
	i++;
      }
      else {
	cerr << "Minimum geocentric distance at observation keyword supplied with no corresponding argument\n";
	show_usage();
	return(1);
      }
    } else if(string(argv[i]) == "-min_obstanvel" || string(argv[i]) == "-minobstanvel" || string(argv[i]) == "-minotv" || string(argv[i]) == "-obstanvel" || string(argv[i]) == "-otv" || string(argv[i]) == "-motv" || string(argv[i]) == "--min_obs_tanvel") {
      if(i+1 < argc) {
	//There is still something to read;
	config.min_obstanvel=stod(argv[++i]);
	default_min_obstanvel = 0;
	i++;
      }
      else {
	cerr << "Keyword for min. tangential velocity w.r.t. observer supplied with no corresponding argument\n";
	show_usage();
	return(1);
      }
    } else if(string(argv[i]) == "-minimpactpar" || string(argv[i]) == "-mip" || string(argv[i]) == "-minimp" || string(argv[i]) == "-minimppar" || string(argv[i]) == "--minimum_impact_parameter" || string(argv[i]) == "--minimpactpar" || string(argv[i]) == "--min_impact_par") {
      if(i+1 < argc) {
	//There is still something to read;
	config.minimpactpar=stod(argv[++i]);
	default_minimpactpar = 0;
	i++;
      }
      else {
	cerr << "Minimum impact parameter keyword supplied with no corresponding argument\n";
	show_usage();
	return(1);
      }
    } else if(string(argv[i]) == "-tanveltol" || string(argv[i]) == "-veltol" || string(argv[i]) == "-tanvel") {
      if(i+1 < argc) {
	//There is still something to read;
	config.tanveltol=stod(argv[++i]);
	default_tanveltol = 0;
	i++;
      }
      else {
	cerr << "Tangential velocity tolerance keyword supplied with no corresponding argument\n";
	show_usage();
	return(1);
      }
    } else if(string(argv[i]) == "-veltol_changerad" || string(argv[i]) == "-vel_changerad" || string(argv[i]) == "-veltol_changerad") {
      if(i+1 < argc) {
	//There is still something to read;
	config.veltol_changerad=stod(argv[++i]);
	default_veltol_changerad = 0;
	i++;
      }
      else {
	cerr << "Velocity tolerance changerad keyword supplied with no corresponding argument\n";
	show_usage();
	return(1);
      }
    } else if(string(argv[i]) == "-vinf" || string(argv[i]) == "-maxvinf" || string(argv[i]) == "-max_v_inf") {
      if(i+1 < argc) {
	//There is still something to read;
	config.max_v_inf=stod(argv[++i]);
	default_max_v_inf = 0;
	i++;
      }
      else {
	cerr << "v_infinity keyword supplied with no corresponding argument\n";
	show_usage();
	return(1);
      }
    } else if(string(argv[i]) == "-nightstep" || string(argv[i]) == "-night_step" || string(argv[i]) == "-trackdiv") {
      if(i+1 < argc) {
	//There is still something to read;
	config.nightstep=stod(argv[++i]);
	default_nightstep = 0;
	i++;
      }
      else {
	cerr << "Nighstep keyword supplied with no corresponding argument\n";
	show_usage();
	return(1);
      }
    } else if(string(argv[i]) == "-use_univar" || string(argv[i]) == "-useunivar" || string(argv[i]) == "-univar" || string(argv[i]) == "-universalvar") {
      if(i+1 < argc) {
	//There is still something to read;
	config.use_univar=stoi(argv[++i]);
	default_use_univar = 0;
	i++;
      }
      else {
	cerr << "Minimum impact parameter keyword supplied with no corresponding argument\n";
	show_usage();
	return(1);
      }
    }  else if(string(argv[i]) == "-use_rr" || string(argv[i]) == "-RR" || string(argv[i]) == "-rr" || string(argv[i]) == "-userr" || string(argv[i]) == "-use_RR") {
      if(i+1 < argc) {
	//There is still something to read;
	config.use_rr=stoi(argv[++i]);
	default_use_rr = 0;
	i++;
      }
      else {
	cerr << "use_rr keyword supplied with no corresponding argument\n";
	show_usage();
	return(1);
      }
    } else if(string(argv[i]) == "-use_dbscan" || string(argv[i]) == "-dbscan" || string(argv[i]) == "-DBSCAN" || string(argv[i]) == "-DBscan") {
      if(i+1 < argc) {
	//There is still something to read;
	config.use_dbscan=stoi(argv[++i]);
	default_use_dbscan = 0;
	i++;
      }
      else {
	cerr << "use_dbscan keyword supplied with no corresponding argument\n";
	show_usage();
	return(1);
      }
    } else if(string(argv[i]) == "-use_taylor" || string(argv[i]) == "-taylor" || string(argv[i]) == "-Taylor" || string(argv[i]) == "-usetaylor") {
      if(i+1 < argc) {
	//There is still something to read;
	config.use_taylor=stoi(argv[++i]);
	default_use_taylor = 0;
	i++;
      }
      else {
	cerr << "use_taylor keyword supplied with no corresponding argument\n";
	show_usage();
	return(1);
      }
    } else if(string(argv[i]) == "-use_smallclust" || string(argv[i]) == "-smallclust" || string(argv[i]) == "-use_lowmem" || string(argv[i]) == "-lowmem") {
      if(i+1 < argc) {
	//There is still something to read;
	config.use_smallclust=stoi(argv[++i]);
	default_use_smallclust = 0;
	i++;
      }
      else {
	cerr << "use_smallclust keyword supplied with no corresponding argument\n";
	show_usage();
	return(1);
      }
    } else if(string(argv[i]) == "-use_smalltrk" || string(argv[i]) == "-smalltrk" || string(argv[i]) == "-use_uint" || string(argv[i]) == "-uint") {
      if(i+1 < argc) {
	//There is still something to read;
	config.use_smalltrk=stoi(argv[++i]);
	default_use_smalltrk = 0;
	i++;
      }
      else {
	cerr << "use_smalltrk keyword supplied with no corresponding argument\n";
	show_usage();
	return(1);
      }
    } else if(string(argv[i]) == "-use_streaming" || string(argv[i]) == "-streaming" || string(argv[i]) == "-usestreaming" || string(argv[i]) == "-Streaming") {
      if(i+1 < argc) {
	//There is still something to read;
	config.use_streaming=stoi(argv[++i]);
	default_use_streaming = 0;
	i++;
      }
      else {
	cerr << "use_streaming keyword supplied with no corresponding argument\n";
	show_usage();
	return(1);
      }
    } else if(string(argv[i]) == "-outsum" || string(argv[i]) == "-sum" || string(argv[i]) == "-rms" || string(argv[i]) == "-outrms" || string(argv[i]) == "-or" || string(argv[i]) == "--outsummaryfile" || string(argv[i]) == "--outsum" || string(argv[i]) == "--sum") {
      if(i+1 < argc) {
	//There is still something to read;
	sumfile=argv[++i];
	default_sumfile = 0;
	i++;
      }
      else {
	cerr << "Output summary file keyword supplied with no corresponding argument\n";
	show_usage();
	return(1);
      }
    } else if(string(argv[i]) == "-clust2det" || string(argv[i]) == "-c2d" || string(argv[i]) == "-clust2detfile" || string(argv[i]) == "--clust2detfile" || string(argv[i]) == "--clust2det" || string(argv[i]) == "--cluster_to_detection") {
      if(i+1 < argc) {
	//There is still something to read;
	clust2detfile=argv[++i];
	default_clust2detfile = 0;
	i++;
      }
      else {
	cerr << "Output cluster-to-detection file keyword supplied with no corresponding argument\n";
	show_usage();
	return(1);
      }
    } else if(string(argv[i]) == "-verbose" || string(argv[i]) == "-verb" || string(argv[i]) == "-VERBOSE" || string(argv[i]) == "-VERB" || string(argv[i]) == "--verbose" || string(argv[i]) == "--VERBOSE" || string(argv[i]) == "--VERB") {
      if(i+1 < argc) {
	//There is still something to read;
	config.verbose=stoi(argv[++i]);
	i++;
      }
      else {
	cerr << "Verbosity keyword supplied with no corresponding argument\n";
	show_usage();
	return(1);
      }
    } else {
      cerr << "ERROR: unrecognized keyword or argument " << argv[i] << "\n";
      return(2);
    }
  }
  if(config.minobsnights > config.cluster_minpts) config.minobsnights = config.cluster_minpts; // Otherwise the low setting of cluster_minpts is not operative.

  if(argc<11)
    {
      cerr << "Too few arguments even for minimalist invocation:\n";
      show_usage();
      return(1);
    }
  
  cout.precision(17);  
  cout << "input image file " << imfile << "\n";
  cout << "input detection file " << pairdetfile << "\n";
  cout << "input tracklet file " << trackletfile << "\n";
  cout << "input trk2det file " << trk2detfile << "\n";
  cout << "input Earth position file " << earthfile << "\n";
  cout << "input heliocentric hypothesis file " << accelfile << "\n";
  cout << "input reference MJD " << config.MJDref << "\n";

  // Catch required parameters if missing
  if(imfile.size()<=0) {
    cout << "\nERROR: input image file is required\n";
    show_usage();
    return(1);
  } else if(pairdetfile.size()<=0) {
    cout << "\nERROR: input detection file is required\n";
    show_usage();
    return(1);
  } else if(trackletfile.size()<=0) {
    cout << "\nERROR: input tracklet file is required\n";
    show_usage();
    return(1);
  } else if(trk2detfile.size()<=0) {
    cout << "\nERROR: input trk2det file is required\n";
    show_usage();
    return(1);
  } else if(earthfile.size()<=0) {
    cout << "\nERROR: input Earth position file is required:\n";
    cout << "e.g. Earth1day2020s_02a.txt\n";
    show_usage();
    return(1);
  } else if(accelfile.size()<=0) {
    cout << "\nERROR: input heliocentric hypothesis file is required\n";
    show_usage();
    return(1);
  }

  if(config.use_streaming>0) {
    cerr << "ERROR: use_streaming is set, but this version of heliolinc\n";
    cerr << "does not support it. Use a heliolinc_omp or a similar multi-threaded version\n";
    return(1);
  }

  // Catch case where max v_inf > 0 but universal variables are not set.
  if(config.max_v_inf>0.0l && config.use_univar<=0) {
    cerr << "ERROR: unbound orbits being probed (max_v_inf = " << config.max_v_inf << "),\n";
    cerr << "but code is not set to use universal variables (use_univar = " << config.use_univar << ").\n";
    show_usage();
    return(1);
  }

  if(config.tanveltol>0.0) use_sad=1;
  if(use_sad>0 && config.use_taylor>0) {
    cerr << "ERROR: non-Keplerian (Taylor Series) interpretation of hypotheses\n";
    cerr << "requested at the same time as intrinsically Keplerian tangential velocity\n";
    cerr << "checking with tolerance " << config.tanveltol << "\n";
    return(0);
  }
   if(use_sad==0 && config.min_obstanvel>0) {
     cerr << "ERROR: min_obstanvel set to operative value (" << config.min_obstanvel << " km/sec without also\n";
     cerr << "setting the heliocentric tangential velocity tolerance tanveltol, which\n";
     cerr << "remains at its default, non-operative value of " << config.tanveltol << " km/sec.\n";
     cerr << "Re-run with -tanveltol set. If you want tanveltol not to exclude any,\n";
     cerr << "tracklets, set it to a very large positive value (e.g. -tanveltol 1000.0).\n";
     cerr << "This, combined with your desired operative setting for min_obstanvel, will render\n";
     cerr << "min_obstanvel operative while tanveltol is not.\n";
     return(0);
   }
 
  if(config.autorun==0) cout << "Automatic reference MJD is disabled: user set MJDref = " << config.MJDref << "\n";
  if(default_clustrad==1) cout << "Defaulting to cluster radius = " << config.clustrad << "km\n";
  else cout << "input clustering radius " << config.clustrad << "km\n";
  if(default_clustchangerad==1) cout << "Defaulting to min. geocentric distance for cluster scaling = " << config.clustchangerad << "AU\n";
  else cout << "Min. geocentric distance for cluster scaling is " << config.clustchangerad << "AU\n";
  cout << "Minimum cluster radius, which will apply for all geocentric distances less\n";
  cout << "than " << config.clustchangerad << "AU, is " << config.clustrad*config.clustchangerad/REF_GEODIST << "km\n";
  if(default_npt==1) cout << "Defaulting to cluster_minpts (min. no. of tracklets in a linkage) = " << config.cluster_minpts << "\n";
  else cout << "input cluster_minpts (min. no. of tracklets in a linkage) is " << config.cluster_minpts << "\n";
  if(default_minobsnights==1) cout << "Defaulting to minimum number of unique nights = " << config.minobsnights << "\n";
  else cout << "minimum number of unique nights is " << config.minobsnights << "\n";
  if(default_mintimespan==1) cout << "Defaulting to minimum time span for a linkage = " << config.mintimespan << " days\n";
  else cout << "minimum time span for a linkage is " << config.mintimespan << " days\n";
  if(default_mingeodist==1) cout << "Defaulting to minimum geocentric distance = " << config.mingeodist << " AU\n";
  else cout << "minimum geocentric distance is " << config.mingeodist << " AU\n";
  if(default_maxgeodist==1) cout << "Defaulting to maximum geocentric distance = " << config.maxgeodist << " AU\n";
  else cout << "maximum geocentric distance is " << config.maxgeodist << " AU\n";
  if(default_geologstep==1) cout << "Defaulting to logarithmic step size for geocentric distance bins = " << config.geologstep << "\n";
  else cout << "logarithmic step size for geocentric distance bins is " << config.geologstep << "\n";
  if(default_mingeoobs==1) cout << "Defaulting to minimum geocentric distance at observation = " << config.mingeoobs << " AU\n";
  else cout << "Minimum geocentric distance at observation = " << config.mingeoobs << " AU\n";
  if(default_minimpactpar==1) cout << "Defaulting to minimum impact parameter = " << config.minimpactpar << " km\n";
  else cout << "Minimum impact parameter is " << config.minimpactpar << " km\n";
  if(default_use_univar==1) cout << "For Keplerian integration, defaulting to f and g functions\nrather than universal variables\n";
  else if(config.use_univar>0) cout << "Using universal variables for Keplerian integration\n";
  else cout << "Using f and g functions for Keplerian integration\n";
  if(default_use_dbscan==1) cout << "Defaulting to k-d tree range query rather than DBSCAN for clustering\n";
  else if(config.use_dbscan>0) cout << "Using DBSCAN for clustering rather than a k-d tree range query.\n";
  else cout << "Using a k-d tree range query rather than DBSCAN for clustering\n";
  if(default_use_rr==1) cout << "Defaulting to using heliolinc3D, r-v parameter space for clustering\n";
  else if(config.use_rr>0) cout << "Using heliolincRR parameter space rather than r-v for clustering.\n";
  else cout << "Using heliolinc3D, r-v parameter space for clustering\n";
  if(default_use_taylor==1) cout << "Defaulting to Keplerian interpretation of heliolinc hypotheses.\n";
  else if(config.use_taylor>0) cout << "Using Taylor Series interpretation of heliolinc hypotheses.\n";
  else cout << "Using Keplerian interpretation of heliolinc hypotheses.\n";
  if(default_use_taylor==1) cout << "Defaulting to Keplerian interpretation of heliolinc hypotheses.\n";
  else if(config.use_taylor>0) cout << "Using Taylor Series interpretation of heliolinc hypotheses.\n";
  else cout << "Using Keplerian interpretation of heliolinc hypotheses.\n";
  if(default_use_smallclust==1) cout << "Defaulting to using normal, full-precision cluster class.\n";
  else if(config.use_smallclust>0) cout << "Using using compressed, memory-efficient cluster class.\n";
  else cout << "Using using normal, full-precision cluster class.\n";
  if(default_use_smalltrk==1) cout << "Defaulting to using normal, full-precision tracklet classes.\n";
  else if(config.use_smalltrk>0) cout << "Using using compressed, memory-efficient tracklet classes.\n";
  else cout << "Using using normal, full-precision tracklet classes.\n";
  if(default_max_v_inf==1) cout << "Defaulting to maximum v_inf relative to the sun = " << config.max_v_inf << " km/sec\n";
  else cout << "Maximum v_inf relative to the sun is " << config.max_v_inf << " km\n";
  if(default_nightstep==1) cout << "Defaulting to nightstep = " << config.nightstep << " days\n";
  else cout << "Nightstep is " << config.nightstep << " days\n";
  if(default_tanveltol==1) cout << "Defaulting to tangential velocity tolerance = " << config.tanveltol << " km/sec\n";
  else cout << "Tangential velocity tolerance is " << config.tanveltol << " km\n";
  if(default_veltol_changerad==1) cout << "Defaulting to veltol_changerad = " << config.veltol_changerad << " AU\n";
  else cout << "veltol_changerad is " << config.veltol_changerad << " km\n";
  if(default_sumfile==1) cout << "WARNING: using default name " << sumfile << " for summary output file\n";
  else cout << "summary output file " << sumfile << "\n";
  if(default_clust2detfile==1) cout << "WARNING: using default name " << clust2detfile << " for output clust2det file\n";
  else cout << "output clust2det file " << clust2detfile << "\n";

  cout << "Heliocentric ephemeris for Earth is named " << earthfile << "\n";
  status = read_horizons_csv(earthfile, earthpos);
  if(status!=0) {
    cerr << "ERROR: could not successfully Earth ephemeris file " << earthfile << "\n";
    cerr << "read_horizons_csv returned status = " << status << ".\n";
   return(1);
  } 

  cout << "File with heliocentric motion hypothesis is named " << accelfile << "\n";
  status=read_radhyp_file(accelfile, radhyp, config.verbose);
  if(status!=0) {
    cerr << "ERROR: could not successfully read accleration file " << accelfile << "\n";
    cerr << "read_radhyp_file returned status = " << status << ".\n";
   return(1);
  }
  
  detvec={};
  int headerlines = 1;
  status=read_hldet_file(pairdetfile, detvec, headerlines, config.verbose);
  if(status!=0) {
    cerr << "ERROR: could not successfully read paired detection file " << pairdetfile << "\n";
    cerr << "read_hldet_file returned status = " << status << ".\n";
   return(1);
  }
  cout << "Read " << detvec.size() << " data lines from paired detection file " << pairdetfile << "\n";
  if(config.use_smalltrk==1 && detvec.size()>=UINT_MAX) {
    cerr << "ERROR: too many detections to work with -use_smalltrk 1\n";
    cerr << "(" << detvec.size() << " >= " << UINT_MAX << "). Re-run with -use_smalltrk 0 (which is the default).\n";
    return(1);
  }

  // Find MJD for first and last detections.
  double minMJD = detvec[0].MJD;
  double maxMJD = detvec[0].MJD;
  for(i=0; i<long(detvec.size()); i++) {
    if(minMJD > detvec[i].MJD) minMJD = detvec[i].MJD;
    if(maxMJD < detvec[i].MJD) maxMJD = detvec[i].MJD;
  }
  // Do something useful in the specific case that there is
  // an input data file, but no reference MJD, and config.autorun==0,
  // prohibiting the code from automatically generating a suitable
  // reference time.
  if(!isnormal(config.MJDref) || config.MJDref < minMJD || config.MJDref > maxMJD) {
    if(config.autorun<=0) {
      show_usage();
      cout << "\nERROR: input positive-valued reference MJD is required\n";
      cout << "MJD range is " << minMJD << " to " << maxMJD << "\n";
      cout << fixed << setprecision(2) << "Suggested reference value is " << minMJD*0.5L + maxMJD*0.5L << "\n";
      cout << "based on your input detection catalog " << pairdetfile << "\n";
      cout << "see usage information above\n";
      return(1);
    } 
  }
  
  image_log={};
  status=read_image_file2(imfile, image_log);
  if(status!=0) {
    cerr << "ERROR: could not successfully read image file " << imfile << "\n";
    cerr << "read_image_file2 returned status = " << status << ".\n";
   return(1);
  }
  cout << "Read " << image_log.size() << " data lines from image file " << imfile << "\n";

  if(config.use_smalltrk==1) {
    vector <uint_tracklet> tracklets;
    vector <uint_pair> trk2det;
    tracklets={};
    trk2det={};
    
    status=read_tracklet_file_uint(trackletfile, tracklets, config.verbose);
    if(status!=0) {
      cerr << "ERROR: could not successfully read tracklet file " << trackletfile << "\n";
      cerr << "read_tracklet_file returned status = " << status << ".\n";
      return(1);
    }
    cout << "Read " << tracklets.size() << " data lines from tracklet file " << trackletfile << "\n";
    if(tracklets.size()>=UINT_MAX) {
      cerr << "ERROR: too many tracklets to work with -use_smalltrk 1\n";
      cerr << "(" << tracklets.size() << " >= " << UINT_MAX << "). Re-run with -use_smalltrk 0 (which is the default).\n";
      return(1);
    }
    status=read_uint_pair_file(trk2detfile, trk2det, config.verbose);
    if(status!=0) {
      cerr << "ERROR: could not successfully read trk2det file " << trk2detfile << "\n";
      cerr << "read_longpair_file returned status = " << status << ".\n";
      return(1);
    }
    cout << "Read " << trk2det.size() << " data lines from trk2det file " << trk2detfile << "\n";
    if(config.use_smallclust==1) {
      status=heliolinc_alg_uint_lowmem(image_log, detvec, tracklets, trk2det, radhyp, earthpos, config, outclust, clust2det);
      if(status!=0) {
	cerr << "ERROR: heliolinc_alg_uint_lowmem failed with status " << status << "\n";
	return(status);
      }
    } else {
      status=heliolinc_alg_uint(image_log, detvec, tracklets, trk2det, radhyp, earthpos, config, outclust, clust2det);
      if(status!=0) {
	cerr << "ERROR: heliolinc_alg_uint failed with status " << status << "\n";
	return(status);
      }
    }
  } else {
    vector <tracklet> tracklets;
    vector <longpair> trk2det;
    tracklets={};
    trk2det={};
    
    status=read_tracklet_file(trackletfile, tracklets, config.verbose);
    if(status!=0) {
      cerr << "ERROR: could not successfully read tracklet file " << trackletfile << "\n";
      cerr << "read_tracklet_file returned status = " << status << ".\n";
      return(1);
    }
    cout << "Read " << tracklets.size() << " data lines from tracklet file " << trackletfile << "\n";
  
    status=read_longpair_file(trk2detfile, trk2det, config.verbose);
    if(status!=0) {
      cerr << "ERROR: could not successfully read trk2det file " << trk2detfile << "\n";
      cerr << "read_longpair_file returned status = " << status << ".\n";
      return(1);
    }
    cout << "Read " << trk2det.size() << " data lines from trk2det file " << trk2detfile << "\n";
    status=heliolinc_alg_all(image_log, detvec, tracklets, trk2det, radhyp, earthpos, config, outclust, clust2det);
    if(status!=0) {
      cerr << "ERROR: heliolinc_alg_all failed with status " << status << "\n";
      return(status);
    }
  }
  
  outstream1.open(sumfile);
  cout << "Writing " << outclust.size() << " lines to output cluster-summary file " << sumfile << "\n";
  outstream1 << "#clusternum,posRMS,velRMS,totRMS,astromRMS,pairnum,timespan,uniquepoints,obsnights,metric,rating,reference_MJD,heliohyp0,heliohyp1,heliohyp2,posX,posY,posZ,velX,velY,velZ,orbit_a,orbit_e,orbit_incl,orbit_MJD,orbitX,orbitY,orbitZ,orbitVX,orbitVY,orbitVZ,orbit_eval_count\n";
  for(clustct=0 ; clustct<long(outclust.size()); clustct++) {
    outstream1 << fixed << setprecision(3) << outclust[clustct].clusternum << "," << outclust[clustct].posRMS << "," << outclust[clustct].velRMS << "," << outclust[clustct].totRMS << ",";
    outstream1 << fixed << setprecision(4) << outclust[clustct].astromRMS << ",";
    outstream1 << fixed << setprecision(6) << outclust[clustct].pairnum << "," << outclust[clustct].timespan << "," << outclust[clustct].uniquepoints << "," << outclust[clustct].obsnights << "," << outclust[clustct].metric << "," << outclust[clustct].rating << ",";
    outstream1 << fixed << setprecision(6) << outclust[clustct].reference_MJD << "," << outclust[clustct].heliohyp0 << "," << outclust[clustct].heliohyp1 << "," << outclust[clustct].heliohyp2 << ",";
    outstream1 << fixed << setprecision(1) << outclust[clustct].posX << "," << outclust[clustct].posY << "," << outclust[clustct].posZ << ",";
    outstream1 << fixed << setprecision(4) << outclust[clustct].velX << "," << outclust[clustct].velY << "," << outclust[clustct].velZ << ",";
    outstream1 << fixed << setprecision(6) << outclust[clustct].orbit_a << "," << outclust[clustct].orbit_e << "," << outclust[clustct].orbit_incl << "," << outclust[clustct].orbit_MJD << ",";
    outstream1 << fixed << setprecision(1) << outclust[clustct].orbitX << "," << outclust[clustct].orbitY << "," << outclust[clustct].orbitZ << ",";
    outstream1 << fixed << setprecision(4) << outclust[clustct].orbitVX << "," << outclust[clustct].orbitVY << "," << outclust[clustct].orbitVZ << "," << outclust[clustct].orbit_eval_count << "\n";
  }
  outstream1.close();
  outstream1.open(clust2detfile);
  cout << "Writing " << clust2det.size() << " lines to output clust2det file " << clust2detfile << "\n";
  outstream1 << "#clusternum,detnum\n";
  for(clustct=0 ; clustct<long(clust2det.size()); clustct++) {
    outstream1 << clust2det[clustct].i1 << "," << clust2det[clustct].i2 << "\n";
  }
  outstream1.close();
  
  return(0);
}
