// December 02, 2026: helio_highgrade:
// Use some of the machinery of heliolinc to attempt to high-grade
// an input data set, for later more efficient processing with the
// full-scale version of heliolinc. This new version avoids the
// inefficient vector_insert operations.

#include "solarsyst_dyn_geo01.h"
#include "cmath"

static void show_usage()
{
  cerr << "Usage: helio_highgrade2 -imgs imfile -pairdets paired detection file -tracklets tracklet file -trk2det tracklet-to-detection file -mjd mjdref -autorun 1=yes_auto-generate_MJDref -obspos observer_position_file -heliodist heliocentric_dist_vel_acc_file -clustrad clustrad -clustchangerad min_distance_for_cluster_scaling -npt dbscan_npt -mintimespan mintimespan -minobs min_unique_obs -mingeodist minimum_geocentric_distance -maxgeodist maximum_geocentric_distance -geologstep logarithmic_step_size_for_geocentric_distance_bins -mingeoobs min_geocentric_dist_at_observation(AU) -minimpactpar min_impact_parameter(km) -useunivar 1_for_univar_0_for_fgfunc -vinf max_v_inf  -outdets output detection file -verbose verbosity\n";
  cerr << "\nor, at minimum:\n\n";
  cerr << "helio_highgrade2 -imgs imfile -pairdets paired detection file -tracklets tracklet file -trk2det tracklet-to-detection file -heliodist heliocentric_dist_vel_acc_file\n\n";
  cerr << "\nNote that the minimum invocation leaves some things set to defaults\n";
  cerr << "that you may well wish to specify: in particular, the output file names\n";
  
}
    
int main(int argc, char *argv[])
{
  vector <hldet> detvec = {};
  vector <hldet> pairdets = {};
  vector <hlimage> image_log;
  vector <tracklet> tracklets;
  vector <longpair> trk2det;
  vector <hlradhyp> radhyp;
  vector <EarthState> earthpos;
  HeliolincConfig config;
  vector <hldet> outdets = {};
  string imfile,pairdetfile,trackletfile,trk2detfile,planetfile,accelfile;
  string outdetsfile = "highgradeout_test.csv";
  int default_clustrad, default_clustchangerad, default_npt;
  int default_mingeodist, default_maxgeodist;
  int default_geologstep,default_outdetsfile;
  int default_mingeoobs, default_minimpactpar;
  int default_use_univar, default_max_v_inf;
  int default_mintimespan, default_minobsnum;
  default_clustrad = default_clustchangerad = default_npt = 1;
  default_mingeodist = default_maxgeodist = default_geologstep = 1;
  default_mintimespan = default_minobsnum = default_outdetsfile = 1;
  default_mingeoobs = default_minimpactpar = 1;
  default_use_univar = default_max_v_inf = 1;
  ofstream outstream1;
  long i=0;
  int status=0;
  config.mintimespan = 1.0;
  long minobsnum = 4;
  
  i=1;
  while(i<argc) {
    cout << "Checking out argv[" << i << "] = " << argv[i] << ".\n";
    if(string(argv[i]) == "-pairdet" || string(argv[i]) == "-pairdets" || string(argv[i]) == "-detpairs") {
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
    }  else if(string(argv[i]) == "-i" || string(argv[i]) == "-imgs" || string(argv[i]) == "-inimgs" || string(argv[i]) == "-img" || string(argv[i]) == "--inimgs" || string(argv[i]) == "--img" || string(argv[i]) == "--image" || string(argv[i]) == "--images") {
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
    } else if(string(argv[i]) == "-tracklets" || string(argv[i]) == "-tracklet" || string(argv[i]) == "-tf" || string(argv[i]) == "-trkfile" || string(argv[i]) == "-trackletfile" || string(argv[i]) == "--trackletfile") {
      if(i+1 < argc) {
	//There is still something to read;
	trackletfile=argv[++i];
	i++;
      }
      else {
	cerr << "Input pair file keyword supplied with no corresponding argument\n";
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
	cerr << "Input pair file keyword supplied with no corresponding argument\n";
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
    } else if(string(argv[i]) == "-obspos" || string(argv[i]) == "-op" || string(argv[i]) == "-obsvec" || string(argv[i]) == "--observer" || string(argv[i]) == "--observer_position" || string(argv[i]) == "--observer_statevec") {
      if(i+1 < argc) {
	//There is still something to read;
	planetfile=argv[++i];
	i++;
      }
      else {
	cerr << "Observer position file keyword supplied with no corresponding argument\n";
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
    } else if(string(argv[i]) == "-npt" || string(argv[i]) == "-npoints" || string(argv[i]) == "-minpts" || string(argv[i]) == "-np" || string(argv[i]) == "--npt" || string(argv[i]) == "--dbscan_npt" || string(argv[i]) == "--DBSCANnpt") {
      if(i+1 < argc) {
	//There is still something to read;
	config.dbscan_npt=stoi(argv[++i]);
	default_npt = 0;
	i++;
      }
      else {
	cerr << "DBSCAN npt keyword supplied with no corresponding argument\n";
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
    } else if(string(argv[i]) == "-minobsnum" || string(argv[i]) == "-minunique" || string(argv[i]) == "-minobs" || string(argv[i]) == "-min_unique_obs" || string(argv[i]) == "--minobs" || string(argv[i]) == "--min_unique_obs" || string(argv[i]) == "--minobsnum") {
      if(i+1 < argc) {
	//There is still something to read;
	minobsnum=stoi(argv[++i]);
	default_minobsnum = 0;
	i++;
      }
      else {
	cerr << "DBSCAN npt keyword supplied with no corresponding argument\n";
	show_usage();
	return(1);
      }
    } else if(string(argv[i]) == "-mingeodist" || string(argv[i]) == "-mingd" || string(argv[i]) == "-mingeo" || string(argv[i]) == "-mingeod" || string(argv[i]) == "--mingeodist" || string(argv[i]) == "--minimum_geocentr_dist") {
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
    }  else if(string(argv[i]) == "-useunivar" || string(argv[i]) == "-use_univar" || string(argv[i]) == "-univar" || string(argv[i]) == "-universalvar") {
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
    }  else if(string(argv[i]) == "-vinf" || string(argv[i]) == "-maxvinf" || string(argv[i]) == "-max_v_inf") {
      if(i+1 < argc) {
	//There is still something to read;
	config.max_v_inf=stod(argv[++i]);
	default_max_v_inf = 0;
	i++;
      }
      else {
	cerr << "Minimum impact parameter keyword supplied with no corresponding argument\n";
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
    } else if(string(argv[i]) == "-outdets" || string(argv[i]) == "-outdetsfile" || string(argv[i]) == "-out" || string(argv[i]) == "-outfile" || string(argv[i]) == "-od" || string(argv[i]) == "--outdetsfile" || string(argv[i]) == "--output_detection_file" || string(argv[i]) == "--outdets") {
      if(i+1 < argc) {
	//There is still something to read;
	outdetsfile=argv[++i];
	default_outdetsfile = 0;
	i++;
      }
      else {
	cerr << "Output summary file keyword supplied with no corresponding argument\n";
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
      cerr << "Warning: unrecognized keyword or argument " << argv[i] << "\n";
      i++;
    }
  }
  if(config.minobsnights > config.dbscan_npt) config.minobsnights = config.dbscan_npt; // Otherwise the low setting of dbscan_npt is not operative.

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
  cout << "input observer position file " << planetfile << "\n";
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
  } else if(planetfile.size()<=0) {
    cout << "\nERROR: input observer position file is required:\n";
    cout << "e.g. Earth1day2020s_02a.txt\n";
    show_usage();
    return(1);
  } else if(accelfile.size()<=0) {
    cout << "\nERROR: input heliocentric hypothesis file is required\n";
    show_usage();
    return(1);
  }

  // Catch case where max v_inf > 0 but universal variables are not set.
  if(config.max_v_inf>0.0l && config.use_univar<=0) {
    cerr << "ERROR: unbound orbits being probed (max_v_inf = " << config.max_v_inf << "),\n";
    cerr << "but code is not set to use universal variables (use_univar = " << config.use_univar << ".\n";
    show_usage();
    return(1);
  }

  if(default_clustrad==1) cout << "Defaulting to cluster radius = " << config.clustrad << "km\n";
  else cout << "input clustering radius " << config.clustrad << "km\n";
  if(default_clustchangerad==1) cout << "Defaulting to min. geocentric distance for cluster scaling = " << config.clustchangerad << "AU\n";
  else cout << "Min. geocentric distance for cluster scaling is " << config.clustchangerad << "AU\n";
  cout << "Minimum cluster radius, which will apply for all geocentric distances less\n";
  cout << "than " << config.clustchangerad << "AU, is " << config.clustrad*config.clustchangerad/REF_GEODIST << "km\n";
  if(default_npt==1) cout << "Defaulting to DBSCAN npt (min. no. of tracklets in a linkage) = " << config.dbscan_npt << "\n";
  else cout << "input DBSCAN npt (min. no. of tracklets in a linkage) is " << config.dbscan_npt << "\n";
  if(default_minobsnum==1) cout << "Defaulting to min obsnum (min. no. of unique observations in a linkage) = " << minobsnum << "\n";
  else cout << "input min obsnum (min. no. of unique observations in a linkage) is " << minobsnum << "\n";
  if(default_mintimespan==1) cout << "Defaulting to mintimespan = " << config.mintimespan << " days\n";
  else cout << "input mintimespan is " << config.mintimespan << " days\n";
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
  if(default_max_v_inf==1) cout << "Defaulting to maximum v_inf relative to the sun = " << config.max_v_inf << " km/sec\n";
  else cout << "Maximum v_inf relative to the sun is " << config.max_v_inf << " km\n";
  if(default_outdetsfile==1) cout << "WARNING: using default name " << outdetsfile << " for output detection file\n";
  else cout << "output detection file " << outdetsfile << "\n";

  cout << "Heliocentric ephemeris for Earth is named " << planetfile << "\n";
  status = read_horizons_csv(planetfile, earthpos);
  if(status!=0) {
    cerr << "ERROR: could not successfully Earth ephemeris file " << planetfile << "\n";
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
  status=read_pairdet_file(pairdetfile, detvec, config.verbose);
  if(status!=0) {
    cerr << "ERROR: could not successfully read paired detection file " << pairdetfile << "\n";
    cerr << "read_pairdet_file returned status = " << status << ".\n";
   return(1);
  }
  cout << "Read " << detvec.size() << " data lines from paired detection file " << pairdetfile << "\n";

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
  
  tracklets={};
  status=read_tracklet_file(trackletfile, tracklets, config.verbose);
  if(status!=0) {
    cerr << "ERROR: could not successfully read tracklet file " << trackletfile << "\n";
    cerr << "read_tracklet_file returned status = " << status << ".\n";
   return(1);
  }
  cout << "Read " << tracklets.size() << " data lines from tracklet file " << trackletfile << "\n";
  
  trk2det={};
  status=read_longpair_file(trk2detfile, trk2det, config.verbose);
  if(status!=0) {
    cerr << "ERROR: could not successfully read trk2det file " << trk2detfile << "\n";
    cerr << "read_longpair_file returned status = " << status << ".\n";
   return(1);
  }
  cout << "Read " << trk2det.size() << " data lines from trk2det file " << trk2detfile << "\n";
  status=heliolinc_highgrade2(image_log, detvec, tracklets, trk2det, radhyp, earthpos, config, minobsnum, outdets);
  if(status!=0) {
    cerr << "ERROR: heliolinc_highgrade failed with status " << status << "\n";
    return(status);
  } 
  
  outstream1.open(outdetsfile);
  cout << "Writing " << outdets.size() << " lines to output high-graded detection file " << outdetsfile << "\n";
  outstream1 << "#MJD,RA,Dec,mag,trail_len,trail_PA,sigmag,sig_across,sig_along,image,idstring,band,obscode,known_obj,det_qual,origindex\n";
  for(i=0;i<long(outdets.size());i++) {
    outstream1 << fixed << setprecision(7) << outdets[i].MJD << "," << outdets[i].RA << "," << outdets[i].Dec << ",";
    outstream1 << fixed << setprecision(4) << outdets[i].mag << ",";
    outstream1 << fixed << setprecision(2) << outdets[i].trail_len << "," << outdets[i].trail_PA << ",";
    outstream1 << fixed << setprecision(4) << outdets[i].sigmag << ",";
    outstream1 << fixed << setprecision(3) << outdets[i].sig_across << "," << outdets[i].sig_along << ",";
    outstream1 << outdets[i].image << "," << outdets[i].idstring << "," << outdets[i].band << ",";
    outstream1 << outdets[i].obscode << "," << outdets[i].known_obj << ","; 
    outstream1 << outdets[i].det_qual << "," << outdets[i].index << "\n"; 
  }
  outstream1.close();
  
  return(0);
}
