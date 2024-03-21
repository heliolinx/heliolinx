/// March 13, 2024: wisetrack01.cpp
// Given barycentric state vectors for the Sun, planets, and a spacecraft,
// calculate the 
// 

#include "solarsyst_dyn_geo01.h"
#include "cmath"
#define NUMPOS 3
#define MJDCOL 2
#define RACOL 3
#define DECCOL 4
#define MINOBSINTERVAL 1.0 // Minimum time-between-images in seconds
#define IMAGETIMETOL 1.0 // Tolerance for matching image time, in seconds
#define MAXVEL 1.5 // Default max angular velocity in deg/day.
#define MAXTIME 1.5 // Default max inter-image time interval
                    // for tracklets, in hours (will be converted
                    // to days before use).
#define IMAGERAD 2.0 // radius from image center to most distant corner (deg)
#define DEBUG 0
#define PHASE_G 0.15
#define IDNUMLEN 7

// Note: configfile contains the masses and ephemerides for all of the planets.
// configfile also contains the name of an observatory code file. The input
// units of the state vector are AU for positions and km/sec for velocities.
// Internally, the program will use AU for positions and AU/timescale for
// velocities. The default value of timescale is one year, but it can be
// set to other values as needed. The observation file must contain
// MJD, RA, Dec, magnitude, astrometric uncertainty, magnitude uncertainty,
// and obscode. The RA and Dec must be in decimal degrees. The astrometric
// uncertainty must be in arcseconds.
static void show_usage()
{
  cerr << "Usage: wisetrack01 -cfg configfile -images imfile -mjdstart mjdstart -mjdend mjdend -statevec x y z vx vy vz -projcen cenRA cenDec -platescale platescale -outfile outfile \n";
}
    
int main(int argc, char *argv[])
{
  int i=0;
  int j=0;
  int status=0;
  ifstream instream1;
  string stest;
  string configfile;
  string imfile;
  string planetfile;
  string spaceobsfile;
  string outfile;
  long double mjdstart = 0.0L;
  long double mjdend = 0.0L;
  int planetfile_startpoint=0;
  int planetfile_endpoint=0;
  int configread=0;
  int polyorder=0;
  int planetnum=0;
  int planetct=0;
  int pctSun=0;
  int pctEarth=0;
  long double ldval=0.0L;
  double dval=0.0L;
  vector <long double> planetmasses;
  vector <long double> mjdtest;
  vector <long double> obsmjd;
  vector <long double> planetmjd;
  vector <point3LD> planetpos;
  vector <point3LD> Earthpos;
  vector <point3LD> Sunpos;
  vector <point3LD> Sunvel;
  vector <point3LD> temppos;
  vector <point3LD> tempvel;
  vector <point3LD> targpos;
  vector <point3LD> targvel;
  vector <point3LD> obspos;
  vector <point3LD> obsvel;
  long double startX,startY,startZ,startVX,startVY,startVZ;
  startX = startY = startZ = startVX = startVY = startVZ = 0L;
  point3LD targ_to_obs = point3LD(0,0,0);
  long double obsdist=1.0L;
  long double light_travel_time;
  vector <point3LD> observerpos;
  vector <long double> targMJD;
  vector <double> imageRA;
  vector <double> imageDec;
  vector <long double> imageMJD;
  point3LD imagepos = point3LD(0,0,0);
  point3LD imagevel = point3LD(0,0,0);
  int imnum=0;
  int imct=0;
  string lnfromfile;
  vector <string> linestringvec;
  int badread=0;
  int reachedeof=0;
  point3LD outpos = point3LD(0,0,0);
  point3LD outvel = point3LD(0,0,0);
  long double outRA,outDec;
  double racenter,deccenter,dist,pa,xpix,ypix,platescale;
  racenter = deccenter = dist = pa = xpix = ypix = 0.0;
  platescale = 1.0;
  int fieldnum=0;
  int MJDcol = 2;
  int RAcol = 3;
  int Deccol = 4;
  ofstream outstream1;
  

  if(argc<19) {
    show_usage();
    return(1);
  }
  
  i=1;
  while(i<argc) {
    cout << "Checking out argv[" << i << "] = " << argv[i] << ".\n";
    if(string(argv[i]) == "-c" || string(argv[i]) == "-cfg" || string(argv[i]) == "-config" || string(argv[i]) == "--config" || string(argv[i]) == "--configfile" || string(argv[i]) == "--configuration" || string(argv[i]) == "--ConfigFile") {
      if(i+1 < argc) {
	//There is still something to read;
	configfile=argv[++i];
	i++;
	// Read the configuration file. This must happen here, so that
	// default values supplied in this file can be overwritten later
	// if the user desires.
	// Read configuration file.
	ifstream instream1 {configfile};
	if(!instream1) {
	  cerr << "ERROR: can't open input config file " << configfile << "\n";
	  return(1);
	}
	// Read polyorder
	status=readconfigint(instream1,&polyorder);
	while(status==1) {
	  // The line we have just read is a pure comment line,
	  // so we just want to skip to the next one.
	  status=readconfigint(instream1,&polyorder);
	}
	if(status<0) {
	  cerr << "Error reading config file\n";
	  return(1);
	} else cout << "Polynomial order for orbit integration read as " << polyorder << "\n";
	// Read the number of planets.
	status=readconfigint(instream1,&planetnum);
	while(status==1) {
	  // The line we have just read is a pure comment line,
	  // so we just want to skip to the next one.
	  status=readconfigint(instream1,&planetnum);
	}
	if(status<0) {
	  cerr << "Error reading config file\n";
	  return(1);
	} else cout << "Number of planets read as " << planetnum << "\n";
	// Read the index of the Sun within the planet vectors
	status=readconfigint(instream1,&pctSun);
	while(status==1) {
	  // The line we have just read is a pure comment line,
	  // so we just want to skip to the next one.
	  status=readconfigint(instream1,&pctSun);
	}
	if(status<0) {
	  cerr << "Error reading config file\n";
	  return(1);
	} else cout << "The Sun is planet number " << pctSun << "\n";
	// Read the index of Earth within the planet vectors
	status=readconfigint(instream1,&pctEarth);
	while(status==1) {
	  // The line we have just read is a pure comment line,
	  // so we just want to skip to the next one.
	  status=readconfigint(instream1,&pctEarth);
	}
	if(status<0) {
	  cerr << "Error reading config file\n";
	  return(1);
	} else cout << "Earth is planet number " << pctEarth << "\n";
	for(planetct=0;planetct<planetnum;planetct++) {
	  // Read the mass for planet number planetct
	  status=readconfigLD(instream1,&ldval);
	  while(status==1) {
	    // The line we have just read is a pure comment line,
	    // so we just want to skip to the next one.
	    status=readconfigLD(instream1,&ldval);
	  }
	  planetmasses.push_back(ldval);
	  if(status<0) {
	    cerr << "Error reading config file\n";
	    return(1);
	  } else {
	    cout << "MG for planet " << planetct << " read as " << planetmasses[planetct] << " km^3/sec^2\n";
	  }
	}
	for(planetct=0;planetct<planetnum;planetct++) {
	  // Read the ephemeris file for planet number planetct
	  status=readconfigstring(instream1,planetfile);
	  while(status==1) {
	    // The line we have just read is a pure comment line,
	    // so we just want to skip to the next one.
	    status=readconfigstring(instream1,planetfile);
	  }
	  if(status<0) {
	    cerr << "Error reading config file\n";
	    return(1);
	  } else cout << "Ephemeris file for planet " << planetct << " is named " << planetfile << "\n";
	  mjdtest={};
	  temppos={};
	  tempvel={};
	  read_horizons_fileLD(planetfile,mjdtest,temppos,tempvel);
	  if(planetct==0) planetmjd=mjdtest;
	  else {
	    for(j=0;j<long(planetmjd.size());j++) {
	      if(mjdtest[j]!=planetmjd[j]) {
		cout << "ERROR: time vectors do not match for input planet files\n";
		cout << planetct+1 << " and 1!\n";
		return(1);
	      }
	    }
	  }
	  for(j=0;j<long(temppos.size());j++) {
	    planetpos.push_back(temppos[j]);
	  }
	  if(planetct == pctEarth) Earthpos = temppos;
	  if(planetct == pctSun) {
	    Sunpos = temppos;
	    Sunvel = tempvel;
	  }
	}
	// Read the ephemeris file for the space observatory (e.g., WISE)
	status=readconfigstring(instream1,spaceobsfile);
	while(status==1) {
	  // The line we have just read is a pure comment line,
	  // so we just want to skip to the next one.
	  status=readconfigstring(instream1,spaceobsfile);
	}
	if(status<0) {
	  cerr << "Error reading config file\n";
	  return(1);
	} else cout << "Ephemeris file for space observatory is named " << spaceobsfile << "\n";
	obsmjd={};
	obspos={};
	obsvel={};
	read_horizons_fileLD(spaceobsfile,obsmjd,obspos,obsvel);
	cout << "Finished reading ephemeris file " << spaceobsfile << "\n";
	// Read default starting MJD
	status=readconfigLD(instream1,&mjdstart);
	while(status==1) {
	  // The line we have just read is a pure comment line,
	  // so we just want to skip to the next one.
	  status=readconfigLD(instream1,&mjdstart);
	}
	if(status<0) {
	  cerr << "Error reading config file\n";
	  return(1);
	} else cout << "Default starting MJD read as " << mjdstart << "\n";
	// Read default ending MJD
	status=readconfigLD(instream1,&mjdend);
	while(status==1) {
	  // The line we have just read is a pure comment line,
	  // so we just want to skip to the next one.
	  status=readconfigLD(instream1,&mjdend);
	}
	if(status<0) {
	  cerr << "Error reading config file\n";
	  return(1);
	} else cout << "Default ending MJD read as " << mjdend << "\n";
	// Close input stream that was reading the config file.
	instream1.close();
	configread=1;
      } else {
	cerr << "Configuration file keyword supplied with no corresponding argument\n";
	show_usage();
	return(1);
      }
    }
    // Stop execution if we haven't read a config file successfully by now.
    if(configread!=1) {
      cerr << "ERROR: configuration file must be supplied and successfully";
      cerr << "read before any of the other parameters. This is because\n";
      cerr << "it contains defaults that can optionally be overridden by\n";
      cerr << "user-supplied values in later arguments.\n";
      return(2);
    } else cout << "Configuration file read successfully\n";
    if(string(argv[i]) == "-i" || string(argv[i]) == "-im" || string(argv[i]) == "-imf" || string(argv[i]) == "-images" || string(argv[i]) == "-imfile" || string(argv[i]) == "--imagefile" || string(argv[i]) == "--imfile" || string(argv[i]) == "--images") {
      if(i+1 < argc) {
	// There is still something to read;
	imfile=argv[++i];
	i++;
      } else {
	cerr << "Image file keyword supplied with no corresponding argument\n";
	show_usage();
	return(1);
      }
    } else if(string(argv[i]) == "-mjdstart" || string(argv[i]) == "-ms" || string(argv[i]) == "-MJDstart" || string(argv[i]) == "--mjdstart" || string(argv[i]) == "--MJDstart" || string(argv[i]) == "--ModifiedJulianDaystart" || string(argv[i]) == "--modifiedjuliandaystart") {
      if(i+1 < argc) {
	//There is still something to read;
	mjdstart=stold(argv[++i]);
	i++;
      }
      else {
	cerr << "MJD start keyword supplied with no corresponding argument\n";
	show_usage();
	return(1);
      }
    } else if(string(argv[i]) == "-mjdend" || string(argv[i]) == "-me" || string(argv[i]) == "-MJDend" || string(argv[i]) == "--mjdend" || string(argv[i]) == "--MJDend" || string(argv[i]) == "--ModifiedJulianDayend" || string(argv[i]) == "--modifiedjuliandayend") {
      if(i+1 < argc) {
	//There is still something to read;
	mjdend=stold(argv[++i]);
	i++;
      }
      else {
	cerr << "MJD end keyword supplied with no corresponding argument\n";
	show_usage();
	return(1);
      }
    } else if(string(argv[i]) == "-mjdend" || string(argv[i]) == "-me" || string(argv[i]) == "-MJDend" || string(argv[i]) == "--mjdend" || string(argv[i]) == "--MJDend" || string(argv[i]) == "--ModifiedJulianDayend" || string(argv[i]) == "--modifiedjuliandayend") {
      if(i+1 < argc) {
	//There is still something to read;
	mjdend=stold(argv[++i]);
	i++;
      }
      else {
	cerr << "MJD end keyword supplied with no corresponding argument\n";
	show_usage();
	return(1);
      }
    } else if(string(argv[i]) == "-statevec" || string(argv[i]) == "-sv" || string(argv[i]) == "-svec" || string(argv[i]) == "--statevector" || string(argv[i]) == "--statevec" || string(argv[i]) == "-statevector" || string(argv[i]) == "--sv") {
      if(i+1 < argc) {
	//There is still something to read;
	startX=stold(argv[++i]);
	cout << "read sv entry as " << startX << "\n";
	if(i+1 >= argc) {
	  cerr << "State vector keyword supplied with fewer than the required 6 arguments\n";
	  show_usage();
	  return(1);
	}
	startY=stold(argv[++i]);
	cout << "read sv entry as " << startY << "\n";
	if(i+1 >= argc) {
	  cerr << "State vector keyword supplied with fewer than the required 6 arguments\n";
	  show_usage();
	  return(1);
	}
	startZ=stold(argv[++i]);
	cout << "read sv entry as " << startZ << "\n";
	if(i+1 >= argc) {
	  cerr << "State vector keyword supplied with fewer than the required 6 arguments\n";
	  show_usage();
	  return(1);
	}
 	startVX=stold(argv[++i]);
	cout << "read sv entry as " << startVX << "\n";
	if(i+1 >= argc) {
	  cerr << "State vector keyword supplied with fewer than the required 6 arguments\n";
	  show_usage();
	  return(1);
	}
	startVY=stold(argv[++i]);
	cout << "read sv entry as " << startVY << "\n";
	if(i+1 >= argc) {
	  cerr << "State vector keyword supplied with fewer than the required 6 arguments\n";
	  show_usage();
	  return(1);
	}
	startVZ=stold(argv[++i]);
	cout << "read sv entry as " << startVZ << "\n";
	i++;
      }
      else {
	cerr << "State vector keyword supplied with no corresponding argument\n";
	show_usage();
	return(1);
      }
    } else if(string(argv[i]) == "-projcen" || string(argv[i]) == "-pcen" || string(argv[i]) == "-projcenter" || string(argv[i]) == "--projcen" || string(argv[i]) == "--projection_center" || string(argv[i]) == "-projcenter" || string(argv[i]) == "-pc") {
      if(i+1 < argc) {
	//There is still something to read;
	racenter=stod(argv[++i]);
	if(i+1 >= argc) {
	  cerr << "projection center keyword supplied with fewer than the required 2 arguments\n";
	  show_usage();
	  return(1);
	}
	deccenter=stod(argv[++i]);
	i++;
      }
      else {
	cerr << "projection center keyword supplied with no corresponding argument\n";
	show_usage();
	return(1);
      }
    } else if(string(argv[i]) == "-platescale" || string(argv[i]) == "-pixscale" || string(argv[i]) == "--platescale" || string(argv[i]) == "--pixscale") {
      if(i+1 < argc) {
	//There is still something to read;
	platescale=stod(argv[++i]);
	i++;
      }
      else {
	cerr << "Output file keyword supplied with no corresponding argument\n";
	show_usage();
	return(1);
      }
    } else if(string(argv[i]) == "-out" || string(argv[i]) == "-outfile" || string(argv[i]) == "-o" || string(argv[i]) == "--outfile" || string(argv[i]) == "--outorb" || string(argv[i]) == "--outorbits") {
      if(i+1 < argc) {
	//There is still something to read;
	outfile=argv[++i];
	i++;
      }
      else {
	cerr << "Output file keyword supplied with no corresponding argument\n";
	show_usage();
	return(1);
      }
    } else {
      cerr << "Warning: unrecognized keyword or argument " << argv[i] << "\n";
      i++;
    }
  }

  cout.precision(17);  
  cout << "input configuration file " << configfile << "\n";
  cout << "input image file " << imfile << "\n";
  cout << "input starting MJD " << mjdstart << "\n";
  cout << "input ending MJD " << mjdend << "\n";
  cout << "output file " << outfile << "\n";


  // Match mjdstart and mjdend to planet file.
  
  planetfile_startpoint = planetfile_endpoint = -99;
  for(j=0;j<long(planetmjd.size());j++) {
    if(fabs(planetmjd[j]-mjdstart) < IMAGETIMETOL/SOLARDAY) planetfile_startpoint = j;
    if(fabs(planetmjd[j]-mjdend) < IMAGETIMETOL/SOLARDAY) planetfile_endpoint = j;
  }
  if(planetfile_startpoint<0 || planetfile_endpoint<0) {
    cerr << "ERROR: mjdstart " << mjdstart << " and/or mjdend " << mjdend << " could not be matched\nto any timestep in the planet files\n";
    return(1);
  } else {
    cout << "mjdstart " << mjdstart << " corresponds to timestep " << planetfile_startpoint << " in the planet files\n";
    cout << "mjdend " << mjdend << " corresponds to timestep " << planetfile_endpoint << " in the planet files\n";
  }
  
  // Read input image file.
  instream1.open(imfile);
  if(!instream1) {
    cerr << "ERROR: unable to open input file " << imfile << "\n";
    return(1);
  }
  badread=0;
  while(reachedeof==0 && !instream1.eof() && !instream1.fail() && !instream1.bad()) {
    getline(instream1,lnfromfile);
    if(!instream1.eof() && !instream1.fail() && !instream1.bad()) ; // Read on.
    else if(instream1.eof()) reachedeof=1; //End of file, fine.
    else if(instream1.fail()) reachedeof=-1; //Something wrong, warn
    else if(instream1.bad()) reachedeof=-2; //Worse problem, warn
    if(lnfromfile.size()>10 && (reachedeof==0 || reachedeof==1)) {
      fieldnum = get_col_vector01(lnfromfile, linestringvec);
      if(fieldnum>=MJDcol) {
	ldval = stold(linestringvec[MJDcol-1]);
	imageMJD.push_back(ldval);
      } else badread=1;
      if(fieldnum>=RAcol) {
	dval = stod(linestringvec[RAcol-1]);
	imageRA.push_back(dval);
      } else badread=1;
      if(fieldnum>=Deccol) {
	dval = stod(linestringvec[Deccol-1]);
	imageDec.push_back(dval);
      } else badread=1;
    }
    cout << imageMJD.size() << " " << imageMJD[imageMJD.size()-1] << " " << imageRA[imageMJD.size()-1] << " " << imageDec[imageMJD.size()-1] << "\n";
  }
  imnum=imageMJD.size();
  if(imnum<2 || long(imageRA.size())!=imnum || long(imageDec.size())!=imnum) {
    cerr  << fixed << setprecision(6) << "Error: image time vectors too short, or of unequal length:\n";
    cerr  << fixed << setprecision(6) << imnum << " " << imageRA.size() << " " << imageDec.size() << "\n";
    return(1);
  }
  if(badread>=1){
    cerr << "Error reading " << imfile << ": apparent short line\n";
    return(1);
  }
  if(reachedeof==1) { 
    cout << "File read successfully to the end.\n";
  }
  else if(reachedeof==-1) cout << "Warning: file read failed\n";
  else if(reachedeof==-2) cout << "Warning: file possibly corrupted\n";
  else cout << "Warning: unknown file read problem\n";
  instream1.close();
  
  // Calculate the exact position of the observer at the time of each image.
  observerpos={};
  for(imct=0;imct<imnum;imct++) {
    //    observer_barycoords01LD(imageMJD[imct], 5, 0L, 0L, 0L, obsmjd, obspos, outpos);
    planetposvel01LD(imageMJD[imct], polyorder, obsmjd, obspos, obsvel, outpos, outvel);
    observerpos.push_back(outpos);
  }

  // Open output file and write header
  outstream1.open(outfile);
  outstream1 << "imageMJD outRA outDec xpix ypix observerX observerY observerZ targX targY targZ targVX targVY targVZ \n";
  outpos = point3LD(startX,startY,startZ);
  outvel = point3LD(startVX,startVY,startVZ);
  // Perform n-body integration beginning with the input state vectors.
  integrate_orbit04LD(polyorder, planetnum, planetmjd, planetmasses, planetpos, outpos, outvel, planetfile_startpoint, planetfile_endpoint, targMJD, targpos, targvel);
  for(imct=0;imct<imnum;imct++)
    {
      planetposvel01LD(imageMJD[imct],polyorder,targMJD,targpos,targvel,imagepos,imagevel);
      startX = imagepos.x;
      startY = imagepos.y;
      startZ = imagepos.z;
      targ_to_obs.x = observerpos[imct].x - imagepos.x;
      targ_to_obs.y = observerpos[imct].y - imagepos.y;
      targ_to_obs.z = observerpos[imct].z - imagepos.z;
      obsdist = sqrt(targ_to_obs.x*targ_to_obs.x + targ_to_obs.y*targ_to_obs.y + targ_to_obs.z*targ_to_obs.z);
      // Instantaneous velocity
      planetpos01LD(imageMJD[imct],polyorder,targMJD,targvel,imagevel);
      // Initial approximation of the coordinates relative to the observer
      imagepos.x -= observerpos[imct].x;
      imagepos.y -= observerpos[imct].y;
      imagepos.z -= observerpos[imct].z;
      light_travel_time = obsdist*1000.0/CLIGHT; // Factor of 1000 converts obsdist to meters
      // Light-travel-time corrected version of coordinates relative to the observer
      imagepos.x -= light_travel_time*imagevel.x;
      imagepos.y -= light_travel_time*imagevel.y;
      imagepos.z -= light_travel_time*imagevel.z;
      // Light-travel-time corrected observer-target distance
      obsdist = sqrt(imagepos.x*imagepos.x + imagepos.y*imagepos.y + imagepos.z*imagepos.z);
      // Calculate unit vector
      imagepos.x /= obsdist;
      imagepos.y /= obsdist;
      imagepos.z /= obsdist;
      // Project onto the celestial sphere.
      stateunitLD_to_celestial(imagepos, outRA, outDec);
      // Carry out stereographic projection
      distradec02(racenter,deccenter,outRA,outDec,&dist,&pa);
      dist = 2.0l*tan(dist/DEGPRAD/2.0l);
      xpix = -dist*sin(pa/DEGPRAD)*ASECPRAD/platescale; // East is left
      ypix = dist*cos(pa/DEGPRAD)*ASECPRAD/platescale;
      outstream1 << fixed << setprecision(6) << imageMJD[imct] << " " << outRA << " " << outDec << " ";
      outstream1 << fixed << setprecision(3) << xpix << " " << ypix << " ";
      outstream1 << fixed << setprecision(6) << observerpos[imct].x << " " << observerpos[imct].y << " " << observerpos[imct].z << " ";
      outstream1 << fixed << setprecision(6) << startX << " " << startY << " " << startZ << " ";
      outstream1 << fixed << setprecision(6) << imagevel.x << " " << imagevel.y << " " << imagevel.z << "\n";
    }
  outstream1.close();
  
  return(0);
}
