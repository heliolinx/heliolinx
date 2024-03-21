/// March 15, 2024: thorhyp01.cpp
// Intended to calculate state vectors for input into successive runs of
// wisetrack01.cpp. Given observer coordinates in the barycentric ecliptic
// system, the corresponding MJD, a celestial position in RA, Dec,
// a distance from the observer in AU, and velocity sampling intervals
// in km/sec in counterclockwise circular speed relative to the Sun,
// the +z dimension of the tangential velocity relative to the Sun,
// and the radial velocity relative to the Sun -- plus the number of
// sample steps to probe in each of these three dimensions, calculate
// ephemerides a state vector at a set earlier time that will pass
// through the exact specified RA, Dec position, as viewed from the
// observer's location, at the required target-observer distance and
// with the specified heliocentric velocity.

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
  cerr << "Usage: thorhyp01 -cfg configfile -obscoords observer barycentric X Y Z -mjdref mjdref -mjdstart mjdstart -radec RA Dec -obsdist dist_from_observer(AU) -circvelsamp circvelsamp -zvelsamp zvelsamp -radvelsamp radvelsamp -circvelnum circvelnum -zvelnum zvelnum -radvelnum radvelnum -outfile outfile \n";
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
  string outfile;
  long double mjdstart = 0.0L;
  long double mjdref = 0.0L;
  long double ldval = 0.0L;
  int configread=0;
  int polyorder=0;
  int planetnum=0;
  int planetct=0;
  int pctSun=0;
  int pctEarth=0;
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
  vector <point3LD> obspos;
  vector <point3LD> obsvel;
  point3LD unitcirc = point3LD(0,0,0);
  point3LD unitrad = point3LD(0,0,0);
  point3LD unitz = point3LD(0,0,0);
  long double startX,startY,startZ,circvel,xynorm;
  startX = startY = startZ = 0L;
  long double circvelsamp,zvelsamp,radvelsamp;
  circvelsamp = zvelsamp = radvelsamp = 0L;
  long circvelnum,zvelnum,radvelnum;
  circvelnum = zvelnum = radvelnum = 0;
  point3LD refpos = point3LD(0,0,0);
  point3LD refvel = point3LD(0,0,0);
  point3LD targpos = point3LD(0,0,0);
  point3LD targvel = point3LD(0,0,0);
  long double obsdist=1.0L;
  long double sundist=1.0L;
  vector <point3LD> observerpos;
  long double racenter,deccenter;
  racenter = deccenter = 0L;
  ofstream outstream1;
  long double a,incl,e;
  

  if(argc<30) {
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
	  status = read_horizons_fileLD(planetfile,mjdtest,temppos,tempvel);
	  if(status!=0) {
	    cerr << "ERROR: attept to read ephemeris file " << planetfile << " failed with status " << status << "\n";
	    return(status);
	  }
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
	// Read default reference MJD
	status=readconfigLD(instream1,&mjdref);
	while(status==1) {
	  // The line we have just read is a pure comment line,
	  // so we just want to skip to the next one.
	  status=readconfigLD(instream1,&mjdref);
	}
	if(status<0) {
	  cerr << "Error reading config file\n";
	  return(1);
	} else cout << "Default ending MJD read as " << mjdref << "\n";
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
    if(string(argv[i]) == "-mjdref" || string(argv[i]) == "-me" || string(argv[i]) == "-Mjdref" || string(argv[i]) == "--mjdref" || string(argv[i]) == "--Mjdref" || string(argv[i]) == "--ModifiedJulianDayend" || string(argv[i]) == "--modifiedjuliandayend") {
      if(i+1 < argc) {
	//There is still something to read;
	mjdref=stold(argv[++i]);
	i++;
      }
      else {
	cerr << "MJD end keyword supplied with no corresponding argument\n";
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
    } else if(string(argv[i]) == "-mjdref" || string(argv[i]) == "-me" || string(argv[i]) == "-Mjdref" || string(argv[i]) == "--mjdref" || string(argv[i]) == "--Mjdref" || string(argv[i]) == "--ModifiedJulianDayref" || string(argv[i]) == "--modifiedjuliandayref") {
      if(i+1 < argc) {
	//There is still something to read;
	mjdref=stold(argv[++i]);
	i++;
      }
      else {
	cerr << "Reference MJD keyword supplied with no corresponding argument\n";
	show_usage();
	return(1);
      }
    } else if(string(argv[i]) == "-obscoords" || string(argv[i]) == "-obspos" || string(argv[i]) == "-obsbary" || string(argv[i]) == "--obscoords" || string(argv[i]) == "--obsbary" || string(argv[i]) == "-obsstatevector" || string(argv[i]) == "--sv") {
      if(i+1 < argc) {
	//There is still something to read;
	startX=stold(argv[++i]);
	cout << "read sv entry as " << startX << "\n";
	if(i+1 >= argc) {
	  cerr << "Observer position keyword supplied with fewer than the required 3 arguments\n";
	  show_usage();
	  return(1);
	}
	startY=stold(argv[++i]);
	cout << "read sv entry as " << startY << "\n";
	if(i+1 >= argc) {
	  cerr << "Observer position keyword supplied with fewer than the required 3 arguments\n";
	  show_usage();
	  return(1);
	}
	startZ=stold(argv[++i]);
	cout << "read sv entry as " << startZ << "\n";
	i++;
      } else {
	cerr << "Observer position keyword supplied with no corresponding argument\n";
	show_usage();
	return(1);
      }
    } else if(string(argv[i]) == "-radec" || string(argv[i]) == "-pcen" || string(argv[i]) == "-projcenter" || string(argv[i]) == "--projcen" || string(argv[i]) == "--projection_center" || string(argv[i]) == "-projcenter" || string(argv[i]) == "-pc") {
      if(i+1 < argc) {
	//There is still something to read;
	racenter=stold(argv[++i]);
	if(i+1 >= argc) {
	  cerr << "projection center keyword supplied with fewer than the required 2 arguments\n";
	  show_usage();
	  return(1);
	}
	deccenter=stold(argv[++i]);
	i++;
      }
      else {
	cerr << "projection center keyword supplied with no corresponding argument\n";
	show_usage();
	return(1);
      }
    } else if(string(argv[i]) == "-obsdist" || string(argv[i]) == "-observerdist" || string(argv[i]) == "--obsdist" || string(argv[i]) == "-geodist") {
      if(i+1 < argc) {
	//There is still something to read;
	obsdist=stold(argv[++i]);
	i++;
      }
      else {
	cerr << "Observer distance keyword supplied with no corresponding argument\n";
	show_usage();
	return(1);
      }
     
    } else if(string(argv[i]) == "-circvelsamp") {
      if(i+1 < argc) {
	//There is still something to read;
	circvelsamp = stold(argv[++i]);
	i++;
      }
      else {
	cerr << "Circular velocity sampling keyword supplied with no corresponding argument\n";
	show_usage();
	return(1);
      }
    } else if(string(argv[i]) == "-zvelsamp") {
      if(i+1 < argc) {
	//There is still something to read;
	zvelsamp = stold(argv[++i]);
	i++;
      }
      else {
	cerr << "Z-direction velocity keyword supplied with no corresponding argument\n";
	show_usage();
	return(1);
      }
    } else if(string(argv[i]) == "-radvelsamp") {
      if(i+1 < argc) {
	//There is still something to read;
	radvelsamp = stold(argv[++i]);
	i++;
      }
      else {
	cerr << "Radial velocity keyword supplied with no corresponding argument\n";
	show_usage();
	return(1);
      }
    } else if(string(argv[i]) == "-circvelnum") {
      if(i+1 < argc) {
	//There is still something to read;
	circvelnum = stol(argv[++i]);
	i++;
      }
      else {
	cerr << "Circular velocity step number keyword supplied with no corresponding argument\n";
	show_usage();
	return(1);
      }
    } else if(string(argv[i]) == "-zvelnum") {
      if(i+1 < argc) {
	//There is still something to read;
	zvelnum = stol(argv[++i]);
	i++;
      }
      else {
	cerr << "Z-direction velocity step number keyword supplied with no corresponding argument\n";
	show_usage();
	return(1);
      }
    } else if(string(argv[i]) == "-radvelnum") {
      if(i+1 < argc) {
	//There is still something to read;
	radvelnum = stol(argv[++i]);
	i++;
      }
      else {
	cerr << "Radial velocity step number keyword supplied with no corresponding argument\n";
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
  cout << "input reference MJD " << mjdstart << "\n";
  cout << "input starting MJD " << mjdref << "\n";
  cout << "output file " << outfile << "\n";

  // Calculate reference position
  celestial_to_statevecLD(racenter,deccenter,obsdist*AU_KM,refpos);
  // Add in observer position
  refpos.x += startX;
  refpos.y += startY;
  refpos.z += startZ;
  // Calculate reference velocity: a circular counterclockwise orbit with no Z-component.
  sundist = vecabs3LD(refpos);
  circvel = sqrt(GMSUN_KM3_SEC2/sundist);
  refvel.z = 0L;
  xynorm = sqrt(refpos.x*refpos.x + refpos.y*refpos.y);
  if(xynorm==0L) {
    // We are on the polar axis of the solar system; not responsible
    // for weird behavior. Just make it go full speed in the +X direction.
    refvel = point3LD(circvel,0L,0L);
  } else {
    // Circular counterclockwise orbit with no Z-component to the velocity at this instant.
    refvel = point3LD(-circvel*refpos.y/xynorm,circvel*refpos.x/xynorm,0L);
  }
  // Calculate orthogonal unit vectors for velocity system.
  if(xynorm==0L) {
    unitcirc = point3LD(1.0L,0L,0L);
    unitz = point3LD(0L,1.0L,0L);
    unitrad = point3LD(0L,0L,1.0L);
  } else {
    unitcirc = point3LD(-refpos.y/xynorm,refpos.x/xynorm,0L);
    unitz = point3LD(-refpos.x/xynorm*refpos.z/sundist,-refpos.y/xynorm*refpos.z/sundist,xynorm/sundist);
    unitrad = point3LD(refpos.x/sundist,refpos.y/sundist,refpos.z/sundist);
  }
  // Check orthogonality and normalization:
  cout << "Check 1 1 1 0 0 0: " << vecabs3LD(unitcirc) << " " << vecabs3LD(unitz) << " " << vecabs3LD(unitrad) << " " << dotprod3LD(unitcirc,unitz) << " " << dotprod3LD(unitcirc,unitrad) << " " << dotprod3LD(unitz,unitrad) << "\n";
  
  // Open output file and write header
  outstream1.open(outfile);
  outstream1 << "circvelnum zvelnum radvelnum refvelX refvelY refvelZ X Y Z VX VY VZ a e incl\n";
  // Loop on velocities
  for(long cvi = -circvelnum ; cvi <= circvelnum ; cvi++) {
    for(long zvi = -zvelnum ; zvi <= zvelnum ; zvi++) {
      for(long rvi = -radvelnum ; rvi <= radvelnum ; rvi++) {
	targvel.x = refvel.x + double(cvi)*circvelsamp*unitcirc.x + double(zvi)*zvelsamp*unitz.x + double(rvi)*radvelsamp*unitrad.x;
	targvel.y = refvel.y + double(cvi)*circvelsamp*unitcirc.y + double(zvi)*zvelsamp*unitz.y + double(rvi)*radvelsamp*unitrad.y;
	targvel.z = refvel.z + double(cvi)*circvelsamp*unitcirc.z + double(zvi)*zvelsamp*unitz.z + double(rvi)*radvelsamp*unitrad.z;
	targpos = refpos;
	vector <long double> targMJD = {};
	vector <point3LD> outpos = {};
	vector <point3LD> outvel = {};
	targMJD.push_back(mjdstart);
	// Perform n-body integration beginning with the input state vector.
	integrate_orbit03LD(polyorder, planetnum, planetmjd, planetmasses, planetpos, targMJD, targpos, targvel, mjdref, outpos, outvel);
	outstream1.precision(12);  
	outstream1 << cvi << " " << zvi << " " << rvi << " " << targvel.x << " " << targvel.y << " " << targvel.z << " ";
	outstream1 << outpos[0].x << " " << outpos[0].y << " " << outpos[0].z << " ";
	outstream1 << outvel[0].x << " " << outvel[0].y << " " << outvel[0].z << " ";

	long double v0 = vecabs3LD(targvel);
	long double u = dotprod3LD(targvel,targpos);
 	a = sundist*GMSUN_KM3_SEC2/(2.0l*GMSUN_KM3_SEC2-v0*v0*sundist);
	long double n = sqrt(GMSUN_KM3_SEC2/a/a/a);
	long double alpha = GMSUN_KM3_SEC2/a;
	if(alpha>0.0l) {
	  // Bound, elliptical orbit
	  long double EC = 1.0L - sundist/a;
	  long double ES = u/n/a/a;
	  e = sqrt(EC*EC + ES*ES);
	} else if (alpha<0.0l) {
	  // Unbound, hyperbolic orbit
	  long double CH = 1.0L - sundist/a;
	  long double SH = u/sqrt(-GMSUN_KM3_SEC2*a);
	  e = sqrt(CH*CH - SH*SH);
	}
	point3LD lvec = crossprod3LD(targpos,targvel);
	vecnorm3LD(lvec);
	incl = DEGPRAD*acos(lvec.z);
	if(lvec.z > 1.0L) incl = 90.0L;
	else if(lvec.z < -1.0L) incl = -90.0L;
	outstream1.precision(6);  
	outstream1 << a/AU_KM << " " << e << " " << incl << "\n";
      }
    }
  }
  outstream1.close();
  return(0);
}
