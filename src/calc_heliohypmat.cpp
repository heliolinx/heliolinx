// June 13, 2024: calc_heliohypmat.cpp
// Like the calc_accel_mat programs, but uses a very different approach
// to the problem of calculating a hypothesis matrix for heliolinc.
// Interpolates an input file giving sampling in each dimension,
// as a function of heliocentric distance.

#include "solarsyst_dyn_geo01.h"
#include "cmath"

#define STEPMAX 10000000

static void show_usage()
{
  cerr << "Usage: calc_heliohypmat -interp interpolation file -mindist minimum distance (AU) -maxdist maximum distance (AU) -distfac factor by which to scale distance sampling -velfac factor by which to scale velocity sampling -accfac factor by which to scale acceleration sampling -middle usemiddle -outfile output file\n";

}

int main(int argc, char *argv[])
{
  double mindist = 1.5l;
  double maxdist = 100.0l;
  double diststep = 1.0l;
  double velstep = 1.0e-2;
  double distfac,velfac,accfac;
  distfac = velfac = accfac = 1.0l;
  long i=0;
  long stepct=0;
  int verbose=0;
  ifstream instream1;
  ofstream outstream1;
  double dist,velkm,velAU,vesc,accelk,accelnorm,g0;
  dist = velkm = velAU = vesc = accelk = accelnorm = g0 = 0.0l;
  double distkm=0.0l;
  double accelstep,minacc,maxacc,vtan,xi;
  accelstep = minacc = maxacc = vtan = xi = 0.0l;
  long accelnum,accelct,velnum,velct;
  string interpfile,outfile;
  int usemiddle=0;
  double ds,vs,as;
  vector <double> distvec,distsampvec,velsampvec,accsampvec;
  long ilo,ihi;
  string lnfromfile;
  
  if(argc<5) {
    show_usage();
    return(1);
  }
  i=1;
  while(i<argc) {
    cout << "Checking out argv[" << i << "] = " << argv[i] << ".\n";
    if(string(argv[i]) == "-interp" || string(argv[i]) == "-interpfile" || string(argv[i]) == "-infile" || string(argv[i]) == "-interpolation" || string(argv[i]) == "--interp" || string(argv[i]) == "--interpfile" || string(argv[i]) == "--interpolation" || string(argv[i]) == "--interpolation_file") {
      if(i+1 < argc) {
	//There is still something to read;
	interpfile = argv[++i];
	i++;
      }
      else {
	cerr << "Interpolation file keyword supplied with no corresponding argument\n";
	show_usage();
	return(1);
      }
    } else if(string(argv[i]) == "-mindist" || string(argv[i]) == "-mind" || string(argv[i]) == "-mindistance" || string(argv[i]) == "-minimumdist" || string(argv[i]) == "--minimum_distance" || string(argv[i]) == "--mindist" || string(argv[i]) == "--mind" || string(argv[i]) == "--mindist") {
      if(i+1 < argc) {
	//There is still something to read;
	mindist = stod(argv[++i]);
	i++;
      }
      else {
	cerr << "Input minimum heliocentric distance keyword supplied with no corresponding argument\n";
	show_usage();
	return(1);
      }
    } else if(string(argv[i]) == "-maxdist" || string(argv[i]) == "-maxd" || string(argv[i]) == "-maxdistance" || string(argv[i]) == "-maximumdist" || string(argv[i]) == "--maximum_distance" || string(argv[i]) == "--maxdist" || string(argv[i]) == "--maxd" || string(argv[i]) == "--maxdist") {
      if(i+1 < argc) {
	//There is still something to read;
	maxdist = stod(argv[++i]);
	i++;
      }
      else {
	cerr << "Input maximum heliocentric distance keyword supplied with no corresponding argument\n";
	show_usage();
	return(1);
      }
    } else if(string(argv[i]) == "-distfac" || string(argv[i]) == "--distfac" ) {
      if(i+1 < argc) {
	//There is still something to read;
	distfac = stod(argv[++i]);
	i++;
      }
      else {
	cerr << "Scaling factor for distance sampling keyword supplied with no corresponding argument\n";
	show_usage();
	return(1);
      }
    } else if(string(argv[i]) == "-velfac" || string(argv[i]) == "--velfac") {
      if(i+1 < argc) {
	//There is still something to read;
	velfac=stod(argv[++i]);
	i++;
      }
      else {
	cerr << "Scaling factor for velocity sampling keyword supplied with no corresponding argument\n";
	show_usage();
	return(1);
      }
    } else if(string(argv[i]) == "-accfac" || string(argv[i]) == "-acfac" || string(argv[i]) == "--acfac" || string(argv[i]) == "--accfac" ) {
      if(i+1 < argc) {
	//There is still something to read;
	accfac = stod(argv[++i]);
	i++;
      }
      else {
	cerr << "Scaling factor for acceleration sampling keyword supplied with no corresponding argument\n";
	show_usage();
	return(1);
      }
    } else if(string(argv[i]) == "-middle" || string(argv[i]) == "-usemiddle" || string(argv[i]) == "-hitmiddle") {
      if(i+1 < argc) {
	//There is still something to read;
	usemiddle = stoi(argv[++i]); 
	i++;
      }
      else {
	cerr << "Middle keyword supplied with no corresponding argument\n";
	show_usage();
	return(1);
      }
    } else if(string(argv[i]) == "-out" || string(argv[i]) == "-outfile" || string(argv[i]) == "--out" || string(argv[i]) == "--outfile") {
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
    }  else if(string(argv[i]) == "-verbose" || string(argv[i]) == "-verb" || string(argv[i]) == "-VERBOSE" || string(argv[i]) == "-VERB" || string(argv[i]) == "--verbose" || string(argv[i]) == "--VERBOSE" || string(argv[i]) == "--VERB") {
      if(i+1 < argc) {
	//There is still something to read;
	verbose=stoi(argv[++i]);
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

  
  cout << "Distances will be probed from " << mindist << " to " << maxdist << "AU\n";
  cout << "Output file to be written is called " << outfile << "\n";

  // Catch required parameters if missing
  if(interpfile.size()<=0) {
    cout << "\nERROR: input interpolation file is required\n";
    show_usage();
    return(1);
  } else if(mindist<=0.0l) {
    cout << "\nERROR: minimum distance must be strictly positive\n";
    show_usage();
    return(1);
  } else if(maxdist<=0.0l || maxdist<=mindist) {
    cout << "\nERROR: maximum distnace must be strictly greater than minimum distance\n";
    show_usage();
    return(1);
  } else if(outfile.size()<=0) {
    cout << "\nERROR: output filename is required\n";
    show_usage();
    return(1);
  }

  // Read interpolation file
  instream1.open(interpfile);
  getline(instream1,lnfromfile); // Skips one-line header
  while(!instream1.bad() && !instream1.fail() && !instream1.eof()) {
    instream1 >> dist >> ds >> vs >> as;
    distvec.push_back(dist);
    distsampvec.push_back(ds);
    velsampvec.push_back(vs);
    accsampvec.push_back(as);
  }
  instream1.close();
  long interpnum = distvec.size();
  if(long(distsampvec.size())!=interpnum || long(velsampvec.size())!=interpnum || long(accsampvec.size())!=interpnum) {
    cerr << "ERROR: vectors read from input file do not all have\n";
    cerr << "the same length: " << interpnum << " " << distsampvec.size() << " " << velsampvec.size() << " " << accsampvec.size() << "\n";
  }
  cout << "Read " << interpnum << " lines from file " << interpfile << "\n";
  for(i=0;i<interpnum;i++) {
    cout << "Interp file: " << distvec[i] << " " << distsampvec[i] << " " << velsampvec[i] << " " << accsampvec[i] << "\n";
  }
  
  outstream1.open(outfile);
  outstream1 << "#r(AU) rdot(AU/day) mean_accel\n";
  dist = mindist;
  while(dist<=maxdist && stepct<=STEPMAX) {
    if(dist<distvec[0]) {
      // We're off the lower end of the interpolation file.
      // Just set all the step sizes to the ones corresponding
      // to the smallest heliocentric distance covered by the file.
      diststep = distsampvec[0];
      velstep = velsampvec[0];
      accelstep = accsampvec[0];
    } else if(dist>distvec[interpnum-1]) {
      // We're off the upper end of the interpolation file.
      // Just set all the step sizes to the ones corresponding
      // to the largest heliocentric distance covered by the file.
      diststep = distsampvec[interpnum-1];
      velstep = velsampvec[interpnum-1];
      accelstep = accsampvec[interpnum-1];
    } else {
      // Find which two points we are interpolating between
      ihi=0;
      while(ihi < interpnum && distvec[ihi]<=dist) ihi++;
      if(ihi<=0) {
	cerr << "Logically impossible case at dist = " << dist << ", distvec[0] = " << distvec[0] << "\n";
	return(1);
      }
      ilo=ihi-1;
      diststep = exp(log(distsampvec[ilo]) + (log(distsampvec[ihi])-log(distsampvec[ilo]))*(log(dist)-log(distvec[ilo]))/(log(distvec[ihi])-log(distvec[ilo])));
      velstep = exp(log(velsampvec[ilo]) + (log(velsampvec[ihi])-log(velsampvec[ilo]))*(log(dist)-log(distvec[ilo]))/(log(distvec[ihi])-log(distvec[ilo])));
      accelstep = exp(log(accsampvec[ilo]) + (log(accsampvec[ihi])-log(accsampvec[ilo]))*(log(dist)-log(distvec[ilo]))/(log(distvec[ihi])-log(distvec[ilo])));
    }
    diststep *= distfac;
    velstep *= velfac;
    accelstep *= accfac;
    
    distkm = dist*AU_KM;
    vesc = sqrt(2.0l*GMSUN_KM3_SEC2/distkm)*SOLARDAY/AU_KM; // AU/day
    g0 = GMSUN_KM3_SEC2/distkm/distkm; // Solar gravity in km/sec^2
    accelstep *= g0;
      
    if(verbose>=1) cout << "Distance is " << dist << " AU = " << distkm << " km, vesc = " << vesc*AU_KM/SOLARDAY << " km/sec, velstep = " << velstep*AU_KM/SOLARDAY << "km/sec, accelstep = " << accelstep << " km/sec^2, GMsun = " << g0 << " km/sec^2\n";
    velnum = ceil(2.0l*vesc/velstep);
    if(usemiddle==1 && velnum>1 && velnum%2==0) {
      // We are currently set to probe an even number of velocity
      // points. This means we won't probe the v=0 midpoint, which
      // is often the most probable. Add one additional point, to
      // ensure that v=0 gets probed.
      velnum+=1;
    }
    if(velnum<=1) {
      // Velocity step is so large, we'll only probe radial velocity=0
      // Since radial velocity is zero, tangential velocity can be anywhere from 0 to vesc.
      // WE TAKE POSITIVE ACCELERATION AS INWARD TOWARD THE SUN
      maxacc = g0;
      minacc = -g0;
      velkm = 0.0;
      velAU = velkm*SOLARDAY/AU_KM;
      if(verbose>=1) cout << "vel = " << velAU << " accelstep = " << accelstep << ", range is " << minacc << " to " << maxacc << "\n";
      accelnum = ceil((maxacc-minacc)/accelstep);
      if(usemiddle==1 && accelnum>1 && accelnum%2==0) {
	// We are currently set to probe an even number of acceleration
	// points. This means we won't probe the center of the range
	// of possible accelerations, which is often the most probable.
	// Add one additional point, to ensure that the midpoint of
	// the possible acceleration range gets probed.
	accelnum+=1;
      }
      if(accelnum<=1) {
	// Acceleration step is so large, we'll only use average acceleration
	outstream1 << dist << " " << velAU << " " << (maxacc+minacc)/2.0l/g0 << "\n";
	if(verbose>=1) cout << "Writing point " << stepct << ": " << dist << " " << velAU << " " << (maxacc+minacc)/2.0l/g0 << "\n";
	stepct++;
      } else {
	// Customize the acceleration sampling to have approximately the nominal
	// step size, but optimally sample the full range with the specified number of steps
	for(accelct=0;accelct<accelnum;accelct++) {
	  accelk = minacc + (double(accelct)+0.5l)*(maxacc-minacc)/double(accelnum);
	  accelnorm = accelk/g0;
	  outstream1 << dist << " " << velAU << " " << accelnorm << "\n";
	  if(verbose>=1) cout << "Writing point " << stepct << ": " << dist << " " << velAU << " " << accelnorm << "\n";
	  stepct++;
	}
      }
    } else {
      // Customize the velocity sampling to have approximately the nominal
      // step size, but optimally sample the full range with the specified number of steps.
      // Now for the loop, probing everything except the endpoints
      for(velct=0;velct<velnum;velct++) {
	velAU = -vesc + (double(velct)+0.5)*2.0l*vesc/double(velnum);
	velkm = velAU*AU_KM/SOLARDAY;
	xi = velAU/vesc;
	maxacc = g0;
	minacc = g0*(2.0l*xi*xi - 1.0l);
	if(verbose>=1) cout << "vel = " << velAU << " accelstep = " << accelstep << ", range is " << minacc << " to " << maxacc << "\n";
	accelnum = ceil((maxacc-minacc)/accelstep);
	if(usemiddle==1 && accelnum>1 && accelnum%2==0) {
	  // We are currently set to probe an even number of acceleration
	  // points. This means we won't probe the center of the range
	  // of possible accelerations, which is often the most probable.
	  // Add one additional point, to ensure that the midpoint of
	  // the possible acceleration range gets probed.
	  accelnum+=1;
	}
	if(accelnum<=1) {
	  // Acceleration step is so large, we'll only use average acceleration
	  outstream1 << dist << " " << velAU << " " << (maxacc+minacc)/2.0l/g0 << "\n";
	  if(verbose>=1) cout << "Writing point " << stepct << ": " << dist << " " << velAU << " " << (maxacc+minacc)/2.0l/g0 << "\n";
	  stepct++;
	} else {
	  // Customize the acceleration sampling to have approximately the nominal
	  // step size, but optimally sample the full range with the specified number of steps
	  for(accelct=0;accelct<accelnum;accelct++) {
	    accelk = minacc + (double(accelct)+0.5l)*(maxacc-minacc)/double(accelnum);
	    accelnorm = accelk/g0;
	    outstream1 << dist << " " << velAU << " " << accelnorm << "\n";
	    if(verbose>=1) cout << "Writing point " << stepct << ": " << dist << " " << velAU << " " << accelnorm << "\n";
	    stepct++;
	  }
	}
      }
    }
    dist += diststep;
  }
	
  return(0);
}
