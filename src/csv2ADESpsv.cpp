// September 19, 2026: csv2ADESpsv:
// Given an input CSV file, a column formatting file, and a file
// listing the correct header filenames to use for different ObsCodes,
// produce an ADES psv file for submission to the MPC.

#include "solarsyst_dyn_geo01.h"
#include "cmath"

static void show_usage()
{
  cerr << "Usage: csv2ADESpsv -csvfile csvfile -colformat column format file -head headerlist -timeshift timeshift -trkdigits trkdigits -obsSubID_digits obsSubID_digits -magprec mag_precision -headerlines headerlines -trkSub trkSub -obsSubID obsSubID -band band -rmsRA rmsRA -rmsDec rmsDec -rmsCorr rmsCorr -rmsMag rmsMag -logSNR logSNR -ObsCode ObsCode -permID permID -provID provID -artSat artSat -uncTime uncTime -rmsTime rmsTime -astCat astCat -photCat photCat -mode mode -force_rmsDec -force_rmsRA -force_rmsMag -force_rmsCorr -out outfile\n";
}
#define COLS_TO_READ 6

int main(int argc, char *argv[])
{
  vector <ADESdet> detvec = {};
  vector <ADESdet> onecodevec = {};
  ADESdet onedet;
  string csvfile,colformatfile,headerlist,outfile,headerfile;
  string lnfromfile,stest;
  vector <string> header_ObsCodes;
  vector <string> header_filenames;
  vector <string> svectemp;
  vector <string> det_ObsCodes;
  long header_filenum;
  long i=0;
  long j=0;
  long lct=0;
  char c='0';
  int reachedeof = 0;
  long detnum=0;
  long detct=0;
  int year,month;
  double day;
  string signstring;
  long double timeshift = 0.0l;
  ofstream outstream1;
  ifstream instream1;
  int colreadct=0;
  int verbose=0;
  int magprec=3;
  int force_rmsRA,force_rmsDec,force_rmsMag,force_rmsCorr;
  force_rmsRA = force_rmsDec = force_rmsMag = force_rmsCorr = 0;
  string trkSub,obsSubID,mode,ObsCode,astCat,photCat,band,permID,provID,artSat;
  double rmsRA = -1.0;
  double rmsDec = -1.0;
  double rmsCorr = -999.9;
  double rmsMag = -1.0;
  double logSNR = -1.0;
  double uncTime = -1.0;
  double rmsTime = -1.0;
  mode="CCD";
  int mode_default=1;
  int trkSub_col,obsSubID_col,mode_col,ObsCode_col,MJD_col,RA_col,Dec_col,rmsRA_col,rmsDec_col,rmsCorr_col,astCat_col,photCat_col,mag_col,rmsMag_col,band_col,logSNR_col,uncTime_col,rmsTime_col,permID_col,provID_col,artSat_col;
  trkSub_col = obsSubID_col = mode_col = ObsCode_col = MJD_col = RA_col = Dec_col = rmsRA_col = rmsDec_col = rmsCorr_col = astCat_col = photCat_col = mag_col = rmsMag_col = band_col = logSNR_col = uncTime_col = rmsTime_col = permID_col = provID_col = artSat_col = 0;
  int mjdread,raread,decread,magread;
  mjdread = raread = decread = magread = 1;
  long ObsCodenum=0;
  long ObsCodect=0;
  int trkdigits=0;
  int obsSubID_digits=3;
  long trkct=-1;
  long allobsct=-1;
  long trkstart=0;
  int leadzero = 0;
  int headerlines=0;
  
  if(argc<9) {
    show_usage();
    return(1);
  }

  i=1;
  while(i<argc) {
    cout << "Checking out argv[" << i << "] = " << argv[i] << ".\n";
    if(string(argv[i]) == "-csvfile" || string(argv[i]) == "-infile" || string(argv[i]) == "-csv"  || string(argv[i]) == "-csvdet" || string(argv[i]) == "-csvdetfile") {
      if(i+1 < argc) {
	//There is still something to read;
	csvfile=argv[++i];
	i++;
      }
      else {
	cerr << "Input hldet file keyword supplied with no corresponding argument\n";
	show_usage();
	return(1);
      }
    } else if(string(argv[i]) == "-colformat" || string(argv[i]) == "-format"  || string(argv[i]) == "-col" || string(argv[i]) == "-cf" || string(argv[i]) == "-colfmt" || string(argv[i]) == "--colformat" || string(argv[i]) == "--columnformat" || string(argv[i]) == "--cformat") {
      if(i+1 < argc) {
	//There is still something to read;
	colformatfile=argv[++i];
	i++;
      }
      else {
	cerr << "Column format file keyword supplied with no corresponding argument\n";
	show_usage();
	return(1);
      }
    } else if(string(argv[i]) == "-head" || string(argv[i]) == "-header" || string(argv[i]) == "-headfile" || string(argv[i]) == "-headerfile" || string(argv[i]) == "-headerlist" || string(argv[i]) == "-headlist" || string(argv[i]) == "--headfile" || string(argv[i]) == "--headerfile") {
      if(i+1 < argc) {
	//There is still something to read;
	headerlist=argv[++i];
	i++;
      }
      else {
	cerr << "Input header file keyword supplied with no corresponding argument\n";
	show_usage();
	return(1);
      }
    } else if(string(argv[i]) == "-timeshift" || string(argv[i]) == "-time" || string(argv[i]) == "--timeshift" || string(argv[i]) == "-timediff") {
      if(i+1 < argc) {
	//There is still something to read;
	timeshift = stod(argv[++i]);
	i++;
      }
      else {
	cerr << "Output file keyword supplied with no corresponding argument\n";
	show_usage();
	return(1);
      }
    } else if(string(argv[i]) == "-trkdigits") {
      if(i+1 < argc) {
	//There is still something to read;
	trkdigits = stoi(argv[++i]);
	i++;
      }
      else {
	cerr << "Input trkdigits keyword supplied with no corresponding argument\n";
	show_usage();
	return(1);
      }
    } else if(string(argv[i]) == "-obsSubID_digits") {
      if(i+1 < argc) {
	//There is still something to read;
	obsSubID_digits = stoi(argv[++i]);
 	i++;
      }
      else {
	cerr << "Input obsSubID_digits keyword supplied with no corresponding argument\n";
	show_usage();
	return(1);
      }
    } else if(string(argv[i]) == "-magprec") {
      if(i+1 < argc) {
	//There is still something to read;
	magprec = stoi(argv[++i]);
	i++;
      }
      else {
	cerr << "Input magnitude precision keyword supplied with no corresponding argument\n";
	show_usage();
	return(1);
      }
    } else if(string(argv[i]) == "-headerlines") {
      if(i+1 < argc) {
	//There is still something to read;
	headerlines = stoi(argv[++i]);
	i++;
      }
      else {
	cerr << "Input header lines keyword supplied with no corresponding argument\n";
	show_usage();
	return(1);
      }
    } else if(string(argv[i]) == "-trkSub" || string(argv[i]) == "-trksub") {
      if(i+1 < argc) {
	//There is still something to read;
	trkSub=argv[++i];
	i++;
      }
      else {
	cerr << "Input trkSub keyword supplied with no corresponding argument\n";
	show_usage();
	return(1);
      }
    } else if(string(argv[i]) == "-obsSubID") {
      if(i+1 < argc) {
	//There is still something to read;
	obsSubID=argv[++i];
	i++;
      }
      else {
	cerr << "Input obsSubID keyword supplied with no corresponding argument\n";
	show_usage();
	return(1);
      }
    } else if(string(argv[i]) == "-mode") {
      if(i+1 < argc) {
	//There is still something to read;
	mode=argv[++i];
	i++;
	mode_default=0;
      }
      else {
	cerr << "Input mode keyword supplied with no corresponding argument\n";
	show_usage();
	return(1);
      }
    } else if(string(argv[i]) == "-ObsCode") {
      if(i+1 < argc) {
	//There is still something to read;
	ObsCode=argv[++i];
	i++;
      }
      else {
	cerr << "Input ObsCode keyword supplied with no corresponding argument\n";
	show_usage();
	return(1);
      }
    } else if(string(argv[i]) == "-astCat") {
      if(i+1 < argc) {
	//There is still something to read;
	astCat=argv[++i];
	i++;
      }
      else {
	cerr << "Input astCat keyword supplied with no corresponding argument\n";
	show_usage();
	return(1);
      }
    } else if(string(argv[i]) == "-photCat") {
      if(i+1 < argc) {
	//There is still something to read;
	photCat=argv[++i];
	i++;
      }
      else {
	cerr << "Input photCat keyword supplied with no corresponding argument\n";
	show_usage();
	return(1);
      }
    } else if(string(argv[i]) == "-band") {
      if(i+1 < argc) {
	//There is still something to read;
	band=argv[++i];
	i++;
      }
      else {
	cerr << "Input band keyword supplied with no corresponding argument\n";
	show_usage();
	return(1);
      }
    } else if(string(argv[i]) == "-rmsRA") {
      if(i+1 < argc) {
	//There is still something to read;
	rmsRA=stod(argv[++i]);
	i++;
      }
      else {
	cerr << "Input rmsRA keyword supplied with no corresponding argument\n";
	show_usage();
	return(1);
      }
    } else if(string(argv[i]) == "-rmsDec") {
      if(i+1 < argc) {
	//There is still something to read;
	rmsDec=stod(argv[++i]);
	i++;
      }
      else {
	cerr << "Input rmsDec keyword supplied with no corresponding argument\n";
	show_usage();
	return(1);
      }
    } else if(string(argv[i]) == "-rmsCorr") {
      if(i+1 < argc) {
	//There is still something to read;
	rmsCorr=stod(argv[++i]);
	i++;
      }
      else {
	cerr << "Input rmsCorr keyword supplied with no corresponding argument\n";
	show_usage();
	return(1);
      }
    } else if(string(argv[i]) == "-rmsMag") {
      if(i+1 < argc) {
	//There is still something to read;
	rmsMag=stod(argv[++i]);
	i++;
      }
      else {
	cerr << "Input rmsMag keyword supplied with no corresponding argument\n";
	show_usage();
	return(1);
      }
    } else if(string(argv[i]) == "-logSNR") {
      if(i+1 < argc) {
	//There is still something to read;
	logSNR=stod(argv[++i]);
	i++;
      }
      else {
	cerr << "Input logSNR keyword supplied with no corresponding argument\n";
	show_usage();
	return(1);
      }
    } else if(string(argv[i]) == "-uncTime") {
      if(i+1 < argc) {
	//There is still something to read;
	uncTime=stod(argv[++i]);
	i++;
      }
      else {
	cerr << "Input uncTime keyword supplied with no corresponding argument\n";
	show_usage();
	return(1);
      }
    } else if(string(argv[i]) == "-rmsTime") {
      if(i+1 < argc) {
	//There is still something to read;
	rmsTime=stod(argv[++i]);
	i++;
      }
      else {
	cerr << "Input rmsTime keyword supplied with no corresponding argument\n";
	show_usage();
	return(1);
      }
    } else if(string(argv[i]) == "-permID" || string(argv[i]) == "-permid") {
      if(i+1 < argc) {
	//There is still something to read;
	permID=argv[++i];
	i++;
      }
      else {
	cerr << "Input permID keyword supplied with no corresponding argument\n";
	show_usage();
	return(1);
      }
    } else if(string(argv[i]) == "-provID" || string(argv[i]) == "-provid") {
      if(i+1 < argc) {
	//There is still something to read;
	provID=argv[++i];
	i++;
      }
      else {
	cerr << "Input provID keyword supplied with no corresponding argument\n";
	show_usage();
	return(1);
      }
    } else if(string(argv[i]) == "-artSat" || string(argv[i]) == "-artsat") {
      if(i+1 < argc) {
	//There is still something to read;
	artSat=argv[++i];
	i++;
      }
      else {
	cerr << "Input artSat keyword supplied with no corresponding argument\n";
	show_usage();
	return(1);
      }
    } else if(string(argv[i]) == "-force_rmsRA") {
      if(i+1 < argc) {
	force_rmsRA=1;
	force_rmsDec=1;
	i++;
      }
    } else if(string(argv[i]) == "-force_rmsDec") {
      if(i+1 < argc) {
	force_rmsRA=1;
	force_rmsDec=1;
	i++;
      }
    } else if(string(argv[i]) == "-force_rmsMag") {
      if(i+1 < argc) {
	force_rmsMag=1;
	i++;
      }
    } else if(string(argv[i]) == "-force_rmsCorr") {
      if(i+1 < argc) {
	force_rmsCorr=1;
	force_rmsRA=1;
	force_rmsDec=1;
	i++;
      }
    } else if(string(argv[i]) == "-out" || string(argv[i]) == "-outfile" || string(argv[i]) == "--outfile" || string(argv[i]) == "--output_ADES_file") {
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
    } else if(string(argv[i]) == "-verbose" || string(argv[i]) == "-verb" || string(argv[i]) == "-VERBOSE" || string(argv[i]) == "-VERB" || string(argv[i]) == "--verbose" || string(argv[i]) == "--VERBOSE" || string(argv[i]) == "--VERB") {
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

  // Catch required parameters if missing
  if(csvfile.size()<=0) {
    cerr << "\nERROR: input csv file is required\n";
    show_usage();
    return(1);
  } else if(outfile.size()<=0) {
    cerr << "\nERROR: output filename is required\n";
    show_usage();
    return(1);
  }

  cout << "input csv file " << csvfile << "\n";
  cout << "time shift " << timeshift << " seconds\n";
  cout << "output file " << outfile << "\n";

  // Read the column formatting file, if any
  if(colformatfile.size()>0) {
    instream1.open(colformatfile);
    if(!instream1)  {
      cerr << "ERROR: unable to open input file " << colformatfile << "\n";
      return(1);
    }
    colreadct=0;
    while(!instream1.eof() && !instream1.fail() && !instream1.bad()) {
      instream1 >> stest;
      if(stest == "trkSub_col") {
	instream1 >> trkSub_col;
	if(!instream1.eof() && !instream1.fail() && !instream1.bad()) colreadct++;
      } else if(stest == "obsSubID_col") {
	instream1 >> obsSubID_col;
	if(!instream1.eof() && !instream1.fail() && !instream1.bad()) colreadct++;
      } else if(stest == "mode_col") {
	instream1 >> mode_col;
	if(!instream1.eof() && !instream1.fail() && !instream1.bad()) colreadct++;
      } else if(stest == "ObsCode_col") {
	instream1 >> ObsCode_col;
	if(!instream1.eof() && !instream1.fail() && !instream1.bad()) colreadct++;
      } else if(stest == "MJD_col") {
	instream1 >> MJD_col;
	if(!instream1.eof() && !instream1.fail() && !instream1.bad()) colreadct++;
      } else if(stest == "RA_col") {
	instream1 >> RA_col;
	if(!instream1.eof() && !instream1.fail() && !instream1.bad()) colreadct++;
      } else if(stest == "Dec_col") {
	instream1 >> Dec_col;
	if(!instream1.eof() && !instream1.fail() && !instream1.bad()) colreadct++;
      } else if(stest == "rmsRA_col") {
	instream1 >> rmsRA_col;
	if(!instream1.eof() && !instream1.fail() && !instream1.bad()) colreadct++;
      } else if(stest == "rmsDec_col") {
	instream1 >> rmsDec_col;
	if(!instream1.eof() && !instream1.fail() && !instream1.bad()) colreadct++;
      } else if (stest == "rmsCorr_col") {
	instream1 >> rmsCorr_col;
	if(!instream1.eof() && !instream1.fail() && !instream1.bad()) colreadct++;
      } else if(stest == "astCat_col") {
	instream1 >> astCat_col;
	if(!instream1.eof() && !instream1.fail() && !instream1.bad()) colreadct++;
      } else if(stest == "photCat_col") {
	instream1 >> photCat_col;
	if(!instream1.eof() && !instream1.fail() && !instream1.bad()) colreadct++;
      } else if(stest == "mag_col") {
	instream1 >> mag_col;
	if(!instream1.eof() && !instream1.fail() && !instream1.bad()) colreadct++;
      } else if(stest == "rmsMag_col") {
	instream1 >> rmsMag_col;
	if(!instream1.eof() && !instream1.fail() && !instream1.bad()) colreadct++;
      } else if(stest == "band_col") {
	instream1 >> band_col;
	if(!instream1.eof() && !instream1.fail() && !instream1.bad()) colreadct++;
      } else if(stest == "logSNR_col") {
	instream1 >> logSNR_col;
	if(!instream1.eof() && !instream1.fail() && !instream1.bad()) colreadct++;
      } else if(stest == "uncTime_col") {
	instream1 >> uncTime_col;
	if(!instream1.eof() && !instream1.fail() && !instream1.bad()) colreadct++;
      } else if(stest == "rmsTime_col") {
	instream1 >> rmsTime_col;
	if(!instream1.eof() && !instream1.fail() && !instream1.bad()) colreadct++;
      } else if(stest == "permID_col") {
	instream1 >> permID_col;
	if(!instream1.eof() && !instream1.fail() && !instream1.bad()) colreadct++;
      } else if(stest == "provID_col") {
	instream1 >> provID_col;
	if(!instream1.eof() && !instream1.fail() && !instream1.bad()) colreadct++;
      } else if(stest == "artSat_col") {
	instream1 >> artSat_col;
	if(!instream1.eof() && !instream1.fail() && !instream1.bad()) colreadct++;
      } else {
	cout << "WARNING: unrecognized string " << stest << " read from column formatting file\n";
      }
    }
    instream1.close();
    if(colreadct<COLS_TO_READ) {
      cerr << "ERROR: only " << colreadct << " column specifications, of " << COLS_TO_READ << " required, were read from column format file " << colformatfile << ".\n";
      return(2);
      // Required columns: MJD, RA, Dec, ObsCode, mag, band.
    }
  } else {
    cerr << "ERROR: no column formatting file supplied\n";
    show_usage();
    return(3);
  }
  cout << colreadct << " columns read from column formatting file:\n";
  if(trkSub_col!=0) cout << "trkSub_col = " << trkSub_col << "\n";
  if(obsSubID_col!=0) cout << "obsSubID_col = " << obsSubID_col << "\n";
  if(mode_col!=0) cout << "mode_col = " << mode_col << "\n";
  if(ObsCode_col!=0) cout << "ObsCode_col = " << ObsCode_col << "\n";
  if(MJD_col!=0) cout << "MJD_col = " << MJD_col << "\n";
  if(RA_col!=0) cout << "RA_col = " << RA_col << "\n";
  if(Dec_col!=0) cout << "Dec_col = " << Dec_col << "\n";
  if(rmsRA_col!=0) cout << "rmsRA_col = " << rmsRA_col << "\n";
  if(rmsDec_col!=0) cout << "rmsDec_col = " << rmsDec_col << "\n";
  if(rmsCorr_col!=0) cout << "rmsCorr_col = " << rmsCorr_col << "\n";
  if(astCat_col!=0) cout << "astCat_col = " << astCat_col << "\n";
  if(photCat_col!=0) cout << "photCat_col = " << photCat_col << "\n";
  if(mag_col!=0) cout << "mag_col = " << mag_col << "\n";
  if(rmsMag_col!=0) cout << "rmsMag_col = " << rmsMag_col << "\n";
  if(band_col!=0) cout << "band_col = " << band_col << "\n";
  if(logSNR_col!=0) cout << "logSNR_col = " << logSNR_col << "\n";
  if(uncTime_col!=0) cout << "uncTime_col = " << uncTime_col << "\n";
  if(rmsTime_col!=0) cout << "rmsTime_col = " << rmsTime_col << "\n";
  if(permID_col!=0) cout << "permID_col = " << permID_col << "\n";
  if(provID_col!=0) cout << "provID_col = " << provID_col << "\n";
  if(artSat_col!=0) cout << "artSat_col = " << artSat_col << "\n";

  if(verbose) cout << "verbose mode\n";
  
 // Read input observation file.
  detvec={};
  instream1.open(csvfile);
  if(!instream1) {
    cerr << "can't open input file " << csvfile << "\n";
    return(1);
  }
  reachedeof = 0;
  for(i=0;i<headerlines;i++) getline(instream1,lnfromfile); // Skip header lines, if any.
  
  while(reachedeof==0) {
    getline(instream1,lnfromfile);
    lct++;
    if(!instream1.eof() && !instream1.fail() && !instream1.bad()) ; // Read on.
    else if(instream1.eof()) reachedeof=1; //End of file, fine.
    else if(instream1.fail()) reachedeof=-1; //Something wrong, warn
    else if(instream1.bad()) reachedeof=-2; //Worse problem, warn
    i=0;
    j = 0;
    c='0';
    mjdread = raread = decread = magread = 1;

    if(trkSub.size()>0) onedet.trkSub = trkSub;
    if(obsSubID.size()>0) onedet.obsSubID = obsSubID;
    if(mode_default==0) onedet.mode = mode;
    if(ObsCode.size()>0) onedet.ObsCode = ObsCode;
    if(astCat.size()>0) onedet.astCat = astCat;
    if(photCat.size()>0) onedet.photCat = photCat;
    if(band.size()>0) onedet.band = band;
    if(rmsRA>0.0) onedet.rmsRA = rmsRA;
    if(rmsDec>0.0) onedet.rmsDec = rmsDec;
    if(rmsCorr>=-1.0 && rmsCorr<=1.0) onedet.rmsCorr = rmsCorr;
    if(rmsMag>0.0) onedet.rmsMag = rmsMag;
    if(logSNR>0.0) onedet.logSNR = logSNR;
    if(uncTime>0.0) onedet.uncTime = uncTime;
    if(rmsTime>0.0) onedet.rmsTime = rmsTime;
    if(permID.size()>0) onedet.permID = permID;
    if(provID.size()>0) onedet.provID = provID;
    if(artSat.size()>0) onedet.artSat = artSat;
  
    while(i<long(lnfromfile.size()) && lnfromfile.size()>=30 && reachedeof == 0) {
      // Note check on line length: it is completely impossible for a
      // line containing all the required quantities at minimum plausible
      // precision to be less than 30 characters long.
      c='0';
      stest="";
      while(i<long(lnfromfile.size()) && c!=',' && c!='\n' && c!=EOF) {
	c=lnfromfile[i];
	if(c!=',' && c!='\n' && c!=EOF) stest.push_back(c);
	i++;
      }
      // We just finished reading something
      j++;
      if(j==trkSub_col) {
	onedet.trkSub = stest;
      } else if(j==obsSubID_col) {
	onedet.obsSubID = stest;
      } else if(j==mode_col) {
	onedet.mode = stest;
      } else if(j==ObsCode_col) {
	onedet.ObsCode = stest;
      } else if(j==MJD_col) {
	mjdread=1;
	try { onedet.MJD = stod(stest); }
	catch(...) { cerr << "ERROR: cannot read MJD string " << stest << " from line " << lnfromfile << "\n"; 
	mjdread=0; }
      } else if(j==RA_col) {
	raread=1;
	try { onedet.RA = stod(stest); }
	catch(...) { cerr << "ERROR: cannot read RA string " << stest << " from line " << lnfromfile << "\n"; 
	raread=0; }
      } else if(j==Dec_col) {
	decread=1;
	try { onedet.Dec = stod(stest); }
	catch(...) { cerr << "ERROR: cannot read Dec string " << stest << " from line " << lnfromfile << "\n"; 
	decread=0; }
      } else if(j==rmsRA_col) {
	try { onedet.rmsRA = stod(stest); }
	catch(...) { cerr << "ERROR: cannot read rmsRA string " << stest << " from line " << lnfromfile << "\n"; }
	cout << "rmsRA string " << stest << " = " << onedet.rmsRA << "\n";
      } else if(j==rmsDec_col) {
	try { onedet.rmsDec = stod(stest); }
	catch(...) { cerr << "ERROR: cannot read rmsDec string " << stest << " from line " << lnfromfile << "\n"; }
      } else if(j==rmsCorr_col) {
	try { onedet.rmsCorr = stod(stest); }
	catch(...) { cerr << "ERROR: cannot read rmsCorr string " << stest << " from line " << lnfromfile << "\n"; }
	cout << "rmsDec string " << stest << " = " << onedet.rmsDec << "\n";
      } else if(j==astCat_col) {
	onedet.astCat = stest;
      } else if(j==photCat_col) {
	onedet.photCat = stest;
      } else if(j==mag_col) {
	magread=1;
	try { onedet.mag=stod(stest); }
	catch(...) { cerr << "ERROR: cannot read mag string " << stest << " from line " << lnfromfile << "\n"; 
	magread=0; }
      } else if(j==band_col) {
	onedet.band = stest;
      } else if(j==rmsMag_col) {
	try { onedet.rmsMag=stod(stest); }
	catch(...) { cerr << "ERROR: cannot read rmsMag string " << stest << " from line " << lnfromfile << "\n"; }
      } else if(j==logSNR_col) {
	try { onedet.logSNR=stod(stest); }
	catch(...) { cerr << "ERROR: cannot read logSNR string " << stest << " from line " << lnfromfile << "\n"; }
      } else if(j==uncTime_col) {
	try { onedet.uncTime=stod(stest); }
	catch(...) { cerr << "ERROR: cannot read uncTime string " << stest << " from line " << lnfromfile << "\n"; }
      } else if(j==rmsTime_col) {
	try { onedet.rmsTime=stod(stest); }
	catch(...) { cerr << "ERROR: cannot read rmsTime string " << stest << " from line " << lnfromfile << "\n"; }
      } else if(j==permID_col) {
	onedet.permID = stest;
      } else if(j==provID_col) {
	onedet.provID = stest;
      } else if(j==artSat_col) {
	onedet.artSat = stest;
      } 
      // cout<<"Column "<< j << " read as " << stest << ".\n";
    }
    if(reachedeof == 0 && lnfromfile.size()>=30) {
      if(!mjdread) {
	cerr << "ERROR: MJD not read from line " << detvec.size()+1 << " of input detection file " << csvfile << "!\n";
	cerr << "Here is the line: " << lnfromfile << "\n";
	return(2);
      }
      if(!raread) {
	cerr << "ERROR: RA not read from line " << detvec.size()+1 << " of input detection file " << csvfile << "!\n";
	cerr << "Here is the line: " << lnfromfile << "\n";
	return(2);
      }
      if(!decread) {
	cerr << "ERROR: Dec not read from line " << detvec.size()+1 << " of input detection file " << csvfile << "!\n";
	cerr << "Here is the line: " << lnfromfile << "\n";
	return(2);
      }
      if(!magread) {
	cerr << "ERROR: mag not read from line " << detvec.size()+1 << " of input detection file " << csvfile << "!\n";
	cerr << "Here is the line: " << lnfromfile << "\n";
	return(2);
      }
      detvec.push_back(onedet);
    }
  }
  instream1.close();
  
  if(reachedeof==1) { 
    cout << "Input file " << csvfile << " read successfully to the end.\n";
  } else if(reachedeof==0) {
    cerr << "ERROR: Stopped reading file " << csvfile << " before the end\n";
    return(1);
  } else if(reachedeof==-1) {
    cerr << "Warning: file read failed\n";
    return(1);
  } else if(reachedeof==-2) {
    cerr << "Warning: file possibly corrupted\n";
    return(2);
  } else return(reachedeof);

  // Done reading input file
  detnum = detvec.size();
  cout << "Read " << detnum << " data lines from observation file " << csvfile << "\n";

  // Read header file. This is a two column file, with ObsCode in column 1
  // and the name of a header file to be used with that ObsCode in column 2.
  // It is fine to have header file lines that don't correspond to any ObsCode
  // that actually occurs in the input data, but it is NOT fine to have ObsCodes
  // in the input data without a corresponding line in the header file.
  
  instream1.open(headerlist);
  if(!instream1) {
    cerr << "can't open input file " << headerlist << "\n";
    return(1);
  }
  reachedeof = 0;
  while(!instream1.eof() && !instream1.fail() && !instream1.bad()) {
    stest="";
    instream1 >> stest;
    if(stest.size()==3) header_ObsCodes.push_back(stest);
    else if(stest.size()>0) {
      cerr << "ERROR: headerlist ObsCode entry is " << stest << ", which does not look like a valid ObsCode\n";
      return(5);
    }
    if(!instream1.eof() && !instream1.fail() && !instream1.bad()) {
      stest="";
      instream1 >> stest;
      if(stest.size()>0) header_filenames.push_back(stest);
    }
  }
  instream1.close();
  header_filenum = header_ObsCodes.size();
  cout << "Read " << header_filenum << " lines from header list file " << headerlist << "\n";
  if(long(header_filenames.size()) != header_filenum) {
    cerr << "ERROR: mismatch between number of header_ObsCodes (" << header_filenum << ") and header filenames (" << header_filenames.size() << ")\n";
    return(6);
  }

  // Create a unique list of ObsCodes used in the input detection file
  svectemp = det_ObsCodes = {};
  for(detct=0;detct<detnum;detct++) {
    svectemp.push_back(detvec[detct].ObsCode);
  }
  sort(svectemp.begin(),svectemp.end());
  det_ObsCodes.push_back(svectemp[0]);
  for(detct=1;detct<detnum;detct++) {
    if(svectemp[detct]!=det_ObsCodes[det_ObsCodes.size()-1]) {
      det_ObsCodes.push_back(svectemp[detct]);
    }
  }
  ObsCodenum=det_ObsCodes.size();
  cout << "Found a total of " << det_ObsCodes.size() << " distinct ObsCodes in the detection file:\n";
  for(i=0;i<long(det_ObsCodes.size());i++) {
    cout << det_ObsCodes[i] << "\n";
  }

   for(i=0;i<long(det_ObsCodes.size());i++) {
     int foundit=0;
     for(j=0;j<header_filenum;j++) {
       if(det_ObsCodes[i]==header_ObsCodes[j]) foundit=1;
     }
     if(foundit==0) {
       cerr << "ERROR: ObsCode " << det_ObsCodes[i] << ", used in the input detection file,\n";
       cerr << "is missing from the input header list file\n";
       return(7);
     }
   }
  outstream1.open(outfile);
  outstream1 << "# version=2022\n";
  for(ObsCodect=0;ObsCodect<ObsCodenum;ObsCodect++) {
    // Create a new vector of type ADESdet with only the detections from this ObsCode.
    onecodevec = {};
    for(detct=0;detct<=detnum;detct++) {
      if(detvec[detct].ObsCode==det_ObsCodes[ObsCodect]) onecodevec.push_back(detvec[detct]);
    }
    // Sort the new vector by MJD
    sort(onecodevec.begin(), onecodevec.end(), early_ADESdet());
    // Find the appropriate header file for this ObsCode
    headerfile="";
    for(j=0;j<header_filenum;j++) {
      if(det_ObsCodes[ObsCodect]==header_ObsCodes[j]) headerfile=header_filenames[j];
    }
    if(headerfile.size()<=0) {
      cerr << "ERROR: unable to find header file for ObsCode " << det_ObsCodes[ObsCodect] << ", though it seemed to be there on the earlier check\n";
      return(7);
    }
    detct=0;
    while(detct<long(onecodevec.size())) {
      if(detct==0 || (detct>0 && (onecodevec[detct].MJD - onecodevec[detct-1].MJD)>NIGHTSTEP)) {
	// We are at the start of a new tracklet.
	trkct++;
	trkstart=detct;
	if(trkct>0) outstream1 << "\n";
	// Write the header to the output file 
	instream1.open(headerfile);
	if(!instream1) {
	  cerr << "can't open input file " << headerfile << "\n";
	  return(1);
	}
	while (!instream1.eof() && !instream1.fail() && !instream1.bad()) {
	  getline(instream1,lnfromfile);
	  if(lnfromfile.size()>0) outstream1 << lnfromfile << "\n";
	}
	instream1.close();
	// Construct the header line appropriate for this tracklet
	if(onecodevec[detct].permID.size()>0) {
	  outstream1 << "permID";
	  for(j=6;j<long(onecodevec[detct].permID.size());j++) outstream1 << " ";
	  outstream1 << "|";
	} else if(onecodevec[detct].provID.size()>0) {
	  outstream1 << "provID";
	  for(j=6;j<long(onecodevec[detct].provID.size());j++) outstream1 << " ";
	  outstream1 << "|";
	} else if(onecodevec[detct].artSat.size()>0) {
	  outstream1 << "artSat";
	  for(j=6;j<long(onecodevec[detct].artSat.size());j++) outstream1 << " ";
	  outstream1 << "|";
	}
	outstream1 << "trkSub ";
	for(j=7;j<long(onecodevec[detct].trkSub.size() + trkdigits);j++) outstream1 << " ";
	outstream1 << "|obsSubID";
	for(j=8;j<long(onecodevec[detct].obsSubID.size() + obsSubID_digits);j++) outstream1 << " ";
	outstream1 << "|mode|stn |obsTime                 |ra         |dec        ";
	if((onecodevec[detct].rmsRA>0.0 && onecodevec[detct].rmsDec>0.0) || force_rmsRA || force_rmsDec) outstream1 << "|rmsRA |rmsDec";
	if((onecodevec[detct].rmsCorr>=-1.0 && onecodevec[detct].rmsCorr<=1.0) | force_rmsCorr) outstream1 << "|rmsCorr";
	outstream1 << "|mag";
	for(j=3;j<3+magprec;j++) outstream1 << " ";
	if(onecodevec[detct].rmsMag>0.0 || force_rmsMag) outstream1 << "|rmsMag";
	outstream1 << "|band";
	if(onecodevec[detct].astCat.size()>0) outstream1 << "|astCat  ";
	if(onecodevec[detct].photCat.size()>0) outstream1 << "|photCat ";
	if(onecodevec[detct].logSNR>0.0) outstream1 << "|logSNR";
	if(onecodevec[detct].uncTime>0.0) outstream1 << "|uncTime";
	if(onecodevec[detct].rmsTime>0.0) outstream1 << "|rmsTime";
        outstream1 << "\n";
	
	// Done with column headers, write out actual quantities.
	allobsct++;
	if(onecodevec[trkstart].permID.size()>0 && onecodevec[detct].permID.size()>0) {
	  outstream1 << onecodevec[detct].permID << "|";
	} else if(onecodevec[trkstart].permID.size()>0 && onecodevec[detct].permID.size()==0) {
	  for(j=0;j<long(onecodevec[trkstart].permID.size());j++) outstream1 << " ";
	  outstream1 << "|";
	} else if(onecodevec[trkstart].provID.size()>0 && onecodevec[detct].provID.size()>0) {
	  outstream1 << onecodevec[detct].provID << "|";
	} else if(onecodevec[trkstart].provID.size()>0 && onecodevec[detct].provID.size()==0) {
	  for(j=0;j<long(onecodevec[trkstart].provID.size());j++) outstream1 << " ";
	  outstream1 << "|";
	} else if(onecodevec[trkstart].artSat.size()>0 && onecodevec[detct].artSat.size()>0) {
	  outstream1 << onecodevec[detct].artSat << "|";
	} else if(onecodevec[trkstart].artSat.size()>0 && onecodevec[detct].artSat.size()==0) {
	  for(j=0;j<long(onecodevec[trkstart].artSat.size());j++) outstream1 << " ";
	  outstream1 << "|";
	}
	outstream1 << onecodevec[detct].trkSub;
	for(j=onecodevec[detct].trkSub.size();j<7;j++) outstream1 << " ";
	if(trkdigits!=0) {
	  leadzero = 0;
	  if(trkct>0) leadzero = trkdigits - floor(log10(trkct)) - 1;
	  if(leadzero<0) leadzero=0;
	  for(j=0;j<leadzero;j++) outstream1 << "0";
	  outstream1 << trkct;
	}
	outstream1 << "|" << onecodevec[detct].obsSubID;
	if(obsSubID_digits!=0) {
	  outstream1 << ".";
	  leadzero = obsSubID_digits-1;
	  if(allobsct>0) leadzero = obsSubID_digits - floor(log10(allobsct)) - 1;
	  if(leadzero<0) leadzero=0;
	  for(j=0;j<leadzero;j++) outstream1 << "0";
	  outstream1 << allobsct;
	}
	outstream1 << "|" << onecodevec[detct].mode << " ";
	outstream1 << "|" << onecodevec[detct].ObsCode;
	for(j=long(onecodevec[detct].ObsCode.size());j<4;j++) outstream1 << " ";
	// Adjust the input MJD by timeshift, which could be positive (e.g., offset to mid-exposure)
	// or negative (e.g., to convert from TAI to UTC).
	mjd2mpcdate(onecodevec[detct].MJD +(timeshift/SOLARDAY), year, month, day);
	long iday = floor(day);
	double fracday = day - double(iday);
	double hour = floor(24.0*fracday);
	double minute = floor(1440.0*fracday - 60.0*hour);
	double second = fracday*SOLARDAY - 3600.0*hour - 60.0*minute;
	outstream1 <<  "|" << year << "-";
	if(month<10) outstream1 << "0";
	outstream1 << month << "-";
	if(iday<10) outstream1 << "0";
	outstream1 << iday << "T";
	if(hour<10.0) outstream1 << "0";
	outstream1 << fixed << setprecision(0) << hour << ":";
	if(minute<10.0) outstream1 << "0";
	outstream1 << fixed << setprecision(0) << minute << ":";
	if(second<9.9995) outstream1 << "0";
	outstream1 << fixed << setprecision(3) << second << "Z";
	
	if(onecodevec[detct].RA>=99.99999995) outstream1  << fixed << setprecision(7) << "|" << onecodevec[detct].RA << "";
	else if(onecodevec[detct].RA>=9.99999995) outstream1  << fixed << setprecision(7) << "| " << onecodevec[detct].RA << "";
	else outstream1  << fixed << setprecision(6) << "|  " << onecodevec[detct].RA;
	
	if(onecodevec[detct].Dec<0.0 && fabs(onecodevec[detct].Dec)>=9.99999995) outstream1  << fixed << setprecision(7) << "|" << onecodevec[detct].Dec << "";
	else if(onecodevec[detct].Dec<0.0) outstream1  << fixed << setprecision(7) << "| " << onecodevec[detct].Dec;
	else if(onecodevec[detct].Dec>=9.99999995) outstream1  << fixed << setprecision(7) << "| " << onecodevec[detct].Dec;
	else outstream1  << fixed << setprecision(7) << "|  " << onecodevec[detct].Dec << "";
      
	if((onecodevec[trkstart].rmsRA>0.0 || force_rmsRA) && onecodevec[detct].rmsRA>0.0) {
	  outstream1 << fixed << setprecision(4) << "|" << onecodevec[detct].rmsRA;
	} else if(onecodevec[trkstart].rmsRA>0.0 || force_rmsRA) {
	  outstream1 << "|      ";
	}
	
	if((onecodevec[trkstart].rmsDec>0.0 || force_rmsDec) && onecodevec[detct].rmsDec>0.0) {
	  outstream1 << fixed << setprecision(4)  << "|" << onecodevec[detct].rmsDec;
	} else if(onecodevec[trkstart].rmsDec>0.0 || force_rmsDec) {
	  outstream1 << "|      ";
	}
	
	if(((fabs(onecodevec[trkstart].rmsCorr)<=1.0 && onecodevec[trkstart].rmsRA>0.0 && onecodevec[trkstart].rmsDec>0.0) || force_rmsCorr) && fabs(onecodevec[detct].rmsCorr)<=1.0 && onecodevec[detct].rmsRA>0.0 && onecodevec[detct].rmsDec>0.0) {
	  if(onecodevec[detct].rmsCorr<0.0) outstream1 << fixed << setprecision(4) << "|" << onecodevec[detct].rmsCorr;
	  else outstream1 << fixed << setprecision(4) << "| " << onecodevec[detct].rmsCorr;
	} else if((fabs(onecodevec[trkstart].rmsCorr)<=1.0 && onecodevec[trkstart].rmsRA>0.0 && onecodevec[trkstart].rmsDec>0.0) || force_rmsCorr) {
	  outstream1 << "|       ";
	} 

	if(magprec<=0) {
	  if(onecodevec[detct].mag>=9.5) outstream1  << fixed << setprecision(magprec)  << "|" << " " << onecodevec[detct].mag;
	  else outstream1  << fixed << setprecision(magprec)  << "|" << "  " << onecodevec[detct].mag;
	} else if(magprec==1) {
	  if(onecodevec[detct].mag>=9.95) outstream1  << fixed << setprecision(magprec)  << "|"  << onecodevec[detct].mag;
	  else outstream1  << fixed << setprecision(magprec)  << "|" << " " << onecodevec[detct].mag << "|";
	} else if(magprec==2) {
	  if(onecodevec[detct].mag>=9.995) outstream1  << fixed << setprecision(magprec)  << "|"<< onecodevec[detct].mag;
	  else outstream1  << fixed << setprecision(magprec)  << "|" << " " << onecodevec[detct].mag << "|";
	} else if(magprec==3) {
	  if(onecodevec[detct].mag>=9.9995) outstream1  << fixed << setprecision(magprec)  << "|" << onecodevec[detct].mag;
	  else outstream1  << fixed << setprecision(magprec)  << "|" << " " << onecodevec[detct].mag << "|";
	} else {
	  if(onecodevec[detct].mag>=9.99995) outstream1  << fixed << setprecision(magprec)  << "|" << onecodevec[detct].mag;
	  else outstream1  << fixed << setprecision(magprec)  << "|" << " " << onecodevec[detct].mag;
	}
	
	if((onecodevec[trkstart].rmsMag>0.0 || force_rmsMag) && onecodevec[detct].rmsMag>0.0) {
	  if(magprec<=0) outstream1  << fixed << setprecision(magprec)  << "|" << onecodevec[detct].rmsMag << "     ";
	  else if(magprec==1) outstream1  << fixed << setprecision(magprec) << "|" << onecodevec[detct].rmsMag << "   ";
	  else if(magprec==2) outstream1  << fixed << setprecision(magprec) << "|" << onecodevec[detct].rmsMag << "  ";
	  else if(magprec==3) outstream1  << fixed << setprecision(magprec) << "|" << onecodevec[detct].rmsMag << " ";
	  else outstream1  << fixed << setprecision(magprec) << "|" << onecodevec[detct].rmsMag << "";
	} else if(onecodevec[trkstart].rmsMag>0.0 || force_rmsMag) {
	  outstream1 << "|      ";
	}
	
	outstream1 << "|" << onecodevec[detct].band;
	for(j=onecodevec[detct].band.size();j<4;j++) outstream1 << " ";
	
	if(onecodevec[trkstart].astCat.size()>0 && onecodevec[detct].astCat.size()>0) {
	  outstream1 << "|" << onecodevec[detct].astCat;
	  for(j=onecodevec[detct].astCat.size();j<8;j++) outstream1 << " ";
	} else if(onecodevec[trkstart].astCat.size()>0 && onecodevec[detct].astCat.size()<=0) {
	  outstream1 << "|        ";
	}

	if(onecodevec[trkstart].photCat.size()>0 && onecodevec[detct].photCat.size()>0) {
	  outstream1 << "|" << onecodevec[detct].photCat;
	  for(j=onecodevec[detct].photCat.size();j<8;j++) outstream1 << " ";
	} else if(onecodevec[trkstart].photCat.size()>0 && onecodevec[detct].photCat.size()<=0) {
	  outstream1 << "|        ";
	}
	
	if(onecodevec[trkstart].logSNR>0.0 && onecodevec[detct].logSNR>0.0) {
	  outstream1 << fixed << setprecision(4) << "|" << onecodevec[detct].logSNR << " ";
	} else if(onecodevec[trkstart].logSNR>0.0 && onecodevec[detct].logSNR<=0.0) {
	  outstream1 << "|      ";
	}

	if(onecodevec[trkstart].uncTime>0.0 && onecodevec[detct].uncTime>0.0) {
	  outstream1 << fixed << setprecision(4) << "|" << onecodevec[detct].uncTime;
	} else if(onecodevec[trkstart].uncTime>0.0 && onecodevec[detct].uncTime<=0.0) {
	  outstream1 << "|       ";
	}
	if(onecodevec[trkstart].rmsTime>0.0 && onecodevec[detct].rmsTime>0.0) {
	  outstream1 << fixed << setprecision(4) << "|" << onecodevec[detct].rmsTime << " ";
	} else if(onecodevec[trkstart].rmsTime>0.0 && onecodevec[detct].rmsTime<=0.0) {
	  outstream1 << "|       ";
	}
	outstream1 << "\n";
	detct++;
      } else {
	// This is a new point, but in the same tracklet
	allobsct++;
	if(onecodevec[trkstart].permID.size()>0 && onecodevec[detct].permID.size()>0) {
	  outstream1 << onecodevec[detct].permID << "|";
	} else if(onecodevec[trkstart].permID.size()>0 && onecodevec[detct].permID.size()==0) {
	  for(j=0;j<long(onecodevec[trkstart].permID.size());j++) outstream1 << " ";
	  outstream1 << "|";
	} else if(onecodevec[trkstart].provID.size()>0 && onecodevec[detct].provID.size()>0) {
	  outstream1 << onecodevec[detct].provID << "|";
	} else if(onecodevec[trkstart].provID.size()>0 && onecodevec[detct].provID.size()==0) {
	  for(j=0;j<long(onecodevec[trkstart].provID.size());j++) outstream1 << " ";
	  outstream1 << "|";
	} else if(onecodevec[trkstart].artSat.size()>0 && onecodevec[detct].artSat.size()>0) {
	  outstream1 << onecodevec[detct].artSat << "|";
	} else if(onecodevec[trkstart].artSat.size()>0 && onecodevec[detct].artSat.size()==0) {
	  for(j=0;j<long(onecodevec[trkstart].artSat.size());j++) outstream1 << " ";
	  outstream1 << "|";
	}
	outstream1 << onecodevec[detct].trkSub;
	for(j=onecodevec[detct].trkSub.size();j<7;j++) outstream1 << " ";
	if(trkdigits!=0) {
	  leadzero = 0;
	  if(trkct>0) leadzero = trkdigits - floor(log10(trkct)) - 1;
	  if(leadzero<0) leadzero=0;
	  for(j=0;j<leadzero;j++) outstream1 << "0";
	  outstream1 << trkct;
	}
	outstream1 << "|" << onecodevec[detct].obsSubID;
	if(obsSubID_digits!=0) {
	  outstream1 << ".";
	  leadzero = obsSubID_digits-1;
	  if(allobsct>0) leadzero = obsSubID_digits - floor(log10(allobsct)) - 1;
	  if(leadzero<0) leadzero=0;
	  for(j=0;j<leadzero;j++) outstream1 << "0";
	  outstream1 << allobsct;
	}
	outstream1 << "|" << onecodevec[detct].mode << " ";
	outstream1 << "|" << onecodevec[detct].ObsCode;
	for(j=long(onecodevec[detct].ObsCode.size());j<4;j++) outstream1 << " ";
	// Adjust the input MJD by timeshift, which could be positive (e.g., offset to mid-exposure)
	// or negative (e.g., to convert from TAI to UTC).
	mjd2mpcdate(onecodevec[detct].MJD +(timeshift/SOLARDAY), year, month, day);
	long iday = floor(day);
	double fracday = day - double(iday);
	double hour = floor(24.0*fracday);
	double minute = floor(1440.0*fracday - 60.0*hour);
	double second = fracday*SOLARDAY - 3600.0*hour - 60.0*minute;
	outstream1 <<  "|" << year << "-";
	if(month<10) outstream1 << "0";
	outstream1 << month << "-";
	if(iday<10) outstream1 << "0";
	outstream1 << iday << "T";
	if(hour<10.0) outstream1 << "0";
	outstream1 << fixed << setprecision(0) << hour << ":";
	if(minute<10.0) outstream1 << "0";
	outstream1 << fixed << setprecision(0) << minute << ":";
	if(second<9.9995) outstream1 << "0";
	outstream1 << fixed << setprecision(3) << second << "Z";

	if(onecodevec[detct].RA>=99.99999995) outstream1  << fixed << setprecision(7) << "|" << onecodevec[detct].RA << "";
	else if(onecodevec[detct].RA>=9.99999995) outstream1  << fixed << setprecision(7) << "| " << onecodevec[detct].RA << "";
	else outstream1  << fixed << setprecision(7) << "|  " << onecodevec[detct].RA;
	
	if(onecodevec[detct].Dec<0.0 && fabs(onecodevec[detct].Dec)>=9.99999995) outstream1  << fixed << setprecision(7) << "|" << onecodevec[detct].Dec << "";
	else if(onecodevec[detct].Dec<0.0) outstream1  << fixed << setprecision(7) << "| " << onecodevec[detct].Dec;
	else if(onecodevec[detct].Dec>=9.99999995) outstream1  << fixed << setprecision(7) << "| " << onecodevec[detct].Dec;
	else outstream1  << fixed << setprecision(7) << "|  " << onecodevec[detct].Dec << "";
      
	if((onecodevec[trkstart].rmsRA>0.0 || force_rmsRA) && onecodevec[detct].rmsRA>0.0) {
	  outstream1 << fixed << setprecision(4) << "|" << onecodevec[detct].rmsRA;
	} else if(onecodevec[trkstart].rmsRA>0.0 || force_rmsRA) {
	  outstream1 << "|      ";
	}
	
	if((onecodevec[trkstart].rmsDec>0.0 || force_rmsDec) && onecodevec[detct].rmsDec>0.0) {
	  outstream1 << fixed << setprecision(4)  << "|" << onecodevec[detct].rmsDec;
	} else if(onecodevec[trkstart].rmsDec>0.0 || force_rmsDec) {
	  outstream1 << "|      ";
	}
	
	if(((fabs(onecodevec[trkstart].rmsCorr)<=1.0 && onecodevec[trkstart].rmsRA>0.0 && onecodevec[trkstart].rmsDec>0.0) || force_rmsCorr) && fabs(onecodevec[detct].rmsCorr)<=1.0 && onecodevec[detct].rmsRA>0.0 && onecodevec[detct].rmsDec>0.0) {
	  if(onecodevec[detct].rmsCorr<0.0) outstream1 << fixed << setprecision(4) << "|" << onecodevec[detct].rmsCorr;
	  else outstream1 << fixed << setprecision(4) << "| " << onecodevec[detct].rmsCorr;
	} else if((fabs(onecodevec[trkstart].rmsCorr)<=1.0 && onecodevec[trkstart].rmsRA>0.0 && onecodevec[trkstart].rmsDec>0.0) || force_rmsCorr) {
	  outstream1 << "|       ";
	} 

	if(magprec<=0) {
	  if(onecodevec[detct].mag>=9.5) outstream1  << fixed << setprecision(magprec)  << "|" << " " << onecodevec[detct].mag;
	  else outstream1  << fixed << setprecision(magprec)  << "|" << "  " << onecodevec[detct].mag;
	} else if(magprec==1) {
	  if(onecodevec[detct].mag>=9.95) outstream1  << fixed << setprecision(magprec)  << "|"  << onecodevec[detct].mag;
	  else outstream1  << fixed << setprecision(magprec)  << "|" << " " << onecodevec[detct].mag << "|";
	} else if(magprec==2) {
	  if(onecodevec[detct].mag>=9.995) outstream1  << fixed << setprecision(magprec)  << "|"<< onecodevec[detct].mag;
	  else outstream1  << fixed << setprecision(magprec)  << "|" << " " << onecodevec[detct].mag << "|";
	} else if(magprec==3) {
	  if(onecodevec[detct].mag>=9.9995) outstream1  << fixed << setprecision(magprec)  << "|" << onecodevec[detct].mag;
	  else outstream1  << fixed << setprecision(magprec)  << "|" << " " << onecodevec[detct].mag << "|";
	} else {
	  if(onecodevec[detct].mag>=9.99995) outstream1  << fixed << setprecision(magprec)  << "|" << onecodevec[detct].mag;
	  else outstream1  << fixed << setprecision(magprec)  << "|" << " " << onecodevec[detct].mag;
	}
	
	if((onecodevec[trkstart].rmsMag>0.0 || force_rmsMag) && onecodevec[detct].rmsMag>0.0) {
	  if(magprec<=0) outstream1  << fixed << setprecision(magprec)  << "|" << onecodevec[detct].rmsMag << "     ";
	  else if(magprec==1) outstream1  << fixed << setprecision(magprec) << "|" << onecodevec[detct].rmsMag << "   ";
	  else if(magprec==2) outstream1  << fixed << setprecision(magprec) << "|" << onecodevec[detct].rmsMag << "  ";
	  else if(magprec==3) outstream1  << fixed << setprecision(magprec) << "|" << onecodevec[detct].rmsMag << " ";
	  else outstream1  << fixed << setprecision(magprec) << "|" << onecodevec[detct].rmsMag << "";
	} else if(onecodevec[trkstart].rmsMag>0.0 || force_rmsMag) {
	  outstream1 << "|      ";
	}
	
	outstream1 << "|" << onecodevec[detct].band;
	for(j=onecodevec[detct].band.size();j<4;j++) outstream1 << " ";
	
	if(onecodevec[trkstart].astCat.size()>0 && onecodevec[detct].astCat.size()>0) {
	  outstream1 << "|" << onecodevec[detct].astCat;
	  for(j=onecodevec[detct].astCat.size();j<8;j++) outstream1 << " ";
	} else if(onecodevec[trkstart].astCat.size()>0 && onecodevec[detct].astCat.size()<=0) {
	  outstream1 << "|        ";
	}

	if(onecodevec[trkstart].photCat.size()>0 && onecodevec[detct].photCat.size()>0) {
	  outstream1 << "|" << onecodevec[detct].photCat;
	  for(j=onecodevec[detct].photCat.size();j<8;j++) outstream1 << " ";
	} else if(onecodevec[trkstart].photCat.size()>0 && onecodevec[detct].photCat.size()<=0) {
	  outstream1 << "|        ";
	}
	
	if(onecodevec[trkstart].logSNR>0.0 && onecodevec[detct].logSNR>0.0) {
	  outstream1 << fixed << setprecision(4) << "|" << onecodevec[detct].logSNR << " ";
	} else if(onecodevec[trkstart].logSNR>0.0 && onecodevec[detct].logSNR<=0.0) {
	  outstream1 << "|      ";
	}

	if(onecodevec[trkstart].uncTime>0.0 && onecodevec[detct].uncTime>0.0) {
	  outstream1 << fixed << setprecision(4) << "|" << onecodevec[detct].uncTime;
	} else if(onecodevec[trkstart].uncTime>0.0 && onecodevec[detct].uncTime<=0.0) {
	  outstream1 << "|       ";
	}
	if(onecodevec[trkstart].rmsTime>0.0 && onecodevec[detct].rmsTime>0.0) {
	  outstream1 << fixed << setprecision(4) << "|" << onecodevec[detct].rmsTime << " ";
	} else if(onecodevec[trkstart].rmsTime>0.0 && onecodevec[detct].rmsTime<=0.0) {
	  outstream1 << "|       ";
	}
	outstream1 << "\n";
	detct++;
      }
    }
  }
  outstream1.close();
  return(0);
}
