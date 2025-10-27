// October 24, 2025: label_hldet03.cpp
// Reads two files, an 'unlabled' file and a
// 'labeling' file. Constructs a 4-dimensional k-d tree of the
// 'unlabeled' file, which is expected to be larger than the
// 'labeling' file. Performs a range-query of the
// 'unlabeled' k-d tree around each point in the 'labeling'
// file. If there is a match within a given tolerance, the
// nearest neighbor is than labeled with string ID from the
// 'labeling' file, which supercedes the string ID originally
// read from the 'unlabeled' file, if any. Finally, the newly
// labeled 'unlabeled' file is written out. This output file
// is guaranteed to have the same number of lines as the
// input 'unlabeled' file, and the same string ID's except
// where matches were found.

#define DAY_TO_DEG_CONV 24.0 // Treats one day as equivalent to 24 degrees:
                             // in other words, one second of time is
                             // equivalent to one arcsecond on the sky.
#define IDCOL 1
#define MJDCOL 2
#define RACOL 3
#define DECCOL 4
#define MAGCOL 5
#define BANDCOL 6
#define OBSCODECOL 7
#define COLS_TO_READ 14

#define ID_LABCOL 1
#define MJD_LABCOL 2
#define RA_LABCOL 3
#define DEC_LABCOL 4
#define COLS_TO_READ2 4

#include "solarsyst_dyn_geo01.h"
#include "cmath"
static void show_usage()
{
  cerr << "Usage: label_hldet -unlabeled unlabeled_file -colformat1 column_format_file1 -label labeling_file -colformat2 column_format_file2 -timeoff time_offset_in_seconds -matchrad match_in_arcsec -timescale day_to_deg_conversion -usenearest usenearest -outfile outfile\n";
}

int main(int argc, char *argv[])
{
  ofstream outstream1;
  ifstream instream1;
  string colformatfile1,colformatfile2,stest,idstring;
  vector <hldet> unlabeled_dets = {};
  vector <double> labeling_MJD;
  vector <double> labeling_RA;
  vector <double> labeling_Dec;
  vector <string> labeling_ID;
  double MJD,RA,Dec;
  double time_offset = 0.0;
  int usenearest=0;
  string unlabeled_file, labeling_file, outfile,lnfromfile;
  double timescale = DAY_TO_DEG_CONV;
  double matchrad = 1.0;
  long ulnum,labnum,ulct,labct,i,j,bandlen,status;
  ulnum = labnum = ulct = labct = i = j = bandlen = status = 0;
  int verbose=0;
  long kdroot=0;
  long splitpoint=0;
  long index=0;
  long colreadct,lct;
  long nearest=0;
  double dist = 0.0;
  int reachedeof = 0;
  int startpoint,endpoint;
  point4d_index onepoint = point4d_index(0,0,0,0,0);
  point4d_index querypoint = point4d_index(0,0,0,0,0);
  vector <point4d_index> poolvec;
  KD_point4d_index kdpoint = KD_point4d_index(onepoint,-1,-1,1,0);
  vector <KD_point4d_index> kdvec;
  vector <long> indexvec;
  vector <double> mjdvec;
  double mjdref = 0.0;
  int idcol=IDCOL;
  int mjdcol = MJDCOL;
  int racol = RACOL;
  int deccol = DECCOL;
  int magcol = MAGCOL;
  int bandcol = BANDCOL;
  int obscodecol = OBSCODECOL;
  int trail_len_col, trail_PA_col, sigmag_col, sig_across_col, sig_along_col, known_obj_col, det_qual_col;
  trail_len_col = trail_PA_col = sigmag_col = sig_across_col = -1;
  sig_along_col = known_obj_col = det_qual_col = -1;
  int id_labcol = ID_LABCOL;
  int mjd_labcol = MJD_LABCOL;
  int ra_labcol = RA_LABCOL;
  int dec_labcol = DEC_LABCOL;
  int colformatfile1_set=0;
  int colformatfile2_set=0;

  if(argc<7) {
    show_usage();
    return(1);
  }
  
  i=1;
  while(i<argc) {
    cout << "Checking out argv[" << i << "] = " << argv[i] << ".\n";
    if(string(argv[i]) == "-unlabeled" || string(argv[i]) == "-unlab" || string(argv[i]) == "-ul" || string(argv[i]) == "-nolab" || string(argv[i]) == "--unlabeled" || string(argv[i]) == "--unlab") {
      if(i+1 < argc) {
	//There is still something to read;
	unlabeled_file=argv[++i];
	i++;
      }
      else {
	cerr << "Unlabeled file keyword supplied with no corresponding argument\n";
	show_usage();
	return(1);
      }
    } else if(string(argv[i]) == "-colformat1" || string(argv[i]) == "-format1"  || string(argv[i]) == "-col1" || string(argv[i]) == "-cf1" || string(argv[i]) == "-colfmt1" || string(argv[i]) == "--colformat1" || string(argv[i]) == "--columnformat1" || string(argv[i]) == "--cformat1") {
      if(i+1 < argc) {
	//There is still something to read;
	colformatfile1=argv[++i];
	colformatfile1_set = 1;
	i++;
      }
      else {
	cerr << "Column format file 1 keyword supplied with no corresponding argument\n";
	show_usage();
	return(1);
      }
    } else if(string(argv[i]) == "-labeled" || string(argv[i]) == "-labeling" || string(argv[i]) == "-lab" || string(argv[i]) == "-labeling_file" || string(argv[i]) == "-label" || string(argv[i]) == "--labeled") {
      if(i+1 < argc) {
	//There is still something to read;
	labeling_file=argv[++i];
	i++;
      }
      else {
	cerr << "Labeling file keyword supplied with no corresponding argument\n";
	show_usage();
	return(1);
      }
    } else if(string(argv[i]) == "-colformat2" || string(argv[i]) == "-format2"  || string(argv[i]) == "-col2" || string(argv[i]) == "-cf2" || string(argv[i]) == "-colfmt2" || string(argv[i]) == "--colformat2" || string(argv[i]) == "--columnformat2" || string(argv[i]) == "--cformat2") {
      if(i+1 < argc) {
	//There is still something to read;
	colformatfile2=argv[++i];
	colformatfile2_set = 1;
	i++;
      }
      else {
	cerr << "Column format file 1 keyword supplied with no corresponding argument\n";
	show_usage();
	return(1);
      }
    } else if(string(argv[i]) == "-matchrad") {
      if(i+1 < argc) {
	//There is still something to read;
	matchrad=stod(argv[++i]);
	i++;
      }
      else {
	cerr << "Matching radius keyword supplied with no corresponding argument\n";
	show_usage();
	return(1);
      }
    } else if(string(argv[i]) == "-timeoff" || string(argv[i]) == "-time_offset") {
      if(i+1 < argc) {
	//There is still something to read;
	time_offset=stod(argv[++i]);
	i++;
      }
      else {
	cerr << "Time offset keyword supplied with no corresponding argument\n";
	show_usage();
	return(1);
      }
    } else if(string(argv[i]) == "-timescale") {
      if(i+1 < argc) {
	//There is still something to read;
	timescale=stod(argv[++i]);
	i++;
      }
      else {
	cerr << "Timescale keyword supplied with no corresponding argument\n";
	show_usage();
	return(1);
      }
    } else if(string(argv[i]) == "-usenearest" || string(argv[i]) == "-use_nearest" || string(argv[i]) == "-nearest") {
      if(i+1 < argc) {
	//There is still something to read;
	usenearest=stoi(argv[++i]);
	i++;
      }
      else {
	cerr << "use_nearest keyword supplied with no corresponding argument\n";
	show_usage();
	return(1);
      }
    } else if(string(argv[i]) == "-outfile" || string(argv[i]) == "-out" || string(argv[i]) == "--outfile") {
      if(i+1 < argc) {
	//There is still something to read;
	outfile=argv[++i];
	i++;
      }
      else {
	cerr << "Labeling file keyword supplied with no corresponding argument\n";
	show_usage();
	return(1);
      }
    } else {
      cerr << "Warning: unrecognized keyword " << argv[i] <<"\n";
      i++;
    }
  }
  cout.precision(17);  

  cout << "Input unlabeled file is called " << unlabeled_file << "\n";
  if(colformatfile1_set==0) {
    cout << "WARNING: No column formatting file supplied for the unlabeled file\n";
  } else cout << "Column formatting file for unlabeled input: " << colformatfile1 << "\n";
  cout << "Input labeling file is called " << labeling_file << "\n";
  if(colformatfile2_set==0) {
    cout << "WARNING: No column formatting file supplied for the labeling file\n";
  } else cout << "Column formatting file for labeling input: " << colformatfile2 << "\n";
  cout << "Matching radius will be " << matchrad << " arcseconds\n";
  cout << "Timescale for converting 1 day of time to equivalent degrees will be " << timescale << "\n";
  cout << "Time offset to be applied to MJD values in the labeling file is " << time_offset << " seconds\n";
  cout << "Output file will be called " << outfile << "\n";
  
  // Read the column formatting file for the unlabeled data
  if(colformatfile1.size()>0)
    {
      instream1.open(colformatfile1);
      if(!instream1)  {
	cerr << "ERROR: unable to open input file " << colformatfile1 << "\n";
	return(1);
      }
      colreadct=0;
      while(!instream1.eof() && !instream1.fail() && !instream1.bad() && colreadct<COLS_TO_READ) {
	instream1 >> stest;
	if(stest == "MJDCOL") {
	  instream1 >> mjdcol;
	  if(!instream1.eof() && !instream1.fail() && !instream1.bad()) colreadct++;
	} else if(stest == "RACOL") {
	  instream1 >> racol;
	  if(!instream1.eof() && !instream1.fail() && !instream1.bad()) colreadct++;
	} else if(stest == "DECCOL") {
	  instream1 >> deccol;
	  if(!instream1.eof() && !instream1.fail() && !instream1.bad()) colreadct++;
	} else if(stest == "MAGCOL") {
	  instream1 >> magcol;
	  if(!instream1.eof() && !instream1.fail() && !instream1.bad()) colreadct++;
	} else if(stest == "TRAILLENCOL") {
	  instream1 >> trail_len_col;
	  if(!instream1.eof() && !instream1.fail() && !instream1.bad()) colreadct++;
	} else if(stest == "TRAILPACOL") {
	  instream1 >> trail_PA_col;
	  if(!instream1.eof() && !instream1.fail() && !instream1.bad()) colreadct++;
	} else if(stest == "SIGMAGCOL") {
	  instream1 >> sigmag_col;
	  if(!instream1.eof() && !instream1.fail() && !instream1.bad()) colreadct++;
	} else if(stest == "SIGACROSSCOL") {
	  instream1 >> sig_across_col;
	  if(!instream1.eof() && !instream1.fail() && !instream1.bad()) colreadct++;
	} else if(stest == "SIGALONGCOL") {
	  instream1 >> sig_along_col;
	  if(!instream1.eof() && !instream1.fail() && !instream1.bad()) colreadct++;
	} else if (stest == "IDCOL") {
	  instream1 >> idcol;
	  if(!instream1.eof() && !instream1.fail() && !instream1.bad()) colreadct++;
	} else if(stest == "BANDCOL") {
	  instream1 >> bandcol;
	  if(!instream1.eof() && !instream1.fail() && !instream1.bad()) colreadct++;
	} else if(stest == "OBSCODECOL") {
	  instream1 >> obscodecol;
	  if(!instream1.eof() && !instream1.fail() && !instream1.bad()) colreadct++;
	} else if(stest == "KNOWNOBJCOL") {
	  instream1 >> known_obj_col;
	  if(!instream1.eof() && !instream1.fail() && !instream1.bad()) colreadct++;
	} else if(stest == "DETQUALCOL") {
	  instream1 >> det_qual_col;
	  if(!instream1.eof() && !instream1.fail() && !instream1.bad()) colreadct++;
	} else {
	  cout << "WARNING: unrecognized string " << stest << " read from column formatting file\n";
	}
      }
      instream1.close();
      if(colreadct<COLS_TO_READ) {
	cout << "WARNING: only " << colreadct << " column specifications, of " << COLS_TO_READ << " expected, were read from column format file " << colformatfile1 << ".\n";
      }
    }

  // Read the column formatting file for the labeling data
  if(colformatfile2.size()>0)
    {
      instream1.open(colformatfile2);
      if(!instream1)  {
	cerr << "ERROR: unable to open input file " << colformatfile2 << "\n";
	return(1);
      }
      colreadct=0;
      while(!instream1.eof() && !instream1.fail() && !instream1.bad() && colreadct<COLS_TO_READ2) {
	instream1 >> stest;
	if(stest == "MJDCOL") {
	  instream1 >> mjd_labcol;
	  if(!instream1.eof() && !instream1.fail() && !instream1.bad()) colreadct++;
	} else if(stest == "RACOL") {
	  instream1 >> ra_labcol;
	  if(!instream1.eof() && !instream1.fail() && !instream1.bad()) colreadct++;
	} else if(stest == "DECCOL") {
	  instream1 >> dec_labcol;
	  if(!instream1.eof() && !instream1.fail() && !instream1.bad()) colreadct++;
	} else if (stest == "IDCOL") {
	  instream1 >> id_labcol;
	  if(!instream1.eof() && !instream1.fail() && !instream1.bad()) colreadct++;
	} else {
	  cout << "WARNING: unrecognized string " << stest << " read from column formatting file\n";
	}
      }
      instream1.close();
      if(colreadct<COLS_TO_READ2) {
	cout << "ERROR: only " << colreadct << " column specifications, of " << COLS_TO_READ2 << " expected, were read from column format file " << colformatfile2 << ".\n";
	return(1);
      }
    }
  
  // Read input unlabeled file
  unlabeled_dets={};
  status = read_detection_filemt2(unlabeled_file, mjdcol, racol, deccol, magcol, idcol, bandcol, obscodecol, trail_len_col, trail_PA_col, sigmag_col, sig_across_col, sig_along_col, known_obj_col, det_qual_col, unlabeled_dets, verbose, 1);  
  if(status==0) { 
    cout << "Input file " << unlabeled_file << " read successfully to the end.\n";
  }
  else if(status==1) {
    cerr << "Warning: file read failed\n";
    return(1);
  } else if(status==2) {
    cerr << "Warning: file possibly corrupted\n";
    return(2);
  } else {
    cerr << "Warning: unknown file read problem\n";
    return(3);
  }
  ulnum = unlabeled_dets.size();
  cout << "Read " << ulnum << " data lines from observation file " << unlabeled_file << "\n";
  
  // Read input labeling file
  labeling_MJD = labeling_RA = labeling_Dec = {};
  labeling_ID = {};
  lct=0;
  instream1.open(labeling_file);
  if(!instream1) {
    cerr << "can't open input file " << labeling_file << "\n";
    return(1);
  }
  // Skip one-line header
  getline(instream1,lnfromfile);
  lct++;
  //cout << "Line " << lct << ": " << lnfromfile << "\n";
  reachedeof = 0;
  while(reachedeof==0) {
    getline(instream1,lnfromfile);
    lct++;
    //cout << "Line " << lct << ": " << lnfromfile << "\n";
    if(!instream1.eof() && !instream1.fail() && !instream1.bad()) ; // Read on.
    else if(instream1.eof()) reachedeof=1; //End of file, fine.
    else if(instream1.fail()) reachedeof=-1; //Something wrong, warn
    else if(instream1.bad()) reachedeof=-2; //Worse problem, warn
    MJD = 0.0l;
    RA = Dec = -999.9;
    idstring = "";
    startpoint = endpoint = j = 0;
    while(endpoint<long(lnfromfile.size()) && startpoint<long(lnfromfile.size()) && lnfromfile.size()>=15 && reachedeof == 0) {
      endpoint = get_csv_string01(lnfromfile,stest,startpoint);
      j++;
      //cout << "j, endpoint, startpoint, stest, cols: " << j << " " << endpoint << " " << startpoint << " " << stest << " " << id_labcol << " " << mjd_labcol << " " << ra_labcol << " " << dec_labcol << "\n";
      if(j==id_labcol) {
	if(endpoint>0) {
	  idstring = stest;
	}
      } else if(j==mjd_labcol) {
	if(endpoint>0) {
	  try { MJD = stod(stest); }
	  catch(...) { cerr << "ERROR: cannot read MJD string " << stest << " from line " << lnfromfile << "\n"; }
	}
      } else if(j==ra_labcol) {
	if(endpoint>0) {
	  try { RA = stod(stest); }
	  catch(...) { cerr << "ERROR: cannot read RA string " << stest << " from line " << lnfromfile << "\n"; }
	}
      } else if(j==dec_labcol) {
	if(endpoint>0) {
	  try { Dec = stod(stest); }
	  catch(...) { cerr << "ERROR: cannot read DEC string " << stest << " from line " << lnfromfile << "\n"; }
	}
      }
      startpoint = endpoint+1;
    }
    if(reachedeof == 0 && lnfromfile.size()>=15) {
      if(MJD==0.0) {
	cerr << "ERROR: MJD not read from line " << labeling_MJD.size()+1 << " of input labeling file " << labeling_file << "!\n";
	cerr << "Here is the line: " << lnfromfile << "\n";
	return(2);
      }
      if(RA==-999.9) {
	cerr << "ERROR: RA not read from line " << labeling_RA.size()+1 << " of input labeling file " << labeling_file << "!\n";
	cerr << "Here is the line: " << lnfromfile << "\n";
	return(2);
      }
      if(Dec==-999.9) {
	cerr << "ERROR: Dec not read from line " << labeling_Dec.size()+1 << " of input labeling file " << labeling_file << "!\n";
	cerr << "Here is the line: " << lnfromfile << "\n";
	return(2);
      }
      if(idstring.size()<=0) {
	cerr << "ERROR: Object ID not read from line " << labeling_ID.size()+1 << " of input labeling file " << labeling_file << "!\n";
	cerr << "Here is the line: " << lnfromfile << "\n";
	return(2);
      }
      labeling_MJD.push_back(MJD + time_offset/SOLARDAY);
      labeling_RA.push_back(RA);
      labeling_Dec.push_back(Dec);
      labeling_ID.push_back(idstring);
      //cout << "MJD, lengths: " << MJD << " " << RA << " " << Dec << " " << idstring << " " << labeling_MJD.size() << " " << labeling_RA.size() << " " << labeling_Dec.size() << " " << labeling_ID.size() << "\n";
    }
  }
  instream1.close();
  
  if(reachedeof==1) { 
    cout << "Input file " << labeling_file << " read successfully to the end.\n";
  } else if(reachedeof==0) {
    cerr << "ERROR: Stopped reading file " << labeling_file << " before the end\n";
    return(1);
  } else if(reachedeof==-1) {
    cerr << "Warning: file read failed\n";
    return(1);
  } else if(reachedeof==-2) {
    cerr << "Warning: file possibly corrupted\n";
    return(2);
  }

  labnum = labeling_MJD.size();
  if(labnum != long(labeling_RA.size()) || labnum != long(labeling_Dec.size()) || labnum != long(labeling_ID.size())) {
    cerr << "ERROR: length mismatch in vectors read from labeling file: " << labnum << " " << labeling_RA.size() << " " << labeling_Dec.size() << " " << labeling_ID.size() << "\n";
  }
  cout << "Read " << labnum << " data lines from observation file " << labeling_file << "\n";

  // Find the median MJD for unlabeled detections
  mjdvec = {};
  for(ulct=0;ulct<ulnum;ulct++) mjdvec.push_back(unlabeled_dets[ulct].MJD);
  mjdref = dmedian(mjdvec);
  cout << "Median MJD of unlabled detections is " << mjdref << "\n";
  
  // Load pool of unlabeled sources
  poolvec = {};
  for(ulct=0;ulct<ulnum;ulct++) {
    // Note that timescale is supposed to indicate the number of
    // degrees on the sky that is the equivalent of one calendar day.
    // Then, since the trigonometric projection of sky-angles converts
    // degrees into radians, it is also necessary to divide the
    // time-quantity by DEGPRAD.
    onepoint = point4d_index((unlabeled_dets[ulct].MJD-mjdref)*timescale/DEGPRAD,cos(unlabeled_dets[ulct].RA/DEGPRAD)*cos(unlabeled_dets[ulct].Dec/DEGPRAD),sin(unlabeled_dets[ulct].RA/DEGPRAD)*cos(unlabeled_dets[ulct].Dec/DEGPRAD), sin(unlabeled_dets[ulct].Dec/DEGPRAD), ulct);
    poolvec.push_back(onepoint);
  }
  cout << "Finished loading pool of unlabeled sources\n";
  // Form KDtree
  kdvec={};
  kdroot = splitpoint = 0;
  splitpoint=medind_4d_index(poolvec,1);
  kdpoint = KD_point4d_index(poolvec[splitpoint],-1,-1,1,0);
  kdvec.push_back(kdpoint);
  kdtree_4d_index(poolvec,1,splitpoint,kdroot,kdvec);
  cout << "Finished constructing k-d tree sources\n";

  // Loop over pool of labeled sources, performing a nearest-neighbor
  // query of the unlabled k-d tree for each of them.
  for(labct=0;labct<labnum;labct++) {
    querypoint = point4d_index((labeling_MJD[labct]-mjdref)*timescale/DEGPRAD,cos(labeling_RA[labct]/DEGPRAD)*cos(labeling_Dec[labct]/DEGPRAD),sin(labeling_RA[labct]/DEGPRAD)*cos(labeling_Dec[labct]/DEGPRAD), sin(labeling_Dec[labct]/DEGPRAD), labct);
    if(usenearest!=1) {
      indexvec={};
      status = kdrange_4d_index(kdvec, querypoint, matchrad/ASECPRAD, indexvec);
      for(i=0;i<long(indexvec.size());i++) {
	// Label the match with an index that maps all the way back to unlabeled_dets
	index = kdvec[indexvec[i]].point.index;
	// Change the label for unlabeled_dets
	stringncopy01(unlabeled_dets[index].idstring,labeling_ID[labct],SHORTSTRINGLEN);
	unlabeled_dets[index].known_obj=999;
      }
    } else {
      nearest = kdnearest_4d_index(kdvec, querypoint);
      dist = ASECPRAD*sqrt(point4d_index_dist2(querypoint, kdvec[nearest].point));
      if(dist<=matchrad) {
	// Label the match with an index that maps all the way back to unlabeled_dets
	index = kdvec[nearest].point.index;
	// Change the label for unlabeled_dets
	stringncopy01(unlabeled_dets[index].idstring,labeling_ID[labct],SHORTSTRINGLEN);
	unlabeled_dets[index].known_obj=999;
      }
    }
    if(labct%1000==0) {
      cout << "Running query for labeling point " << labct << "\n";
      //cout << querypoint.t << " " << querypoint.x << " " << querypoint.y << " " << querypoint.z << ", " << indexvec.size() << " matches found\n";
    }
  }
  
  outstream1.open(outfile);
  outstream1 << "#MJD,RA,Dec,mag,trail_len,trail_PA,sigmag,sig_across,sig_along,image,idstring,band,obscode,known_obj,det_qual,origindex\n";
  for(ulct=0;ulct<ulnum;ulct++) {
    outstream1 << fixed << setprecision(7) << unlabeled_dets[ulct].MJD << "," << unlabeled_dets[ulct].RA << "," << unlabeled_dets[ulct].Dec << ",";
    outstream1 << fixed << setprecision(4) << unlabeled_dets[ulct].mag << ",";
    outstream1 << fixed << setprecision(2) << unlabeled_dets[ulct].trail_len << "," << unlabeled_dets[ulct].trail_PA << ",";
    outstream1 << fixed << setprecision(4) << unlabeled_dets[ulct].sigmag << ",";
    outstream1 << fixed << setprecision(3) << unlabeled_dets[ulct].sig_across << "," << unlabeled_dets[ulct].sig_along << ",";
    outstream1 << unlabeled_dets[ulct].image << "," << unlabeled_dets[ulct].idstring << "," << unlabeled_dets[ulct].band << ",";
    outstream1 << unlabeled_dets[ulct].obscode << "," << unlabeled_dets[ulct].known_obj << ","; 
    outstream1 << unlabeled_dets[ulct].det_qual << "," << unlabeled_dets[ulct].index << "\n"; 
  }
  outstream1.close();

  return(0);
}
