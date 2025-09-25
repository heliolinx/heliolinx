// September 23, 2025: merge_tracklet_files.cpp:
// Given a file listing outputs from multiple different runs of make_tracklets,
// merge all of the runs into one consistent set of make_tracklets output files.
// Note that each set of make_tracklets file contains the following, in order:
// 1. Image file with MJD, RA, Dec, and observer coordinates at the time of each image
// 2. Paired detection file with indices to the image file
// 3. tracklet summary file, with indices to the image file
// 4. trk2det file, with indices to the paired detection file
// Hence, when merging files, the indices in the paired detection file
// and the tracklet summary file have to be changed to reflect the
// new, expanded set of images, and the indices in the trk2det file
// have to be changed to reflect the new, expanded set of paired
// detections.

#include "solarsyst_dyn_geo01.h"
#include "cmath"

#define MAXCOMB 100
#define ARCSEC_MATCHTOL 0.1
#define NEG_MJD_MARK -999999999.9

static void show_usage()
{
  cerr << "Usage: merge_tracklet_files -inlist input_file_list -outputs output_image file output_paired_detection_file output_tracklet_file output_tracklet-to-detection_file -matchtol detection_match_tol_arcsec\n";
}
    
int main(int argc, char *argv[])
{
  string inlist;
  string in_image_files[MAXCOMB], in_pairdet_files[MAXCOMB];
  string in_tracklet_files[MAXCOMB], in_trk2det_files[MAXCOMB];
  string imagefile,pairdetfile,trackletfile,trk2detfile;
  string imagefile_out,pairdetfile_out,trackletfile_out,trk2detfile_out;
  vector <hlimage> inimage;
  vector <hldet> inpairdets;
  vector <tracklet> intracklets;
  vector <longpair> intrk2det;
  vector <long_index> tracklet_hash;
  vector <hlimage> masterimage;
  vector <hlimage> tempimage;
  vector <hldet> masterpairdets;
  vector <hldet> temppairdets;
  vector <hldet> oldpairdets;
  vector <tracklet> mastertracklets;
  vector <longpair> master_trk2det;
  int verbose=0;
  ofstream outstream1;
  int status=0;
  ifstream instream1;
  int filenum,filect,tmpim;
  filenum = filect = tmpim = 0;
  vector <long> temp_imindex;
  vector <long> in_imindex;
  vector <long> in_detindex;
  vector <long> old_detindex;
  long masterimnum,inimnum,tempimnum,masterimct,inimct,tempimct;
  masterimnum = inimnum = tempimnum = masterimct = inimct = tempimct = 0;
  int matchfound=0;
  long masterdetnum,indetnum,olddetnum,detct,i,j;
  masterdetnum = indetnum = olddetnum = detct = i = j = 0;
  long master_trknum,in_trknum,trkct;
  long_index one_lindex = long_index(0,0);
  vector <long_index> partnervec;
  double arcsec_matchtol = ARCSEC_MATCHTOL;
  double deg_matchtol = arcsec_matchtol/3600.0;
  vector <long> trkvec;
  longpair onepair = longpair(0,0);
  
  if(argc!=8 && argc!=10) {
    cout << "Need argc=8 or 10, got " << argc << "\n";
    show_usage();
    return(1);
  }
  
  i=1;
  while(i<argc) {
    cout << "Checking out argv[" << i << "] = " << argv[i] << ".\n";
    if(string(argv[i]) == "-inlist") {
      if(i+3 < argc) {
	//There is still something to read;
        inlist=argv[++i];
	i++;
      }
      else {
	cerr << "Input file list keyword supplied with no corresponding argument\n";
	show_usage();
	return(1);
      }
    } else if(string(argv[i]) == "-outputs" || string(argv[i]) == "-output" || string(argv[i]) == "-outfile" || string(argv[i]) == "-outfiles" || string(argv[i]) == "-out") {
      if(i+3 < argc) {
	//There is still something to read;
        imagefile_out=argv[++i];
        pairdetfile_out=argv[++i];
        trackletfile_out=argv[++i];
        trk2detfile_out=argv[++i];
	i++;
      }
      else {
	cerr << "Output file keyword supplied with fewer than the required four arguments\n";
	show_usage();
	return(1);
      }
    } else if(string(argv[i]) == "-matchdist" || string(argv[i]) == "-matchtol" || string(argv[i]) == "-matchrad") {
      if(i+1 < argc) {
	//There is still something to read;
	arcsec_matchtol=stod(argv[++i]);
	deg_matchtol = arcsec_matchtol/3600.0;
	i++;
      }
      else {
	cerr << "Match tolerance keyword supplied with no corresponding argument\n";
	show_usage();
	return(1);
      }
    } else {
      cerr << "Warning: unrecognized keyword or argument " << argv[i] << "\n";
      i++;
    }
  }

  instream1.open(inlist);
  if(!instream1) {
    cerr << "can't open input file " << inlist << "\n";
    return(1);
  }
  filect=0;
  while(!instream1.eof() && !instream1.fail() && !instream1.bad()) {
    imagefile = pairdetfile = trackletfile = trk2detfile = "";
    instream1 >> imagefile;
    instream1 >> pairdetfile;
    instream1 >> trackletfile;
    instream1 >> trk2detfile;
    if(imagefile.size()>0 && pairdetfile.size()>0 && trackletfile.size()>0 && trk2detfile.size()>0) {
      in_image_files[filect] = imagefile;
      in_pairdet_files[filect] = pairdetfile;
      in_tracklet_files[filect] = trackletfile;
      in_trk2det_files[filect] = trk2detfile;
      filect++;
    }
  }
  filenum = filect;
  cout << "Read " << filenum << " lines from input file list " << inlist << "\n";
  instream1.close();
    
  // Read first set of input files into the master vectors
  masterimage={};
  status=read_image_file2(in_image_files[0], masterimage);
  if(status!=0) {
    cerr << "ERROR: could not successfully read image file " << in_image_files[0] << "\n";
    cerr << "read_image_file2 returned status = " << status << ".\n";
   return(1);
  }
  cout << "Read " << masterimage.size() << " data lines from image file " << in_image_files[0] << "\n";
  masterpairdets={};
  status=read_pairdet_file(in_pairdet_files[0], masterpairdets , verbose);
  if(status!=0) {
    cerr << "ERROR: could not successfully read paired detection file " << in_pairdet_files[0] << "\n";
    cerr << "read_pairdet_file returned status = " << status << ".\n";
   return(1);
  }
  cout << "Read " << masterpairdets.size() << " data lines from paired detection file " << in_pairdet_files[0] << "\n";
  mastertracklets={};
  status=read_tracklet_file(in_tracklet_files[0], mastertracklets, verbose);
  if(status!=0) {
    cerr << "ERROR: could not successfully read tracklet file " << in_tracklet_files[0] << "\n";
    cerr << "read_tracklet_file returned status = " << status << ".\n";
   return(1);
  }
  cout << "Read " << mastertracklets.size() << " data lines from tracklet file " << in_tracklet_files[0] << "\n";
  master_trk2det={};
  status=read_longpair_file(in_trk2det_files[0], master_trk2det, verbose);
  if(status!=0) {
    cerr << "ERROR: could not successfully read trk2det file " << in_trk2det_files[0] << "\n";
    cerr << "read_longpair_file returned status = " << status << ".\n";
   return(1);
  }
  cout << "Read " << master_trk2det.size() << " data lines from trk2det file " << in_trk2det_files[0] << "\n";

  // Loop over remaining file sets
  for(filect=1;filect<filenum;filect++) {
    inimage={};
    status=read_image_file2(in_image_files[filect], inimage);
    if(status!=0) {
      cerr << "ERROR: could not successfully read image file " << in_image_files[filect] << "\n";
      cerr << "read_image_file2 returned status = " << status << ".\n";
      return(1);
    }
    cout << "Read " << inimage.size() << " data lines from image file " << in_image_files[filect] << "\n";
    inpairdets={};
    status=read_pairdet_file(in_pairdet_files[filect], inpairdets, verbose);
    if(status!=0) {
      cerr << "ERROR: could not successfully read paired detection file " << in_pairdet_files[filect] << "\n";
      cerr << "read_pairdet_file returned status = " << status << ".\n";
      return(1);
    }
    cout << "Read " << inpairdets.size() << " data lines from paired detection file " << in_pairdet_files[filect] << "\n";
    intracklets={};
    status=read_tracklet_file(in_tracklet_files[filect], intracklets, verbose);
    if(status!=0) {
      cerr << "ERROR: could not successfully read tracklet file " << in_tracklet_files[filect] << "\n";
      cerr << "read_tracklet_file returned status = " << status << ".\n";
      return(1);
    }
    cout << "Read " << intracklets.size() << " data lines from tracklet file " << in_tracklet_files[filect] << "\n";
    intrk2det={};
    status=read_longpair_file(in_trk2det_files[filect], intrk2det, verbose);
    if(status!=0) {
      cerr << "ERROR: could not successfully read trk2det file " << in_trk2det_files[filect] << "\n";
      cerr << "read_longpair_file returned status = " << status << ".\n";
      return(1);
    }
    cout << "Read " << intrk2det.size() << " data lines from trk2det file " << in_trk2det_files[filect] << "\n";
    // Copy masterimage into tempimage
    tempimage = masterimage;
    inimnum = inimage.size();
    tempimnum = tempimage.size();
    make_lvec(masterimnum, temp_imindex);
    make_lvec(inimnum, in_imindex);
    // Concatenate the new image catalog
    for(i=0;i<inimnum;i++) masterimage.push_back(inimage[i]);
    // Sort the concatenated catalog
    sort(masterimage.begin(), masterimage.end(), early_hlimage());
    // De-duplicate the masterimage vector
    cout << "Master image set has " << masterimage.size() << " images prior to de-duplication\n";
    i=1;
    while(i<long(masterimage.size()) && masterimage.size()>0) {
      j=i-1;
      int deleted_one = 0;
      while(masterimage[j].MJD >= masterimage[i].MJD-IMAGETIMETOL/SOLARDAY && j>=0) {
	if(fabs(masterimage[i].MJD - masterimage[j].MJD)<IMAGETIMETOL/SOLARDAY && stringnmatch01(masterimage[i].obscode,masterimage[j].obscode,3)==0) {
	  // Image i is a duplicate of image j. Erase image i from the vector
	  masterimage.erase(masterimage.begin()+i);
	  deleted_one=1;
	  break;
	} else j--;
      }
      if(!deleted_one) i++;
    }
    masterimnum = masterimage.size();
    cout << "Master image set has " << masterimnum << " images after de-duplication\n";

    // Load index vector connecting the old masterimage vector with the new
    masterimct=0;
    tempimct=0;
    make_lvec(tempimnum, temp_imindex);
    for(tempimct=0;tempimct<tempimnum;tempimct++) {
      // Find a match for tempimage[tempimct] in the master image vector
      // Go backwards through masterimage until we are certain to be before tempimage[tempimct]
      while(masterimct>0 && masterimage[masterimct].MJD >= tempimage[tempimct].MJD - IMAGETIMETOL/SOLARDAY) masterimct--;
      // Now go forward until we find a match
      matchfound=0;
      while(masterimct<masterimnum && masterimage[masterimct].MJD <= tempimage[tempimct].MJD + IMAGETIMETOL/SOLARDAY) {
	if(fabs(tempimage[tempimct].MJD - masterimage[masterimct].MJD)<IMAGETIMETOL/SOLARDAY && stringnmatch01(tempimage[tempimct].obscode, masterimage[masterimct].obscode,3)==0) {
	  matchfound=1;
	  temp_imindex[tempimct] = masterimct;
	  break;
	}
	masterimct++;
      }
      if(!matchfound) {
	cerr << "ERROR: no match found to image " << tempimct << " " << tempimage[tempimct].MJD << " " << tempimage[tempimct].RA << " " << tempimage[tempimct].Dec << " " << tempimage[tempimct].obscode << "\n";
	return(3);
      }
    }
    cout << "Created index-matching vector for tempimage to masterimage\n";
    masterimct=0;
    inimct=0;
    make_lvec(inimnum, in_imindex);
    for(inimct=0;inimct<inimnum;inimct++) {
      // Find a match for inimage[inimct] in the master image vector
      // Go backwards through masterimage until we are certain to be before inimage[inimct]
      while(masterimct>0 && masterimage[masterimct].MJD >= inimage[inimct].MJD - IMAGETIMETOL/SOLARDAY) masterimct--;
      // Now go forward until we find a match
      matchfound=0;
      while(masterimct<masterimnum && masterimage[masterimct].MJD <= inimage[inimct].MJD + IMAGETIMETOL/SOLARDAY) {
	if(fabs(inimage[inimct].MJD - masterimage[masterimct].MJD)<IMAGETIMETOL/SOLARDAY && stringnmatch01(inimage[inimct].obscode, masterimage[masterimct].obscode,3)==0) {
	  matchfound=1;
	  in_imindex[inimct] = masterimct;
	  break;
	}
	masterimct++;
      }
      if(!matchfound) {
	cerr << "ERROR: no match found to image " << inimct << " " << inimage[inimct].MJD << " " << inimage[inimct].RA << " " << inimage[inimct].Dec << " " << inimage[inimct].obscode << "\n";
	return(3);
      }
    }
    cout << "Created index-matching vector for inimage to masterimage\n";
    // Update the image indices in the old pairdets file
    masterdetnum = masterpairdets.size();
    for(detct=0;detct<masterdetnum;detct++) {
      tmpim = masterpairdets[detct].image;
      masterpairdets[detct].image = temp_imindex[tmpim];
    }
    // Update the image indices in the latest input pairdets file
    indetnum = inpairdets.size();
    for(detct=0;detct<indetnum;detct++) {
      tmpim = inpairdets[detct].image;
      inpairdets[detct].image = in_imindex[tmpim];
    }
    // Update the image indices in the tracklet files
    master_trknum = mastertracklets.size();
    for(trkct=0;trkct<master_trknum;trkct++) {
      tmpim = mastertracklets[trkct].Img1;
      mastertracklets[trkct].Img1 = temp_imindex[tmpim];
      tmpim = mastertracklets[trkct].Img2;
      mastertracklets[trkct].Img2 = temp_imindex[tmpim];
    }
    in_trknum = intracklets.size();
    for(trkct=0;trkct<in_trknum;trkct++) {
      tmpim = intracklets[trkct].Img1;
      intracklets[trkct].Img1 = in_imindex[tmpim];
      tmpim = intracklets[trkct].Img2;
      intracklets[trkct].Img2 = in_imindex[tmpim];
    }

    // Concatenate the pairdets files
    oldpairdets = temppairdets = masterpairdets;
    olddetnum = oldpairdets.size();
    indetnum = inpairdets.size();
    for(detct=0;detct<indetnum;detct++) temppairdets.push_back(inpairdets[detct]);
    // Sort temppairdets by image (then Dec, then RA)
    sort(temppairdets.begin(), temppairdets.end(), imgsort_hldet());
    partnervec = {};
    i=0;
    for(detct=0;detct<long(temppairdets.size());detct++) {
      // Decrement i until we are certain to be below any duplicates
      while(i>0 && (temppairdets[i].image > temppairdets[detct].image || (temppairdets[i].image == temppairdets[detct].image && temppairdets[i].Dec >= temppairdets[detct].Dec - deg_matchtol))) i--;
      long matchnum=0;
      // Find all possible duplicates for temppairdets[detct]
      while(i<long(temppairdets.size()) && (temppairdets[i].image < temppairdets[detct].image || (temppairdets[i].image == temppairdets[detct].image && temppairdets[i].Dec <= temppairdets[detct].Dec + deg_matchtol))) {
	if(i==detct) matchnum++;
	else if(temppairdets[i].image == temppairdets[detct].image) {
	  // Calculate the distance from detection i to detection detct
	  double dist =  distradec01(temppairdets[detct].RA, temppairdets[detct].Dec, temppairdets[i].RA, temppairdets[i].Dec);
	  if(dist<=deg_matchtol) matchnum++;
	}
	i++;
      }
      one_lindex = long_index(matchnum,detct);
      partnervec.push_back(one_lindex);
    }
    if(long(partnervec.size()) != long(temppairdets.size())) {
      cerr << "ERROR: partnervec improperly loaded, size mismatch " << partnervec.size() << " vs. " << temppairdets.size() << "\n";
      return(2);
    }
    // Sort partnervec by number of matches
    sort(partnervec.begin(), partnervec.end(), lower_long_index());
    // Loop over temppairdets in order of decreasing numbers of pair-partners,
    // and load a de-duplicated version into masterpairdets.
    masterpairdets = {};
    for(i=partnervec.size()-1 ; i>=0 ; i--) {
      detct = partnervec[i].index;
      long matchnum = partnervec[i].lelem;
      if(matchnum==1 && temppairdets[detct].MJD>=0.0) {
	// There were no duplicates, and this detection was not previously
	// marked as a duplicate of some other detection. Load it to masterpairdets.
	masterpairdets.push_back(temppairdets[detct]);
      } else if(matchnum>1 && temppairdets[detct].MJD>=0.0) {
	// This detection was not previously marked as a duplicate of some other
	// detection, but it has some duplicates. Push back the current (central) one.
	masterpairdets.push_back(temppairdets[detct]);
	// Find all possible duplicates for temppairdets[detct], and mark them
	// as duplicates by setting MJD to NEG_MJD_MARK
	j=detct;
	// Decrement j until we are certain to be below any duplicates
	while(j>0 && (temppairdets[j].image > temppairdets[detct].image || (temppairdets[j].image == temppairdets[detct].image && temppairdets[j].Dec >= temppairdets[detct].Dec - deg_matchtol))) j--;
	// Search upward through temppairdets until all duplicates have been found
	while(j<long(temppairdets.size()) && (temppairdets[j].image < temppairdets[detct].image || (temppairdets[j].image == temppairdets[detct].image && temppairdets[j].Dec <= temppairdets[detct].Dec + deg_matchtol))) {
	  if(j!=detct && temppairdets[j].image == temppairdets[detct].image) {
	    // Calculate the distance from detection j to detection detct
	    double dist =  distradec01(temppairdets[detct].RA, temppairdets[detct].Dec, temppairdets[j].RA, temppairdets[j].Dec);
	    if(dist<=deg_matchtol) temppairdets[j].MJD = NEG_MJD_MARK;
	  }
	  j++;
	}
      }
    }
    masterdetnum = masterpairdets.size();
    cout << "The latest pairdets catalog had " << temppairdets.size() << " detections prior to\n";
    cout << "de-duplication, and " << masterdetnum << " detections after\n";
    // Sort the master pairdets vector
    sort(masterpairdets.begin(), masterpairdets.end(), imgsort_hldet());
    
    // Perform index-matching from oldpairdets to masterpairdets.
    olddetnum = oldpairdets.size();
    make_lvec(olddetnum, old_detindex);
    for(detct=0;detct<olddetnum;detct++) {
      // Perform a binary search to identify a match to oldpairdets[detct] in masterpairdets.
      i=0;
      long istep = masterdetnum;
      int search_direction=1;
      while(search_direction!=0 && istep>=1) {
	int foundit=0;
	if(i<0) i=0;
	else if(i>=masterdetnum) i = masterdetnum-1;
	if(masterpairdets[i].image == oldpairdets[detct].image && fabs(masterpairdets[i].Dec - oldpairdets[detct].Dec) <= deg_matchtol) foundit=1;
	else if(masterpairdets[i].image > oldpairdets[detct].image || (masterpairdets[i].image == oldpairdets[detct].image && masterpairdets[i].Dec > oldpairdets[detct].Dec)) {
	  search_direction=-1;
	  istep/=2;
	  if(istep<1) istep=1;
	  i -= istep;
	  if(i<0) i=0;
	} else if(masterpairdets[i].image < oldpairdets[detct].image || (masterpairdets[i].image == oldpairdets[detct].image && masterpairdets[i].Dec < oldpairdets[detct].Dec)) {
	  search_direction=1;
	  istep/=2;
	  if(istep<1) istep=1;
	  i += istep;
	} else {
	  cerr << "Logically forbidden case in binary search\n";
	  cerr << "Images " << oldpairdets[detct].image << " and " << masterpairdets[i].image << "; Dec values " << oldpairdets[detct].Dec << " and " << masterpairdets[i].Dec << "\n";
	}
	if(foundit==1) {
	  // Decrement i until we are certain to be below the match
	  while(i>0 && (masterpairdets[i].image > oldpairdets[detct].image || (masterpairdets[i].image == oldpairdets[detct].image && masterpairdets[i].Dec >= oldpairdets[detct].Dec - deg_matchtol))) i--;
	  // Now search upward until the match is found
	  while(search_direction!=0 && i<masterdetnum && (masterpairdets[i].image < oldpairdets[detct].image || (masterpairdets[i].image == oldpairdets[detct].image && masterpairdets[i].Dec <= oldpairdets[detct].Dec + deg_matchtol))) {
	    // Calculate the distance from detection i to detection detct
	    double dist =  distradec01(oldpairdets[detct].RA, oldpairdets[detct].Dec, masterpairdets[i].RA, masterpairdets[i].Dec);
	    if(dist<=deg_matchtol && masterpairdets[i].image == oldpairdets[detct].image) {
	      old_detindex[detct]=i;
	      search_direction=0;
	    }
	    i++;
	  }
	}
      }
      if(search_direction!=0) {
	cerr << "ERROR: failed to match oldpairdets detection " << detct << " " << oldpairdets[detct].image << " " << oldpairdets[detct].MJD << " " << oldpairdets[detct].RA << " " << oldpairdets[detct].Dec << "\n";
	return(3);
      }
    }
    cout << "Finished index-matching from oldpairdets to masterpairdets.\n";
    // Update detection indices in master_trk2det
    for(i=0;i<long(master_trk2det.size());i++) {
      detct = master_trk2det[i].i2;
      master_trk2det[i].i2 = old_detindex[detct];
    }
 
    // Perform index-matching from inpairdets to masterpairdets.
    indetnum = inpairdets.size();
    make_lvec(indetnum, in_detindex);
    for(detct=0;detct<indetnum;detct++) {
      // Perform a binary search to identify a match to inpairdets[detct] in masterpairdets.
      i=0;
      long istep = masterdetnum;
      int search_direction=1;
      while(search_direction!=0 && istep>=1) {
	int foundit=0;
	if(i<0) i=0;
	else if(i>=masterdetnum) i = masterdetnum-1;
	if(masterpairdets[i].image == inpairdets[detct].image && fabs(masterpairdets[i].Dec - inpairdets[detct].Dec) <= deg_matchtol) foundit=1;
	else if(masterpairdets[i].image > inpairdets[detct].image || (masterpairdets[i].image == inpairdets[detct].image && masterpairdets[i].Dec > inpairdets[detct].Dec)) {
	  search_direction=-1;
	  istep/=2;
	  if(istep<1) istep=1;
	  i -= istep;
	  if(i<0) i=0;
	} else if(masterpairdets[i].image < inpairdets[detct].image || (masterpairdets[i].image == inpairdets[detct].image && masterpairdets[i].Dec < inpairdets[detct].Dec)) {
	  search_direction=1;
	  istep/=2;
	  if(istep<1) istep=1;
	  i += istep;
	} else {
	  cerr << "Logically forbidden case in binary search\n";
	  cerr << "Images " << inpairdets[detct].image << " and " << masterpairdets[i].image << "; Dec values " << inpairdets[detct].Dec << " and " << masterpairdets[i].Dec << "\n";
	}
	if(foundit==1) {
	  // Decrement i until we are certain to be below the match
	  while(i>0 && (masterpairdets[i].image > inpairdets[detct].image || (masterpairdets[i].image == inpairdets[detct].image && masterpairdets[i].Dec >= inpairdets[detct].Dec - deg_matchtol))) i--;
	  // Now search upward until the match is found
	  while(search_direction!=0 && i<masterdetnum && (masterpairdets[i].image < inpairdets[detct].image || (masterpairdets[i].image == inpairdets[detct].image && masterpairdets[i].Dec <= inpairdets[detct].Dec + deg_matchtol))) {
	    // Calculate the distance from detection i to detection detct
	    double dist =  distradec01(inpairdets[detct].RA, inpairdets[detct].Dec, masterpairdets[i].RA, masterpairdets[i].Dec);
	    if(dist<=deg_matchtol && masterpairdets[i].image == inpairdets[detct].image) {
	      in_detindex[detct]=i;
	      search_direction=0;
	    }
	    i++;
	  }
	}
      }
      if(search_direction!=0) {
	cerr << "ERROR: failed to match inpairdets detection " << detct << " " << inpairdets[detct].image << " " << inpairdets[detct].MJD << " " << inpairdets[detct].RA << " " << inpairdets[detct].Dec << "\n";
	return(3);
      }
    }
    cout << "Finished index-matching from inpairdets to masterpairdets.\n";
    // Update detection indices in intrk2det
    for(i=0;i<long(intrk2det.size());i++) {
      detct = intrk2det[i].i2;
      intrk2det[i].i2 = in_detindex[detct];
    }
    // Update tracklet indices in intracklets and intrk2det
    for(trkct=0;trkct<in_trknum;trkct++) intracklets[trkct].trk_ID += master_trknum;
    for(i=0;i<long(intrk2det.size());i++) intrk2det[i].i1 += master_trknum;
    // Append intracklets to mastertracklets
    for(trkct=0;trkct<in_trknum;trkct++) mastertracklets.push_back(intracklets[trkct]);
    // Append intrk2det to master_trk2det
    for(i=0;i<long(intrk2det.size());i++) master_trk2det.push_back(intrk2det[i]);

    cout << "Finished adding tracklet set " << filect << "\n";
    cout << "Total images now: " << masterimage.size() << "\n";
    cout << "Total detections now: " << masterpairdets.size() << "\n";
    cout << "Total tracklets now: " << mastertracklets.size() << "\n";
    cout << "Length of trk2det vector now : " << master_trk2det.size() << "\n";
  }

  // DE_DUPLICATE TRACKLETS
  // Hash all of the tracklets
  intracklets = mastertracklets;
  intrk2det = master_trk2det;
  tracklet_hash = {};
  i=0;
  for(trkct=0;trkct<long(intracklets.size());trkct++) {
    if(intracklets[trkct].trk_ID != trkct) {
      cerr << "ERROR: tracklet index mismatch " << trkct << " vs " << intracklets[trkct].trk_ID << "\n";
      return(4);
    }
    trkvec={};
    while(i<long(intrk2det.size()) && intrk2det[i].i1==trkct) {
      trkvec.push_back(intrk2det[i].i2);
      i++;
    }
    if(trkvec.size()>1) {
      // Hash trkvec
      long trkhash = blend_vector(trkvec);
      one_lindex = long_index(trkhash,trkct);
      tracklet_hash.push_back(one_lindex);
    } else {
      cerr << "ERROR: tracklet " << trkct << " found to have only " << trkvec.size() << " points\n";
      return(5);
    }
  }
  if(intracklets.size() != tracklet_hash.size()) {
    cerr << "ERROR: tracklet_hash vector length doesn't match tracklets: " << tracklet_hash.size() << " vs " << intracklets.size() << "\n";
    return(6);
  }
  // Sort the trkvec hashes
  sort(tracklet_hash.begin(), tracklet_hash.end(), lower_long_index());
  mastertracklets = {};
  master_trk2det = {};
  // Load the first tracklet
  i=0;
  trkct = tracklet_hash[i].index;
  trkvec = tracklet_lookup(intrk2det, trkct);
  for(j=0;j<long(trkvec.size());j++) {
    onepair = longpair(mastertracklets.size(),trkvec[j]);
    master_trk2det.push_back(onepair);
  }
  intracklets[trkct].trk_ID = mastertracklets.size();
  mastertracklets.push_back(intracklets[trkct]);
  for(i=1;i<long(intracklets.size());i++) {
    if(tracklet_hash[i].lelem != tracklet_hash[i-1].lelem) {
      // This tracklet is not a duplicate. Add it to the master output vectors
      trkct = tracklet_hash[i].index;
      trkvec = tracklet_lookup(intrk2det, trkct);
      for(j=0;j<long(trkvec.size());j++) {
	onepair = longpair(mastertracklets.size(),trkvec[j]);
	master_trk2det.push_back(onepair);
      }
      intracklets[trkct].trk_ID = mastertracklets.size();
      mastertracklets.push_back(intracklets[trkct]);
    }
  }
  cout << "De-duplicated tracklets total: " << mastertracklets.size() << "\n";
  cout << "De-duplicated trk2det total: " << master_trk2det.size() << "\n";
  
  // Write image file
  cout << "Writing output image catalog " << imagefile_out << " with " << masterimage.size() << " lines\n";
  outstream1.open(imagefile_out);
  for(masterimct=0;masterimct<long(masterimage.size());masterimct++) {
    outstream1 << fixed << setprecision(8) << masterimage[masterimct].MJD << " " << masterimage[masterimct].RA;
    outstream1 << fixed << setprecision(8) << " " << masterimage[masterimct].Dec << " " << masterimage[masterimct].obscode << " ";
    outstream1 << fixed << setprecision(1) << masterimage[masterimct].X << " " << masterimage[masterimct].Y << " " << masterimage[masterimct].Z << " ";
    outstream1 << fixed << setprecision(4) << masterimage[masterimct].VX << " " << masterimage[masterimct].VY << " " << masterimage[masterimct].VZ << " ";
    outstream1 << masterimage[masterimct].startind << " " << masterimage[masterimct].endind << " " << masterimage[masterimct].exptime << "\n";
  }
  outstream1.close();

  // Write paired detection file
  cout << "Writing paired detection file " << pairdetfile_out << " with " << masterpairdets.size() << " lines\n";
  outstream1.open(pairdetfile_out);
  outstream1 << "#MJD,RA,Dec,mag,trail_len,trail_PA,sigmag,sig_across,sig_along,image,idstring,band,obscode,known_obj,det_qual,origindex\n";
  for(i=0;i<long(masterpairdets.size());i++) {
    outstream1 << fixed << setprecision(7) << masterpairdets[i].MJD << "," << masterpairdets[i].RA << "," << masterpairdets[i].Dec << ",";
    outstream1 << fixed << setprecision(4) << masterpairdets[i].mag << ",";
    outstream1 << fixed << setprecision(2) << masterpairdets[i].trail_len << "," << masterpairdets[i].trail_PA << ",";
    outstream1 << fixed << setprecision(4) << masterpairdets[i].sigmag << ",";
    outstream1 << fixed << setprecision(3) << masterpairdets[i].sig_across << "," << masterpairdets[i].sig_along << ",";
    outstream1 << masterpairdets[i].image << "," << masterpairdets[i].idstring << "," << masterpairdets[i].band << ",";
    outstream1 << masterpairdets[i].obscode << "," << masterpairdets[i].known_obj << ","; 
    outstream1 << masterpairdets[i].det_qual << "," << masterpairdets[i].index << "\n"; 
  }
  outstream1.close();

  // Write tracklet file
  cout << "Writing tracklet file " << trackletfile_out << " with " << mastertracklets.size() << " lines\n";
  outstream1.open(trackletfile_out);
  outstream1 << "#Image1,RA1,Dec1,Image2,RA2,Dec2,npts,trk_ID\n";
  for(i=0;i<long(mastertracklets.size());i++) {
    outstream1 << fixed << setprecision(7) << mastertracklets[i].Img1 << "," << mastertracklets[i].RA1 << "," << mastertracklets[i].Dec1 << ",";
    outstream1 << fixed << setprecision(7) << mastertracklets[i].Img2 << "," << mastertracklets[i].RA2 << "," << mastertracklets[i].Dec2 << ",";
    outstream1 << mastertracklets[i].npts << "," << mastertracklets[i].trk_ID << "\n"; 
  }
  outstream1.close();

  // Write trk2det file
  cout << "Writing trk2det file " << trk2detfile_out << " with " << master_trk2det.size() << " lines\n";
  outstream1.open(trk2detfile_out);
  outstream1 << "#trk_ID,detnum\n";
  for(i=0;i<long(master_trk2det.size());i++) {
    outstream1 << master_trk2det[i].i1 << "," << master_trk2det[i].i2 << "\n"; 
  }
  outstream1.close();

  return(0);
}
