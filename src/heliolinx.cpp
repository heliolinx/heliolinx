// March 21, 2023: heliohypy (heliocentric hypothesis code for python)
// Implementation of make_tracklets for python.

#include "solarsyst_dyn_geo01.h"
#include "cmath"

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>
#include <pybind11/stl_bind.h>

namespace py = pybind11;

void fill_struct(hldet & out, hldet const& in) {
    out.MJD = in.MJD;
    out.RA = in.RA;
    out.Dec = in.Dec;
    out.mag = in.mag;
    out.trail_len = in.trail_len;
    out.trail_PA = in.trail_PA;
    out.sigmag = in.sigmag;
    out.sig_across = in.sig_across;
    out.sig_along = in.sig_along;
    out.image = in.image;
    memcpy(out.idstring, in.idstring, sizeof(in.idstring));
    memcpy(out.band, in.band, sizeof(in.band));
    memcpy(out.obscode, in.obscode, sizeof(in.obscode));
    out.known_obj = in.known_obj;
    out.det_qual = in.det_qual;
    out.index = in.index;
}

void fill_struct(EarthState & out, EarthState const& in) {
    out.MJD = in.MJD;
    out.x = in.x;
    out.y = in.y;
    out.z = in.z;
    out.vx = in.vx;
    out.vy = in.vy;
    out.vz = in.vz;
}

void fill_struct(hlimage & out, hlimage const& in) {
    out.MJD = in.MJD;
    out.RA = in.RA;
    out.Dec = in.Dec;
    memcpy(out.obscode, in.obscode, sizeof(in.obscode));
    out.X = in.X;
    out.Y = in.Y;
    out.Z = in.Z;
    out.VX = in.VX;
    out.VY = in.VY;
    out.VZ = in.VZ;
    out.exptime = in.exptime;
}

void fill_struct(tracklet & out, tracklet const& in) {
    out.Img1 = in.Img1;
    out.RA1 = in.RA1;
    out.Dec1 = in.Dec1;
    out.Img2 = in.Img2;
    out.RA2 = in.RA2;
    out.Dec2 = in.Dec2;
    out.npts = in.npts;
    out.trk_ID = in.trk_ID;
}

void fill_struct(longpair & out, longpair const& in) {
    out.i1 = in.i1;
    out.i2 = in.i2;
}

void fill_struct(hlradhyp & out, hlradhyp const& in) {
    out.HelioRad = in.HelioRad;
    out.R_dot = in.R_dot;
    out.R_dubdot = in.R_dubdot;
}

void fill_struct(hlclust & out, hlclust const& in) {
    out.clusternum = in.clusternum;
    out.posRMS = in.posRMS;
    out.velRMS = in.velRMS;
    out.totRMS = in.totRMS;
    out.astromRMS = in.astromRMS;
    out.pairnum = in.pairnum;
    out.timespan = in.timespan;
    out.uniquepoints = in.uniquepoints;
    out.obsnights = in.obsnights;
    out.metric = in.metric;
    memcpy(out.rating, in.rating, sizeof(in.rating));
    out.reference_MJD = in.reference_MJD;
    out.heliohyp0 = in.heliohyp0;
    out.heliohyp1 = in.heliohyp1;
    out.heliohyp2 = in.heliohyp2;
    out.posX = in.posX;
    out.posY = in.posY;
    out.posZ = in.posZ;
    out.velX = in.velX;
    out.velY = in.velY;
    out.velZ = in.velZ;
    out.orbit_a = in.orbit_a;
    out.orbit_e = in.orbit_e;
    out.orbit_MJD = in.orbit_MJD;
    out.orbitX = in.orbitX;
    out.orbitY = in.orbitY;
    out.orbitZ = in.orbitZ;
    out.orbitVX = in.orbitVX;
    out.orbitVY = in.orbitVY;
    out.orbitVZ = in.orbitVZ;
    out.orbit_eval_count = in.orbit_eval_count;
}

void fill_struct(glint_trail & out, glint_trail const& in) {
    out.x  = in.x;
    out.y  = in.y;
    out.length  = in.length;
    out.PA  = in.PA;
    out.linrms  = in.linrms;
    out.eqrms  = in.eqrms;
    out.magmean  = in.magmean;
    out.magrms  = in.magrms;
    out.stepsize  = in.stepsize;
    out.qc1  = in.qc1;
    out.npt  = in.npt;
    out.flashnum  = in.flashnum;
};

void fill_struct(point3d_index & out, point3d_index const& in) {
    out.x  = in.x;
    out.y  = in.y;
    out.z  = in.z;
    out.index  = in.index;
}; 


template<typename T>
std::vector<T> ndarray_to_vec(py::array_t<T> py_vec) {
    py::dtype expected_dtype = py::dtype::of<T>();
    if (!py_vec.dtype().is(expected_dtype)) {
	throw std::runtime_error("Array dtype does not match expected type 'foo'. See mjuric@astro.washington.edu for more details");
    }
    std::vector<T> vec = {};

    // Get a reference to the py_array data
    auto data_ref = py_vec.unchecked();

    // Place the numpy data into the c++ type
    for (long int i = 0; i < data_ref.size(); i++) {
        T data_out;
        auto &data_in = data_ref[i];

        fill_struct(data_out, data_in);

        vec.push_back(data_out);
    }

    return vec;
}

template<typename T>
py::array vec_to_ndarray(std::vector<T> const& vec) {
    // Allocate a structured numpy array of type T
    auto py_vec = py::array_t<T>(vec.size());

    // Get a mutable reference to the ndarray data
    auto data_ref = py_vec.mutable_unchecked();

    // Place vector data into numpy array
    for (long int i = 0; i < data_ref.size(); i++) {
        auto data_in = vec[i];
        auto &data_out = data_ref[i];

        fill_struct(data_out, data_in);
    }
    return py_vec;
}

py::array iotest02(py::array_t<hldet> py_ioin)
{
  std::vector <hldet> iovec = ndarray_to_vec(py_ioin);
  long unsigned int i=0;
  for(i=0; i<iovec.size(); i++) {
    std::cout << iovec[i].MJD << " " << iovec[i].RA << " " << iovec[i].Dec << " " << iovec[i].mag << " " << iovec[i].trail_len << " " << iovec[i].trail_PA << " " << iovec[i].sigmag  << " " << iovec[i].sig_across << " " << iovec[i].sig_along << " " << iovec[i].image << " " << iovec[i].idstring << " " << iovec[i].band << " " << iovec[i].obscode << " " << iovec[i].index << "\n";
    iovec[i].MJD += MJDOFF;
  }

  auto py_ioout = vec_to_ndarray<hldet>(iovec);
  return(py_ioout);
}

py::array observer_coords(double detmjd, double lon, double obscos, double obssine, py::array_t<EarthState> earthin)
{
  std::vector <EarthState> earthpos = ndarray_to_vec(earthin);
  int polyorder=5;
  std::vector <double> posmjd;
  std::vector <point3d> planetpos;
  int i=0;
  int earthnum = earthpos.size();
  point3d earthnow = point3d(0,0,0);
  std::vector <double> outvec;
  
  posmjd={};
  planetpos={};
  for(i=0; i<earthnum; i++) {
    posmjd.push_back(earthpos[i].MJD);
    earthnow = point3d(earthpos[i].x, earthpos[i].y, earthpos[i].z);
    planetpos.push_back(earthnow);
  }  
  observer_barycoords01(detmjd, polyorder, lon, obscos, obssine, posmjd, planetpos, earthnow);
  outvec={};
  outvec.push_back(earthnow.x);
  outvec.push_back(earthnow.y);
  outvec.push_back(earthnow.z);
  
  auto py_outvec = py::array_t<double>(outvec.size());
  auto data_ref = py_outvec.mutable_unchecked();
  // Place vector data into numpy array
  for (long int i = 0; i < data_ref.size(); i++) {
    auto data_in = outvec[i];
    auto &data_out = data_ref[i];
    data_out = data_in;
  }

  return(py_outvec);
}

py::array observer_vel(double detmjd, double lon, double obscos, double obssine, py::array_t<EarthState> earthin)
{

  std::vector <EarthState> earthpos = ndarray_to_vec(earthin);
  int polyorder=5;
  std::vector <double> posmjd;
  std::vector <point3d> planetpos;
  std::vector <point3d> planetvel;
  int i=0;
  int earthnum = earthpos.size();
  point3d earthnow = point3d(0,0,0);
  point3d earthvel = point3d(0,0,0);
  std::vector <double> outvec;

  
  posmjd={};
  planetpos={};
  planetvel={};
  for(i=0; i<earthnum; i++) {
    posmjd.push_back(earthpos[i].MJD);
    earthnow = point3d(earthpos[i].x, earthpos[i].y, earthpos[i].z);
    planetpos.push_back(earthnow);
    earthvel = point3d(earthpos[i].vx, earthpos[i].vy, earthpos[i].vz);
    planetvel.push_back(earthvel);
  } 
  observer_baryvel01(detmjd, polyorder, lon, obscos, obssine, posmjd, planetpos, planetvel, earthnow, earthvel);

  outvec={};
  outvec.push_back(earthnow.x);
  outvec.push_back(earthnow.y);
  outvec.push_back(earthnow.z);
  outvec.push_back(earthvel.x);
  outvec.push_back(earthvel.y);
  outvec.push_back(earthvel.z);
 
  auto py_outvec = py::array_t<double>(outvec.size());
  auto data_ref = py_outvec.mutable_unchecked();
  // Place vector data into numpy array
  for (long int i = 0; i < data_ref.size(); i++) {
    auto data_in = outvec[i];
    auto &data_out = data_ref[i];
    data_out = data_in;
  }

  return(py_outvec);
}

// makeTracklets: April 07, 2023:
// Minimalist wrapper, only handles the python <-> C++ translations.
// All the interesting algorithmic stuff happens in functions called from solarsyst_dyn_geo01.

std::tuple<py::array, py::array, py::array> makeTracklets(
    MakeTrackletsConfig config,
    py::array_t<hldet> py_detvec,
    py::array_t<hlimage> py_imglog
  ) {
  cout << "C++ wrapper for make_tracklets, now fully functional\n";
  
  std::vector <hldet> detvec = ndarray_to_vec(py_detvec);
  std::vector <hlimage> image_log = ndarray_to_vec(py_imglog);
  std::vector <hldet> pairdets;
  std::vector <tracklet> tracklets;
  std::vector <longpair> trk2det;  
  
  make_tracklets3(detvec,image_log,config,pairdets,tracklets,trk2det);
  
  auto py_detout1 = vec_to_ndarray<hldet>(pairdets);
  cout << "loaded pairdets\n";
  auto py_detout2 = vec_to_ndarray<tracklet>(tracklets);
  cout << "loaded tracklets\n";
  auto py_detout3 = vec_to_ndarray<longpair>(trk2det);
  cout << "loaded trk2det\n";

  return(std::make_tuple(py_detout1, py_detout2, py_detout3));
}

// makeTrailedTracklets: May 27, 2025:
// Minimalist wrapper, only handles the python <-> C++ translations.
// All the interesting algorithmic stuff happens in functions called from solarsyst_dyn_geo01.

std::tuple<py::array, py::array, py::array> makeTrailedTracklets(
    MakeTrackletsConfig config,
    py::array_t<hldet> py_detvec,
    py::array_t<hlimage> py_imglog
  ) {
  cout << "C++ wrapper for make_trailed_tracklets, now fully functional\n";
  
  std::vector <hldet> detvec = ndarray_to_vec(py_detvec);
  std::vector <hlimage> image_log = ndarray_to_vec(py_imglog);
  std::vector <hldet> pairdets;
  std::vector <tracklet> tracklets;
  std::vector <longpair> trk2det;  

  // Note: make_trailed_tracklets2 is an advance over make_trailed_tracklets
  // in that it culls the ouput 'pairdets' array down to only those
  // detections (sources) that were actually included in a tracklet.
  make_trailed_tracklets2(detvec,image_log,config,pairdets,tracklets,trk2det);
  
  auto py_detout1 = vec_to_ndarray<hldet>(pairdets);
  cout << "loaded pairdets\n";
  auto py_detout2 = vec_to_ndarray<tracklet>(tracklets);
  cout << "loaded tracklets\n";
  auto py_detout3 = vec_to_ndarray<longpair>(trk2det);
  cout << "loaded trk2det\n";

  return(std::make_tuple(py_detout1, py_detout2, py_detout3));
}

// heliolinc: April 11, 2023:
// Minimalist wrapper, only handles the python <-> C++ translations.
// All the interesting algorithmic stuff happens in functions called from solarsyst_dyn_geo01.

std::tuple<py::array, py::array>heliolinc(
    HeliolincConfig config,
    py::array_t<hlimage> py_imglog,
    py::array_t<hldet> py_detvec,
    py::array_t<tracklet> py_tracklets,
    py::array_t<longpair> py_trk2det,
    py::array_t<hlradhyp> py_radhyp,
    py::array_t<EarthState> earthin
  ) {
  cout << "C++ wrapper for heliolinc\n";
  
  std::vector <hlimage> image_log = ndarray_to_vec(py_imglog);
  std::vector <hldet> detvec = ndarray_to_vec(py_detvec);
  std::vector <tracklet> tracklets = ndarray_to_vec(py_tracklets);
  std::vector <longpair> trk2det = ndarray_to_vec(py_trk2det);
  std::vector <hlradhyp> radhyp = ndarray_to_vec(py_radhyp);
  std::vector <EarthState> earthpos = ndarray_to_vec(earthin);
  int status = 0;
  std::vector <hlclust> outclust;
  std::vector <longpair> clust2det;
     
  status = heliolinc_alg_all(image_log, detvec, tracklets, trk2det, radhyp, earthpos, config, outclust, clust2det);
  if(status!=0) {
    cerr << "ERROR: heliolinc returned failure status " << status << "\n";
    auto py_clustout = vec_to_ndarray<hlclust>({});
    return(std::make_tuple(py_clustout, py_clustout));
  }
      
  auto py_detout1 = vec_to_ndarray<hlclust>(outclust);
  auto py_detout2 = vec_to_ndarray<longpair>(clust2det);

  return(std::make_tuple(py_detout1, py_detout2));
}


// link_Purify: March 08, 2024:
// Minimalist wrapper, only handles the python <-> C++ translations.
// All the interesting algorithmic stuff happens in functions called from solarsyst_dyn_geo01.

std::tuple<py::array, py::array>linkPurify(
    LinkPurifyConfig config,
    py::array_t<hlimage> py_imglog,
    py::array_t<hldet> py_detvec,
    py::array_t<hlclust> py_inclust,
    py::array_t<longpair> py_inclust2det
  ) {
  cout << "C++ wrapper for link_refine_Herget\n";
  
  std::vector <hlimage> image_log = ndarray_to_vec(py_imglog);
  std::vector <hldet> detvec = ndarray_to_vec(py_detvec);
  std::vector <hlclust> inclust = ndarray_to_vec(py_inclust);
  std::vector <longpair> inclust2det = ndarray_to_vec(py_inclust2det);
  int status = 0;
  std::vector <hlclust> outclust;
  std::vector <longpair> outclust2det;
     
  status = link_purify(image_log, detvec, inclust, inclust2det, config, outclust, outclust2det);
  if(status!=0) {
    cerr << "ERROR: link_purify returned failure status " << status << "\n";
    auto py_clustout = vec_to_ndarray<hlclust>({});
    return(std::make_tuple(py_clustout, py_clustout));
  }
      
  auto py_detout1 = vec_to_ndarray<hlclust>(outclust);
  auto py_detout2 = vec_to_ndarray<longpair>(outclust2det);

  return(std::make_tuple(py_detout1, py_detout2));
}

// link_Planarity: March 22, 2024:
// Minimalist wrapper, only handles the python <-> C++ translations.
// All the interesting algorithmic stuff happens in functions called from solarsyst_dyn_geo01.

std::tuple<py::array, py::array>linkPlanarity(
    LinkPurifyConfig config,
    py::array_t<hlimage> py_imglog,
    py::array_t<hldet> py_detvec,
    py::array_t<hlclust> py_inclust,
    py::array_t<longpair> py_inclust2det
  ) {
  cout << "C++ wrapper for link_planarity\n";
  
  std::vector <hlimage> image_log = ndarray_to_vec(py_imglog);
  std::vector <hldet> detvec = ndarray_to_vec(py_detvec);
  std::vector <hlclust> inclust = ndarray_to_vec(py_inclust);
  std::vector <longpair> inclust2det = ndarray_to_vec(py_inclust2det);
  int status = 0;
  std::vector <hlclust> outclust;
  std::vector <longpair> outclust2det;
     
  status = link_planarity(image_log, detvec, inclust, inclust2det, config, outclust, outclust2det);
  if(status!=0) {
    cerr << "ERROR: link_planarity returned failure status " << status << "\n";
    auto py_clustout = vec_to_ndarray<hlclust>({});
    return(std::make_tuple(py_clustout, py_clustout));
  }
      
  auto py_detout1 = vec_to_ndarray<hlclust>(outclust);
  auto py_detout2 = vec_to_ndarray<longpair>(outclust2det);

  return(std::make_tuple(py_detout1, py_detout2));
}


// findGlints: July 17, 2024:
// Minimalist wrapper, only handles the python <-> C++ translations.
// All the interesting algorithmic stuff happens in functions called from solarsyst_dyn_geo01.
// This is the pixel x,y version.

std::tuple<py::array, py::array>findGlints(
    FindGlintsConfig config,
    py::array_t<point3d_index> py_detvec
  ) {
  cout << "C++ wrapper for find_glints_xypix\n";
  
  std::vector <point3d_index> detvec = ndarray_to_vec(py_detvec);
  int status = 0;
  std::vector <glint_trail> trailvec;
  std::vector <longpair> trail2det;
     
  status = find_glints_xypix(detvec, config, trailvec, trail2det);
  if(status!=0) {
    cerr << "ERROR: find_glints_xypix returned failure status " << status << "\n";
    auto py_empty = vec_to_ndarray<point3d_index>({});
    return(std::make_tuple(py_empty, py_empty));
  }
      
  auto py_detout1 = vec_to_ndarray<glint_trail>(trailvec);
  auto py_detout2 = vec_to_ndarray<longpair>(trail2det);

  return(std::make_tuple(py_detout1, py_detout2));
}


// findGlintsRadec: July 17, 2024:
// Minimalist wrapper, only handles the python <-> C++ translations.
// All the interesting algorithmic stuff happens in functions called from solarsyst_dyn_geo01.
// This is the RA, Dec version.
std::tuple<py::array, py::array>findGlintsRadec(
    FindGlintsConfig config,
    py::array_t<point3d_index> py_sourcecat
  ) {
  cout << "C++ wrapper for find_glints_radec\n";
  
  std::vector <point3d_index> detvec = ndarray_to_vec(py_sourcecat);
  int status = 0;
  std::vector <glint_trail> trailvec;
  std::vector <longpair> trail2det;
     
  status = find_glints_radec(detvec, config, trailvec, trail2det);
  if(status!=0) {
    cerr << "ERROR: find_glints_radec returned failure status " << status << "\n";
    auto py_empty = vec_to_ndarray<point3d_index>({});
    return(std::make_tuple(py_empty, py_empty));
  }
      
  auto py_detout1 = vec_to_ndarray<glint_trail>(trailvec);
  auto py_detout2 = vec_to_ndarray<longpair>(trail2det);

  return(std::make_tuple(py_detout1, py_detout2));
}


template <typename S>
py::array_t<S> create_recarray(size_t n) {
    return py::array_t<S>(n);
}
#define NDARRAY_FACTORY(S) m.def("create_" #S, &create_recarray<S>);

PYBIND11_MODULE(heliolinx, m) {
    m.doc() = "pybind11 I/O test"; // optional module docstring
    
    PYBIND11_NUMPY_DTYPE(hldet, MJD, RA, Dec, mag, trail_len, trail_PA, sigmag, sig_across, sig_along, image, idstring, band, obscode, known_obj, det_qual, index);
    PYBIND11_NUMPY_DTYPE(EarthState, MJD, x, y, z, vx, vy, vz);
    PYBIND11_NUMPY_DTYPE(hlimage, MJD, RA, Dec, obscode, X, Y, Z, VX, VY, VZ, startind, endind, exptime);
    PYBIND11_NUMPY_DTYPE(longpair, i1, i2);
    PYBIND11_NUMPY_DTYPE(tracklet, Img1, RA1, Dec1, Img2, RA2, Dec2, npts, trk_ID);
    PYBIND11_NUMPY_DTYPE(hlradhyp, HelioRad, R_dot, R_dubdot);
    PYBIND11_NUMPY_DTYPE(hlclust, clusternum, posRMS, velRMS, totRMS, astromRMS, pairnum, timespan, uniquepoints, obsnights, metric, rating, reference_MJD, heliohyp0, heliohyp1, heliohyp2, posX, posY, posZ, velX, velY, velZ, orbit_a, orbit_e, orbit_MJD, orbitX, orbitY, orbitZ, orbitVX, orbitVY, orbitVZ, orbit_eval_count);
    PYBIND11_NUMPY_DTYPE(point3d_index, x, y, z, index);
    PYBIND11_NUMPY_DTYPE(glint_trail, x, y, length, PA, linrms, eqrms, magmean, magrms, stepsize, qc1, npt, flashnum);
		     
    NDARRAY_FACTORY(hldet)
    NDARRAY_FACTORY(EarthState)
    NDARRAY_FACTORY(hlimage)
    NDARRAY_FACTORY(longpair)
    NDARRAY_FACTORY(tracklet)
    NDARRAY_FACTORY(hlradhyp)
    NDARRAY_FACTORY(hlclust)
    NDARRAY_FACTORY(point3d_index)
    NDARRAY_FACTORY(glint_trail)

/*    py::class_<hldet>(m, "hldet")
      .def(py::init<double &, double &, double &, float &, float &, float &, float &, float &, float &, int &, std::string &, std::string &, std::string &, long &, long &, long &>());
    
    py::class_<EarthState>(m, "EarthState")
      .def(py::init<double &, double &, double &, double &, double &, double &, double &>());
    
    py::class_<hlimage>(m, "hlimage")
      .def(py::init<double &, double &, double &, std::string &, double &, double &, double &, double &, double &, double &, int &, int &, double &>());

    py::class_<longpair>(m, "longpair")
      .def(py::init<long &, long &>());

    py::class_<tracklet>(m, "tracklet")
      .def(py::init<long &, double &, double &, long &, double &, double &, int &, long &>());

    py::class_<hlradhyp>(m, "hlradhyp")
      .def(py::init<double &, double &, double &>());
    
    py::class_<hlclust>(m, "hlclust")
      .def(py::init<long &, double &, double &, double &, double &, int &, double &, int &, int &, double &, std::string &, double &, double &, double &, double &, double &, double &, double &, double &, double &, double &, double &, double &, double &, double &, double &, double &, double &, double &, long &>());
*/

    // Config class for MakeTracklets
    py::class_<MakeTrackletsConfig>(m, "MakeTrackletsConfig")
      .def(py::init<>())
      .def_readwrite("mintrkpts", &MakeTrackletsConfig::mintrkpts)
      .def_readwrite("imagetimetol", &MakeTrackletsConfig::imagetimetol)
      .def_readwrite("maxvel", &MakeTrackletsConfig::maxvel)
      .def_readwrite("minvel", &MakeTrackletsConfig::minvel)
      .def_readwrite("minarc", &MakeTrackletsConfig::minarc)
      .def_readwrite("maxtime", &MakeTrackletsConfig::maxtime)
      .def_readwrite("mintime", &MakeTrackletsConfig::mintime)
      .def_readwrite("imagerad", &MakeTrackletsConfig::imagerad)
      .def_readwrite("maxgcr", &MakeTrackletsConfig::maxgcr)
      .def_readwrite("exptime", &MakeTrackletsConfig::exptime)
      .def_readwrite("siglenscale", &MakeTrackletsConfig::siglenscale)
      .def_readwrite("sigpascale", &MakeTrackletsConfig::sigpascale)
      .def_readwrite("max_netl", &MakeTrackletsConfig::max_netl)
      .def_readwrite("time_offset", &MakeTrackletsConfig::time_offset)
      .def_readwrite("forcerun", &MakeTrackletsConfig::forcerun)
      .def_readwrite("verbose", &MakeTrackletsConfig::verbose);

    py::class_<HeliolincConfig>(m, "HeliolincConfig")
      .def(py::init<>())
      .def_readwrite("MJDref", &HeliolincConfig::MJDref)
      .def_readwrite("autorun", &HeliolincConfig::autorun)
      .def_readwrite("clustrad", &HeliolincConfig::clustrad) 
      .def_readwrite("clustchangerad", &HeliolincConfig::clustchangerad) 
      .def_readwrite("dbscan_npt", &HeliolincConfig::dbscan_npt) 
      .def_readwrite("minobsnights", &HeliolincConfig::minobsnights) 
      .def_readwrite("mintimespan", &HeliolincConfig::mintimespan)
      .def_readwrite("mingeodist", &HeliolincConfig::mingeodist)
      .def_readwrite("maxgeodist", &HeliolincConfig::maxgeodist)
      .def_readwrite("geologstep", &HeliolincConfig::geologstep)
      .def_readwrite("mingeoobs", &HeliolincConfig::mingeoobs)
      .def_readwrite("minimpactpar", &HeliolincConfig::minimpactpar)
      .def_readwrite("use_univar", &HeliolincConfig::use_univar)
      .def_readwrite("max_v_inf", &HeliolincConfig::max_v_inf)
      .def_readwrite("verbose", &HeliolincConfig::verbose);

    // Config class for LinkRefine    
    py::class_<LinkRefineConfig>(m, "LinkRefineConfig")
      .def(py::init<>())
      .def_readwrite("MJDref", &LinkRefineConfig::MJDref)
      .def_readwrite("simptype", &LinkRefineConfig::simptype) 
      .def_readwrite("ptpow", &LinkRefineConfig::ptpow) 
      .def_readwrite("nightpow", &LinkRefineConfig::nightpow) 
      .def_readwrite("timepow", &LinkRefineConfig::timepow)
      .def_readwrite("rmspow", &LinkRefineConfig::rmspow)
      .def_readwrite("maxrms", &LinkRefineConfig::maxrms)
      .def_readwrite("verbose", &LinkRefineConfig::verbose);
    
    // Config class for LinkPurify
    py::class_<LinkPurifyConfig>(m, "LinkPurifyConfig")
      .def(py::init<>())
      .def_readwrite("useorbMJD", &LinkPurifyConfig::useorbMJD) 
      .def_readwrite("simptype", &LinkPurifyConfig::simptype) 
      .def_readwrite("ptpow", &LinkPurifyConfig::ptpow) 
      .def_readwrite("nightpow", &LinkPurifyConfig::nightpow) 
      .def_readwrite("timepow", &LinkPurifyConfig::timepow)
      .def_readwrite("rmspow", &LinkPurifyConfig::rmspow)
      .def_readwrite("maxrms", &LinkPurifyConfig::maxrms)
      .def_readwrite("max_oop", &LinkPurifyConfig::max_oop)
      .def_readwrite("rejfrac", &LinkPurifyConfig::rejfrac)
      .def_readwrite("maxrejnum", &LinkPurifyConfig::maxrejnum)
      .def_readwrite("max_astrom_rms", &LinkPurifyConfig::max_astrom_rms)
      .def_readwrite("minobsnights", &LinkPurifyConfig::minobsnights) 
      .def_readwrite("minpointnum", &LinkPurifyConfig::minpointnum) 
      .def_readwrite("use_heliovane", &LinkPurifyConfig::use_heliovane) 
      .def_readwrite("verbose", &LinkPurifyConfig::verbose);

    // Config class for FindGlints
    py::class_<FindGlintsConfig>(m, "FindGlintsConfig")
      .def(py::init<>())
      .def_readwrite("minpoints", &FindGlintsConfig::minpoints) 
      .def_readwrite("maxgcr", &FindGlintsConfig::maxgcr)
      .def_readwrite("maxrange", &FindGlintsConfig::maxrange)
      .def_readwrite("centerknown", &FindGlintsConfig::centerknown)
      .def_readwrite("incenRA", &FindGlintsConfig::incenRA)
      .def_readwrite("incenDec", &FindGlintsConfig::incenDec)
      .def_readwrite("freq_downscale", &FindGlintsConfig::freq_downscale)
      .def_readwrite("freq_upscale", &FindGlintsConfig::freq_upscale)
      .def_readwrite("max_phase_err", &FindGlintsConfig::max_phase_err);

    
    m.def("iotest02", &iotest02, "A function to test python I/O");
    m.def("observer_coords", &observer_coords, "calculate position of an observer on Earth");
    m.def("observer_vel", &observer_vel, "calculate velocity of an observer on Earth");
    m.def("makeTracklets", &makeTracklets,  "Make tracklets from set of detections.");
    m.def("makeTrailedTracklets", &makeTrailedTracklets,  "Make tracklets from set of trailed detections.");
    m.def("heliolinc", &heliolinc,  "Link input tracklets into candidate discoveries.");
    m.def("linkPurify", &linkPurify, "Purify linkages, eliminating duplicates, and rejecting astrometric outliers.");
    m.def("linkPlanarity", &linkPlanarity, "Purify linkages, eliminating duplicates, and rejecting astrometric outliers.");
    m.def("findGlints", &findGlints, "Identify glint trails produced by space junk, using pixel x,y coordinates");
    m.def("findGlintsRadec", &findGlintsRadec, "Identify glint trails produced by space junk, using RA, Dec coordinates");
  }


