#ifndef CSV_DATA_H_GUARD
#define CSV_DATA_H_GUARD

#include <fstream>
#include <string>

struct csv_data {
  short electron_sector;
  float w;
  float q2;
  float electron_p;
  float electron_theta;
  float electron_phi;
  float proton_p;
  float proton_theta;
  float proton_phi;
  std::string type;

  // Static functions can be called without making a new struct
  static std::string header() {
    // Make a string for the header of the csv file
    return "electron_sector,w_uncorr,q2_uncorr,electron_p,electron_theta,electron_phi,proton_p,proton_theta,proton_"
           "phi,type";
  }

  friend std::ostream &operator<<(std::ostream &os, const csv_data &data) {
    os << std::setprecision(10);
    os << data.electron_sector << ",";
    os << data.w << ",";
    os << data.q2 << ",";
    os << data.electron_p << ",";
    os << data.electron_theta << ",";
    os << data.electron_phi << ",";
    os << data.proton_p << ",";
    os << data.proton_theta << ",";
    os << data.proton_phi << ",";
    os << data.type;

    return os;
  }
};

#endif