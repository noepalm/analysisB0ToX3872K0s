#ifndef LumiConstants_h
#define LumiConstants_h

#include <iostream> 
#include <map> 
#include "TString.h"

namespace LumiConstants{

    std::map<TString, double> LumiPerYear_Run2{ // https://twiki.cern.ch/twiki/bin/viewauth/CMS/PdmVRun2LegacyAnalysis
      {"2016preVFP", 19.50},
      {"2016postVFP", 16.80},
      {"2017", 41.48},
      {"2018", 59.83}
    };

    std::map<TString, double> LumiPerYear_Run3 { //average for each era from https://docs.google.com/spreadsheets/d/1GVdrnFBLlWNC6ZGb8z7qOCGzgcGGb4q5FcpcNqVg2FY/edit#gid=1770029759
      {"2022", 29.608},
      {"2022D", 2.161}, 
      {"2022E", 6.078},
      {"2022F", 18.271},
      {"2022G", 3.098}
    };


    std::map<TString, double> McNormPerYear_X3872{ 
      {"2016preVFP", 0.048},
      {"2016postVFP", 0.039},
      {"2017", 0.049},
      {"2018", 0.072},
      {"2022", 29.608 * 0.0047837 * 4./7.}, // [PLACEHOLDER] TO BE UPDATED!!!
    };

    std::map<TString, double> McNormPerYear_Psi2S{ 
      {"2016preVFP", 1.41},
      {"2016postVFP", 1.21},
      {"2017", 1.50},
      {"2018", 2.25},
      {"2022", 29.608 * 0.1897 * 4./7.},
    };

    // D: 2.161 / 415437. * 5.5e11 * 0.0239 * 5.8e-4 / 2. * 0.3468 * 0.05961  = 0.4099
    // F: 18.271 / 415437. * 5.5e11 * 0.0239 * 5.8e-4 / 2. * 0.3468 * 0.05961 = 3.46588
    // all: 29.608 / 415437. * 5.5e11 * 0.0239 * 5.8e-4 / 2. * 0.3468 * 0.05961 = 5.6164344

    // [data] all: 29.608 / 471983. * 5.426e11 * 0.0334 * 1.1e-4 / 2 * 0.038 * 0.05961
}

#endif 
