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


    std::map<TString, double> McNormPerYear_X3872{ 
       {"2016preVFP", 0.048},
       {"2016postVFP", 0.039},
       {"2017", 0.049},
       {"2018", 0.072}
    };

    std::map<TString, double> McNormPerYear_Psi2S{ 
       {"2016preVFP", 1.41},
       {"2016postVFP", 1.21},
       {"2017", 1.50},
       {"2018", 2.25}
    };

}

#endif 
