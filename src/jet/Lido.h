// -----------------------------------------
// JetScape (modular/task) based framework
// Intial Design: Joern Putschke (2017)
//                (Wayne State University)
// -----------------------------------------
// License and Doxygen-like Documentation to be added ...

// quick and dirty test class for Eloss modules ...
// can be used as a user template ...


#ifndef LIDO_H
#define LIDO_H

#include <fstream>
#include <math.h>
#include "JetEnergyLossModule.h"
#include "JetScapeConstants.h"
#include "workflow.h"

using namespace Jetscape;

class LidoParticleInfo: public fjcore::PseudoJet::UserInfoBase
  {
  public:
    std::vector<particle> radlist_;
    double T0_, mfp0_, Tf_; // production temperature, local mfp
	  fourvec x0_; // production location
	  fourvec p0_; // production momentum
	  fourvec mother_p_;
    std::vector<double> vcell_;
    
  LidoParticleInfo(double T0, double mfp0, double Tf, fourvec x0, fourvec p0, std::vector<double> vcell, vector<particle> radlist, fourvec mother_p) 
    : T0_(T0), mfp0_(mfp0), Tf_(Tf), x0_(x0), p0_(p0), vcell_(vcell), radlist_(radlist), mother_p_(mother_p) {};
    std::vector<particle> radlist() const {return radlist_;}
  };

class Lido : public JetEnergyLossModule<Lido>
{
 private:

 std::map<int, std::vector<Process> > AllProcesses;
 double hydro_Tc;      // critical temperature
 double hydro_tStart;  // initilization time of hydro
 double Q0;

 public:

  Lido();
  virtual ~Lido();

  //main//
  void Init();
  void DoEnergyLoss(double deltaT, double Time, double Q2, vector<Parton>& pIn, vector<Parton>& pOut);
  void WriteTask(weak_ptr<JetScapeWriter> w) {};
};



#endif
