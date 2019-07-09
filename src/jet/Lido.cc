// -----------------------------------------
// JetScape (modular/task) based framework
// Intial Design: Joern Putschke (2017)
//                (Wayne State University)
// -----------------------------------------
// License and Doxygen-like Documentation to be added ...

// quick and dirty test class implementation for Eloss modules ...
// can be used as a user template ...

#include "Lido.h"
#include "JetScapeLogger.h"
//#include "JetScapeXML.h"
#include <string>
#include <random>
#include <chrono>

#define MAGENTA "\033[35m"

unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
std::default_random_engine generator(seed);
std::normal_distribution<double> ZeroOneDistribution(0.0, 1.0);

//double fmc_to_GeV_m1 = 5.026;

//#include "tinyxml2.h"
#include <iostream>

#include "FluidDynamics.h"

using namespace Jetscape;

using std::ifstream;
using std::ios;
using std::ofstream;
using std::ostream;

Lido::Lido()
{
  SetId("Lido");
  //VERBOSE(8);
}

Lido::~Lido()
{
  VERBOSE(8);
}

void Lido::Init()
{
  JSINFO << "Initialize Lido ...";
  std::string mode = "old";
  std::string path = "settings.xml";
  std::string table_path = "table.h5";
  boost::property_tree::ptree config;
  double mu = 2.0;
  double const_alphas = 0.3;

  initialize(mode, path, table_path, mu, -const_alphas, 0., 0., 0., 0., 0., 0.,
             1., 0.75);

  tinyxml2::XMLElement *eloss= JetScapeXML::Instance()->GetXMLRoot()->FirstChildElement("Eloss" );
  if ( !eloss )     throw std::runtime_error("Eloss not properly initialized in XML file ...");

  Q0=1.0;
  hydro_Tc = 0.16;
  //martini->FirstChildElement("hydro_Tc")->QueryDoubleText(&hydro_Tc);
  hydro_tStart = 0.0;

}

void Lido::DoEnergyLoss(double deltaT, double Time, double Q2, vector<Parton> &pIn, vector<Parton> &pOut)
{

  //VERBOSESHOWER(0)<< MAGENTA << "SentInPartons Signal received : "<<deltaT<<" "<<Q2<<" "<<&pIn;

  // particle info
  int Id, newId;
  FourVector pin;     // 4 vector of incoming parton
  FourVector pinrest; // 4 vector of incoming parton in rest frame of fluid cell
  FourVector pout;
  FourVector poutrest;

  FourVector xin;  // 4 vector for incoming position
  FourVector xout; // 4 vector for outgoing position (for next time step!)
  double eta;      // pseudo-rapidity

  // flow info
  double vx, vy, vz; // 3 components of flow velocity
  double T;          // Temperature of fluid cell

  for (int i = 0; i < pIn.size(); i++)
  {
    
    Id = pIn[i].pid();

    //VERBOSE(8) << "--------------------------------particle id: " << Id << " channel " << pIn[i].user_info<HQInfoBase>().hq_channel() << " mother id: " << pIn[i].user_info<HQInfoBase>().hq_mother_id();

    VERBOSE(0)<<pIn[i].px()<<" "<<pIn[i].py()<<" "<<pIn[i].pz()<<" "<<pIn[i].e()<<" ";    
    pin = FourVector(pIn[i].px(), pIn[i].py(), pIn[i].pz(), pIn[i].e());
    xin = FourVector(pIn[i].x_in().x(), pIn[i].x_in().y(), pIn[i].x_in().z(), Time);

    std::unique_ptr<FluidCellInfo> check_fluid_info_ptr;

    double tt = xin.t();
    double xx = xin.x() + (Time-tt)*pin.x()/pin.t();
    double yy = xin.y() + (Time-tt)*pin.y()/pin.t();
    double zz = xin.z() + (Time-tt)*pin.z()/pin.t();

    GetHydroCellSignal(Time, xx, yy, zz, check_fluid_info_ptr);
    //VERBOSE(0)<<"Inputs: "<<Time<<" "<<xx<<" "<< yy<<" "<< zz;
    //VERBOSE(0) << "Temperature (Signal) = "
    //           << check_fluid_info_ptr->temperature << " " << check_fluid_info_ptr->vx << " " << check_fluid_info_ptr->vy << " " << check_fluid_info_ptr->vz;

    vx = check_fluid_info_ptr->vx;
    vy = check_fluid_info_ptr->vy;
    vz = check_fluid_info_ptr->vz;
    T = check_fluid_info_ptr->temperature;

    //if (pIn[i].t() > Q0*Q0 + 0.00001 || Time <=hydro_tStart || T < hydro_Tc) continue;
    if (pIn[i].t() > Q0*Q0 + 0.0000001) continue;
    TakeResponsibilityFor ( pIn[i] ); // Generate error if another module already has responsibility.

    std::vector<fourvec> FS;
    std::vector<fourvec> FSG;
    particle Lido_particle;
    Lido_particle.mass = pIn[i].restmass();                                                                                        // mass
    Lido_particle.pid = Id;                                                                                                        // charm quark
    Lido_particle.x = fourvec{xin.t() * fmc_to_GeV_m1, xin.x() * fmc_to_GeV_m1, xin.y() * fmc_to_GeV_m1, xin.z() * fmc_to_GeV_m1}; //position
    Lido_particle.p = fourvec{pin.t(), pin.x(), pin.y(), pin.z()};
    Lido_particle.weight = 1.;
    if (pIn[i].has_user_info<LidoParticleInfo>())
    {
      Lido_particle.radlist = pIn[i].user_info<LidoParticleInfo>().radlist();
      Lido_particle.x0 = pIn[i].user_info<LidoParticleInfo>().x0_;
      Lido_particle.p0 = pIn[i].user_info<LidoParticleInfo>().p0_;
      Lido_particle.T0 = pIn[i].user_info<LidoParticleInfo>().T0_;
      Lido_particle.mfp0 = pIn[i].user_info<LidoParticleInfo>().mfp0_;
      Lido_particle.Tf = pIn[i].user_info<LidoParticleInfo>().Tf_;
      Lido_particle.vcell = pIn[i].user_info<LidoParticleInfo>().vcell_;
    }
    else
    {
      vector<particle> rad;
      Lido_particle.radlist = rad;
      Lido_particle.x0 = Lido_particle.x;
      Lido_particle.p0 = Lido_particle.p;
      Lido_particle.Tf=0.;
      Lido_particle.vcell.resize(3);
			Lido_particle.vcell[0] = vx; 
			Lido_particle.vcell[1] = vy; 
			Lido_particle.vcell[2] = vz; 
      Lido_particle.is_vac = false;
			Lido_particle.is_virtual = false;
    }

    //JSINFO<<" position0: "<<Lido_particle.x0.t()<<" position: "<<Lido_particle.x.t()<<" momentum0: "<<Lido_particle.p0.t()<<" momentum: "<<Lido_particle.p.t();
    //JSINFO<<" T0: "<<Lido_particle.T0<<" "<<Lido_particle.mfp0<<" "<<Lido_particle.Tf;
    //JSINFO<<" vcell:" <<Lido_particle.vcell[0]<<" "<<Lido_particle.vcell[1]<<" "<<Lido_particle.vcell[2];
    //JSINFO<<" radlist #:" <<Lido_particle.radlist.size();

    std::vector<particle> Lido_pOut;
    int channel = update_particle_momentum_Lido(deltaT * fmc_to_GeV_m1, T, {vx, vy, vz}, Lido_particle, Lido_pOut);
    if(channel>1)
    {
      JSINFO << "channel: " <<channel; 
      
    } 

    for (int j = 0; j < Lido_pOut.size(); j++)
    {
      particle par = Lido_pOut[j];
      FourVector x = FourVector(par.x.x()/fmc_to_GeV_m1, par.x.y()/fmc_to_GeV_m1, par.x.z()/fmc_to_GeV_m1, par.x.t()/fmc_to_GeV_m1);
      FourVector p = FourVector(par.p.x(), par.p.y(), par.p.z(), par.p.t());
      pOut.push_back(Parton(0, par.pid, 0, p, x));
      pOut[pOut.size()-1].set_user_info(new LidoParticleInfo(par.T0, par.mfp0, par.Tf, 
      par.x0, par.p0, par.vcell, par.radlist));
      pOut[pOut.size()-1].set_form_time(0.);
    }
  }

  //JSINFO <<"pOut #: "<<pOut.size(); 

  return;   
}
