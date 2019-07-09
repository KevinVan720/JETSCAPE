#include <iostream>
#include <fstream>
#include <memory>
#include <chrono>
#include <thread>

#include "gzstream.h"
#include "PartonShower.h"
#include "JetScapeLogger.h"
#include "JetScapeReader.h"
#include "JetScapeBanner.h"
#include "fjcore.hh"

#include <GTL/dfs.h>

using namespace std;
//using namespace fjcore;

using namespace Jetscape;

template <typename T>
std::string to_string_with_precision(const T a_value, const int n = 2)
{
	std::ostringstream out;
	out << std::setprecision(n) << a_value;
	return out.str();
}

void jet_ana(double jetR, double ymin, double ymax)
{
	int n_events=100000;
	int count=0;
	int eCM = 7000;
	std::fstream fs;
	fs.open("jet_"+to_string_with_precision(eCM)+"_"+to_string_with_precision(jetR)+"_"+to_string_with_precision(ymin)+"_"+to_string_with_precision(ymax)+".dat", std::fstream::out);
	std::fstream fs2;
	fs2.open("bjet_"+to_string_with_precision(eCM)+"_"+to_string_with_precision(jetR)+"_"+to_string_with_precision(ymin)+"_"+to_string_with_precision(ymax)+".dat", std::fstream::out);
	std::fstream fs3;
	fs3.open("cjet_"+to_string_with_precision(eCM)+"_"+to_string_with_precision(jetR)+"_"+to_string_with_precision(ymin)+"_"+to_string_with_precision(ymax)+".dat", std::fstream::out);
	vector<double> pTHatBin({5, 15, 30, 40, 50, 60, 70, 80, 90, 100, 110, 120, 130, 140, 150, 160, 180, 220, 250});
	size_t nBin = pTHatBin.size() - 1;

	int power = -1; // power=-1: anti-kT
	double jetRadius = jetR, jetpTMin = 2.;
	double jetyMin = ymin, jetyMax = ymax;

	std::vector<double> pTBin({18, 21, 24, 28, 32, 37, 43, 49, 56, 64, 74, 84, 97, 114, 133, 153, 174, 196});
	std::vector<double> jet_cs(pTBin.size(), 0.);
	std::vector<double> err(pTBin.size(), 0.);
	std::vector<double> sqSum(pTBin.size(), 0.);
	std::vector<double> bjet_cs(pTBin.size(), 0.);
	std::vector<double> berr(pTBin.size(), 0.);
	std::vector<double> bsqSum(pTBin.size(), 0.);
	std::vector<double> cjet_cs(pTBin.size(), 0.);
	std::vector<double> cerr(pTBin.size(), 0.);
	std::vector<double> csqSum(pTBin.size(), 0.);

	fjcore::JetDefinition jetDef(fjcore::genkt_algorithm, jetRadius, power);
	std::vector<fjcore::PseudoJet> fjInputs;
	fjcore::Selector select_rapidity = fjcore::SelectorRapRange(jetyMin, jetyMax);

	for (int i = 0; i < nBin; i++)
	{
		count=0;
		std::vector <int> jet_ct = std::vector<int>(pTBin.size() - 1, 0);
		std::vector <int> cjet_ct = std::vector<int>(pTBin.size() - 1, 0);
		std::vector <int> bjet_ct = std::vector<int>(pTBin.size() - 1, 0);
		std::string file_name = "(" + to_string_with_precision(pTHatBin[i], 3) + ", " + to_string_with_precision(pTHatBin[i + 1], 3) + ")test_out.dat";
		auto reader = make_shared<JetScapeReaderAscii>(file_name);
		while (!reader->Finished() && count<n_events)
		{
			count++;
			fjInputs.resize(0);
			reader->Next();

			//fjInputs=reader->GetHadronsForFastJet();
			std::vector<shared_ptr<PartonShower>> mShowers;
			mShowers = reader->GetPartonShowers();
			for (int i = 0; i < mShowers.size(); i++)
			{
				std::vector<fjcore::PseudoJet> temp = mShowers[i]->GetFinalPartonsForFastJet();
				for (int j = 0; j < temp.size(); j++)
				{
					fjInputs.push_back(temp[j]);
				}
			}

			//cout<<"parton number: "<<fjInputs.size()<<endl;
			double sigmapb_weight = reader->GetSigmaGen()* 1.0e9/n_events;
			vector<fjcore::PseudoJet> inclusiveJets, jets;
			fjcore::ClusterSequence clustSeq(fjInputs, jetDef);
			inclusiveJets = clustSeq.inclusive_jets(jetpTMin);
			std::vector<fjcore::PseudoJet> selected_jets = select_rapidity(inclusiveJets);
			jets = sorted_by_pt(selected_jets);
			for (int k = 0; k < jets.size(); k++)
			{
				vector<fjcore::PseudoJet> constituents = jets[k].constituents();
				
				for (int p = 0; p < constituents.size(); p++)
				{
					if (constituents[p].has_user_info<HQInfoBase>() 
					&& constituents[p].user_info<HQInfoBase>().hq_channel()==5)
					{
						
						for (unsigned int j = 0; j < pTBin.size() - 1; j++)
						{
							if (jets[k].pt() >= pTBin[j] && jets[k].pt() < pTBin[j + 1])
							{
								bjet_cs[j] +=  sigmapb_weight;
								bsqSum[j] +=  pow(sigmapb_weight, 2);
								break;
							}
						}
						break;
					}
					if (constituents[p].has_user_info<HQInfoBase>() 
					&& constituents[p].user_info<HQInfoBase>().hq_channel()==4)
					{
						
						for (unsigned int j = 0; j < pTBin.size() - 1; j++)
						{
							if (jets[k].pt() >= pTBin[j] && jets[k].pt() < pTBin[j + 1])
							{
								cjet_cs[j] +=  sigmapb_weight;
								csqSum[j] +=  pow(sigmapb_weight, 2);
								break;
							}
						}
						break;
					}
				}
				for (unsigned int j = 0; j < pTBin.size() - 1; j++)
				{
					if (jets[k].pt() >= pTBin[j] && jets[k].pt() < pTBin[j + 1])
					{
						jet_cs[j] += sigmapb_weight;
						sqSum[j] += pow(sigmapb_weight, 2);
						break;
					}
				}
			}

			
		}
		reader->Close();
	}
	for (unsigned int j = 0; j < pTBin.size() - 1; j++)
	{
		err[j] = jet_cs[j] / sqrt(pow(jet_cs[j], 2) / sqSum[j]);
		cerr[j] = cjet_cs[j] / sqrt(pow(cjet_cs[j], 2) / csqSum[j]);
		berr[j] = bjet_cs[j] / sqrt(pow(bjet_cs[j], 2) / bsqSum[j]);
		cout << (pTBin[j] + pTBin[j + 1]) / 2 << " " << jet_cs[j] / (pTBin[j + 1] - pTBin[j]) << " " << err[j] / (pTBin[j + 1] - pTBin[j]) << endl;
		fs << (pTBin[j] + pTBin[j + 1]) / 2 << " " << jet_cs[j] / (pTBin[j + 1] - pTBin[j]) << " " << err[j] / (pTBin[j + 1] - pTBin[j]) << endl;
		fs2 << (pTBin[j] + pTBin[j + 1]) / 2 << " " << bjet_cs[j] / (pTBin[j + 1] - pTBin[j]) << " " << berr[j] / (pTBin[j + 1] - pTBin[j]) << endl;
		fs3 << (pTBin[j] + pTBin[j + 1]) / 2 << " " << cjet_cs[j] / (pTBin[j + 1] - pTBin[j]) << " " << cerr[j] / (pTBin[j + 1] - pTBin[j]) << endl;
	}

	fs.close();
	fs2.close();
	fs3.close();
}

int main(int argc, char **argv)
{
	jet_ana(0.5, 0, 0.5);
	jet_ana(0.5, 0.5, 1.0);
	jet_ana(0.5, 1.0, 1.5);
	jet_ana(0.5, 1.5, 2.0);
	jet_ana(0.5, 2.0, 2.2);

	return 0;
}