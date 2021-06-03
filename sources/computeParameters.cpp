//g++ -o DisplayPulses displayPulse.cpp `root-config --cflags --glibs`

#include <iostream>
#include <fstream>
#include <algorithm>
#include <cstdlib>
#include <string>
#include "TTree.h"
#include "TDirectoryFile.h"
#include "TFile.h"
#include "TTreeIndex.h"
#include "TKey.h"
#include "TH1F.h"
#include "TGraph.h"
#include "TMath.h"
#include "TTree.h"

using namespace std;

TGraph* GetGraph(vector<string> pulse, vector<string> header, bool debug)
{

  double sampling = 40;
  TGraph* graph = new TGraph(pulse.size()-1);
  
    for(int f = 1; f < pulse.size() ; f++) {
        double y;
        double x = f/sampling;
        y = stoi(pulse[f]);
	graph->SetPoint(f-1,x,y);
    }

    if (debug){
    for(int i = 0; i < header.size(); i++){
      cout<<header[i]<<" ";
    }
    cout<<endl;
    }
   
    return graph;
}

double GetParameter(vector<string> pulse, string parameter)
{

  double sampling = 40;
  TGraph* graph = new TGraph(pulse.size()-1);

  double baseline = 0;

    for(int f = 1; f < pulse.size() ; f++) {
        double y;
        double x = f/sampling;
        y = stoi(pulse[f]);
        if (f>99 && f<300) baseline+=y/200;
        graph->SetPoint(f-1,x,y);
    }

    int n=4000;
    double max = TMath::MaxElement(n,graph->GetY());
    double amplitude  = max - baseline;
    //cout<<"Baseline: "<<baseline<<"\t Max: "<<max<<"\t Amplitude: "<<amplitude<<"\t Rise: ";

    int t1 = 360;
    int t2 = 440;

    for(int t = 360; t < 440; t++) {
      if((stod(pulse[t]) - baseline) < 0.1*amplitude) t1 = t;
      if((stod(pulse[t]) - baseline) < 0.9*amplitude){
        t2 = t;
      }else{
        t=440;
      }
    }

    //cout<<" "<<t1<<" "<<t2<<" ";

    double rise = (t2-t1)/40.;

    //cout<<rise<<endl;
    //graph->SetTitle(Form("Amplitude: %f Baseline: %f Risetime: %f",amplitude,baseline,rise));

    double param;
    
    if (parameter == "Amplitude") param = amplitude;
    if (parameter == "Baseline") param = baseline;
    if (parameter == "Risetime") param = 1000*rise;
    
    return param;

}

int main(int argc, char** argv){

  double amplitude, baseline,rise,event;
  bool debug;
  
  TFile* oFile = new TFile(argv[2],"RECREATE");

  TTree* outTree = new TTree("outTree", "GeMPI-2 Waveform Analysis");
  outTree->Branch("Amplitude", &amplitude, "Amplitude/D");
  outTree->Branch("Baseline", &baseline, "Baseline/D");
  outTree->Branch("Risetime", &rise, "Risetime/D");
  outTree->Branch("EventNumber", &event, "EventNumber/D");
  
  fstream iFile;
  iFile.open(argv[1]);

  if(!iFile){
    cout<<"File not found!";
  }
  
  if(argc==4){
    debug = argv[3];
  }else{
    debug =false;
  }
    
  /* varaibili di suporto*/
  string lineNumber, CurrentData;
  vector<string> data;
  vector<string> header;
  vector<string> pulse;

  int pulse_index=0;
  int header_index=18;

  bool isPulse= false;
  bool isHeader= true;
  int counter =0;
  event = 0;
  
  while(1)
    {
      
      iFile>>CurrentData;
      if(debug) cout<<CurrentData<<" ";
      if (iFile.eof()) break;

      if(stoi(CurrentData)==8030 && isHeader) header_index=18;
      
      if(isHeader){
	header.push_back(CurrentData);
	header_index-=1;
	if(debug) cout<<header_index<<endl;;

	if(header_index==0){
	  isHeader=false;
	  isPulse=true;
	  header_index=12;
	}
	continue;
      }

      if(isPulse){
        pulse.push_back(CurrentData);
        pulse_index+=1;
      
	if(pulse_index==4000){

	  TGraph* graph = GetGraph(pulse,header,debug);
          oFile->WriteTObject(graph);

	  amplitude = GetParameter(pulse,"Amplitude");
	  baseline = GetParameter(pulse,"Baseline");
          rise = GetParameter(pulse,"Risetime");
	  event+=1;
	  outTree->Fill();
	  cout<<"Event: "<<event<<"\t Baseline: "<<baseline<<"\t Amplitude: "<<amplitude<<"\t Rise: "<<rise<<endl; 
	  pulse.clear();
          header.clear();
          pulse_index=0;
	  isPulse=false;
	  isHeader=true;
        }
	continue;
      }
  
    }
     
  iFile.close();
  outTree->Write();
  oFile->Close();

}
