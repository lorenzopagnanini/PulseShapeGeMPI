#include <atomic>
#include <chrono>
#include <csignal>
#include <ctime>
#include <iomanip>
#include <iostream>
#include <iterator>
#include <fstream>
#include <sstream>
#include <vector>
#include <math.h>

#include <chrono>
#include <ctime>

// ROOT dependencies

#include "TGraph.h"
#include "TRint.h"
#include "TTreeReader.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TCanvas.h"
#include "TSystem.h"
#include "TSpectrum.h"
#include "TPaveText.h"
#include "TFile.h"
#include "TParameter.h"
#include "TMath.h"
#include "TF1.h"
#include "TArrayL64.h"

//*** function which describes the Pulse ***
Double_t fit_func(Double_t *x, Double_t *par)
{
    Double_t result;
    Double_t A0 = par[0];
    Double_t t0 = par[1];
    Double_t risetime = par[2];
    Double_t timeconst = par[3];
    Double_t P0 = par[4];
    
    result = (A0*(TMath::TanH((x[0]-t0)/risetime)+1))*TMath::Exp(-timeconst*(x[0]-t0))+P0;
    
    return result;
}

//*** function which gives the postion i of the maxium value of a vector
Int_t maximum_value(Int_t _array_length, std::vector<double> Sample)
{
    Double_t maximum=0;
    Int_t maximum_i=0;
    for(Int_t j=0;j<_array_length;j++){
        if(Sample[j]>maximum){
            maximum=Sample[j];
            maximum_i=j;
        }
        
    }
    return maximum_i;
}

//*** function to find peaks in a TGraph, returns n_peaks ***
std::vector<Double_t> get_starting_values(Int_t _array_length, std::vector<Double_t> Sample, TGraph *graph){
    
    // CAEN will always only record a certain time before trigger -> triggered event will be at the beginning -> only search for peaks in the first 2/3 of the samples
    // get the maximal value of the Sample
    _array_length=Int_t(_array_length*(0.66));
    Int_t max_val_i=maximum_value(_array_length,Sample);
    Double_t max_val=Sample[max_val_i];
    
    // roughly estimate the baseline
    Double_t baseline=0;
    for(Int_t i=0; i<50;i++){baseline=baseline+Sample[i];}
    baseline=baseline/50;
    
    // search for the postion i of t0
    Int_t t0_i, count=0;
    while(count<max_val_i){
        if(Sample[count]<((max_val-baseline)/2+baseline)){t0_i=count;}
        else{}
        ++count;
    }

    Double_t high_x, high_y, t0_x, t0_y;
    graph->GetPoint(max_val_i,high_x,high_y);
    graph->GetPoint(t0_i,t0_x,t0_y);
    
    // setting the starting values for the fit
    std::vector<Double_t> starting_values(5,0);
    starting_values[0]=0.5*(max_val-baseline);
    starting_values[1]=t0_x;
    starting_values[2]=high_x-t0_x;
    starting_values[3]=0.00014;
    starting_values[4]=baseline;


    return starting_values;
}

std::vector<Double_t> fitting(std::vector<Double_t> starting_values, Int_t _x_range, TGraph *fit_Pulse)
{   
    bool converge = false;
    // setting the starting values
    std::vector<Double_t> fit_parameter(10,0);
    Double_t starting_values_arr[5];
       for(Int_t i=0; i<5;i++){starting_values_arr[i]=starting_values[i];}
    
    // the fitting process
    Int_t n_runs=0, fitStatus;
    while(!converge){
       TF1 *func = new TF1("func",fit_func,0,_x_range,5);
       func->SetParameters(starting_values_arr);
       func->SetNpx(5000);
       
       fitStatus = fit_Pulse->Fit("func","Q"); 
       
       TF1 *fitted_func = fit_Pulse->GetFunction("func");
       
       for (Int_t p=0; p<5;p++){
           fit_parameter[2*p]=fitted_func->GetParameter(p);
           fit_parameter[2*p+1]=fitted_func->GetParError(p);
       }
       
       if(TMath::IsNaN(fit_parameter[1])==1 || fitStatus==4 || fit_parameter[1]>fit_parameter[0]*10 || fit_parameter[2]>50000){converge=false;}
       else{converge=true;}
       if(n_runs>10){
           for(Int_t i=0;i<2*5;i++){fit_parameter[i]=-2;}
           converge=true;
           std::cout << "not converged after 10 runs" << std::endl; 
       } 
       ++n_runs;
       delete func;
       delete fitted_func;
    }
    std::cout << "fit status: " << fitStatus << std::endl;
    std::cout<< "fit complete" << std::endl;
    return fit_parameter;
}

void fit_one_file(std::string _path, std::string _save_path,std::string TreeName, ULong64_t File ,Int_t energy_cut,UShort_t _sampling_ns = 10)
{
    TFile fileIn(_path.data());
    
//     std::string TreeName = "Data";
    TTreeReader theReader(TreeName.data(), &fileIn);
   
    TTreeReaderValue<UShort_t> Energy(theReader, "Energy");
    TTreeReaderValue<ULong64_t> Timestamp(theReader, "Timestamp");
    TTreeReaderValue<UShort_t> Board(theReader, "Board");
    TTreeReaderValue<UShort_t> Channel(theReader, "Channel");
    TTreeReaderValue<UInt_t> Flags(theReader, "Flags");
    TTreeReaderValue<TArrayS> Samples(theReader, "Samples");

    std::vector<Double_t> fit_para_risetime;
    ULong64_t event_time;
    UShort_t energy_uncalib;
    UShort_t board_id;
    UShort_t channel_number;
    UInt_t flags_num;
    std::vector<UShort_t> event_ID;
    TFile *fileOut = TFile::Open(_save_path.data(),"RECREATE");
    TTree *tree = new TTree("Data","Risetime");    
    tree-> Branch("Risetime",&fit_para_risetime);
    TBranch *Time = tree-> Branch("Timestamp",&event_time);
    TBranch *Board_ID = tree-> Branch("Board",&board_id);
    TBranch *Channel_num = tree-> Branch("Channel",&channel_number);
    TBranch *Flags_num = tree-> Branch("Flags",&flags_num);
    TBranch *Energy_uncalibrated = tree-> Branch("Energy",&energy_uncalib);
    TBranch *Event_ID = tree-> Branch("Event_ID",&event_ID);
    
    std::cout << "Number of Events: " << theReader.GetEntries(1) << std::endl; 
    
    ULong64_t Event = 0;
    while(theReader.Next()){
        
        event_ID.resize(2);
        event_ID[0]=File;
        event_ID[1]=Event;
        
        std::cout << Event << ": ";
        TGraph * Pulse = nullptr;
        // fill the corresponding vectors with the Data points
        Int_t n = (*Samples).GetSize();
        std::vector<double> SampleX(n,0), SampleY(n,0); //creates a vector with n entries and put in each entry a 0
        for (Int_t i = 0 ; i < n ; ++i ){
            SampleX[i] = i*_sampling_ns; // alle 10 ns ein datenpunkt gespeichert
            SampleY[i] = (*Samples).At(i);
        }
        // Create a graph from the vectors
        Pulse = new TGraph(n, SampleX.data(), SampleY.data());
        Int_t maximum_pos=maximum_value(n,SampleY);
        Double_t maximum = SampleY[maximum_pos];
        std::vector<Double_t> fit_para;
    
        // energy_cut: energy cut in ADC channels(only fit events with higher energy
        Int_t begin_cut= 1000; // to avoid events with the highst point at the begining
        Int_t end_cut= 8000; // need to be changed if pretrigger changes
        if(*Energy>energy_cut && maximum_pos*_sampling_ns>begin_cut){
            std::vector<Double_t> par = get_starting_values(n,SampleY,Pulse);
            fit_para = fitting(par,n*_sampling_ns,Pulse);
            std::cout << "fitting" << std::endl;
        }
        else{
            fit_para.clear();
            fit_para.resize(10);
            for(Int_t i=0; i<10; i++){fit_para[i]=-1;}
            std::cout << "below threshold or no peak" << std::endl;
        }
        

        fit_para_risetime.resize(10);
        for (Int_t i=0; i<10; ++i){
            fit_para_risetime[i]=fit_para[i];
        }
        Event=Event+1;
        
        event_time=*Timestamp;
    
        energy_uncalib=*Energy;
        
        board_id=*Board;
        
        channel_number=*Channel;
        
        flags_num=*Flags;

        tree->Fill();

        delete Pulse;
    }
    
    fileIn.Close(); // close "reading" file
    std::cout<<" ----- END READING " << _path << " -----"<<std::endl;

    tree->Write();
    fileOut->Close(); //
    std::cout<<" ----- END WRITING "<< _save_path <<" -----"<<std::endl;    
}

void fit_multiple_files(std::string _path,std::string TreeName, Int_t data_min, Int_t data_max ,Int_t energy_cut,UShort_t _sampling_ns = 10)
{
    ULong64_t file=0;
    for(Int_t i=data_min;i<=data_max;i++){
        std::string path = _path,save_path=_path;
        Int_t length_path=path.length();
        path.erase(length_path-5,5);
        save_path.erase(length_path-5,5);
        Int_t slash_pos=save_path.find_last_of("/");
        std::string begin_of_save_path=save_path.substr(0,slash_pos);
        save_path=save_path.substr(slash_pos+1);
        if(i==0){
            path=path + ".root";
            save_path=save_path+"_risetime" +".root";
        }
        else{
            path=path + Form("_%d",i) + ".root";
            save_path=save_path + Form("_%d",i) + "_risetime" + ".root";
        }
        
        slash_pos=begin_of_save_path.find_last_of("/"); // make sure that the file is saved in the Risetime folder (which is not in UNFILTERED,...)
        begin_of_save_path=begin_of_save_path.substr(0,slash_pos);
        
        save_path=begin_of_save_path+"/Risetime/"+save_path;
        file=i;
        fit_one_file(path,save_path,TreeName,file,energy_cut);
        
        std::cout << "----- END FITTING RISETIME " << Form("File-%d",i) << " -----" << std::endl;
    }
}

void fit_risetime(std::string _path, Int_t data_min, Int_t data_max,Int_t energy_cut)
{
    Int_t Datapos= _path.find("Data");
    std::string data_name=_path.substr(Datapos);
    Int_t CHpos = data_name.find("_CH");
    data_name=data_name.substr(0,CHpos);
    std::string Treename;
    if(data_name=="Data"){Treename="Data";}
    else if(data_name=="DataR"){Treename="Data_R";}
    else if(data_name=="DataF"){Treename="Data_F";}
    else {std::cerr << "------ Problems with the TreeName ------" << std::endl;}
    fit_multiple_files(_path,Treename,data_min,data_max,energy_cut);
}

