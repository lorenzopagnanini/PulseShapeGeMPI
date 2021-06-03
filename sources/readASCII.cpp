
//g++ -o ReadASCII readASCII.cpp `root-config --cflags --glibs`

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
#include "TGraph.h"
#include "TH1F.h"

using namespace std;

int main(int argc, char** argv){

  fstream iFile;
  iFile.open(argv[1]);

  if(!iFile){
    cout<<"File not found!";
  }

  ofstream oFile(argv[2]);

  /* varaibili di suporto*/
  string lineNumber;
  int data0,data1,data2,data3,data4,data5,data6,data7;
  vector<string> data;
  vector<string> header;
  vector<string> pulse;

  while(1)
    {

      iFile >> lineNumber >> data0 >> data1 >> data2 >> data3 >> data4 >> data5 >> data6 >> data7; 
      oFile<<data0<<endl;
      oFile<<data1<<endl;
      oFile<<data2<<endl;
      oFile<<data3<<endl;
      oFile<<data4<<endl;
      oFile<<data5<<endl;
      oFile<<data6<<endl;
      oFile<<data7<<endl;
     
      if (iFile.eof()) break;
      
    }
  
  iFile.close();
  oFile.close();
}


