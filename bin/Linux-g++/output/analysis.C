#include <string.h>
#include <fstream>
#include <cstdlib>
#include <iostream>

void Help(){
  cout << "********************************************************" << endl;
  cout << "Available Funcitons  :  " << endl;
  cout  << "-void Draw()" << endl;
  cout  << "********************************************************" << endl;
}

void Draw(){     
  ifstream InFile;
  double energy; int en;
  int pix;
  char * name = new char[12];

  TCanvas *canv = new TCanvas("canv","canvas",1);
  TH1F *histo = new TH1F[20];
  canv->Divide(4,5);
  
  for (int i =0; i<20 ; i++ ){
    name = "channel";
    name[7] = " " ;
    name[8] = 49+i;
    histo[i]=new TH1F(name,name,10000,0,140);
  }
  
  for (int num=0; num <20 ; num++){
    InFile.open("pix_energy.txt");
    while (!InFile.eof())
    {
      InFile >> pix >> energy;
      en = energy*1000;
      if(num == pix) histo[num].Fill(energy*1000);
    }
    canv->cd(num+1);
    histo[num].Draw();
    InFile.close();
  }
} 
