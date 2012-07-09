// Program to read 
// Last modified 14/09/10.
// Created by Burcin Donmez
{
#include <string.h>
	
  gROOT->Reset();
  gStyle->SetPalette(1);
  
  ifstream InFile;
  TString ReadFile;
  ULong64_t NumberOfLines = 0;
  float X, Y, Depth, Energy;

  // Open data file to read.
  cout << "Data file to read:\n";
  cin >> ReadFile;
 
 
  InFile.open(ReadFile);
 
  //....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
  
  // Create a new Root file.
  TFile *File = new TFile("Atten.root","RECREATE");
  
  // Create Histograms.
  // 1-D Histograms.
  TH1F *h1000 = new TH1F("H1000","Depth",50, 0, 5);     
  TH1F *h1010 = new TH1F("H1010","Energy",1000, 0, 999);     
  TH2F *h2000 = new TH2F("H2000","Depth",50, 0, 10, 50, 0, 10);     
 
	
  // Create new Tree called Event.
  TTree *Event = new TTree("Event","Events Tree");
  
  // Create Branches.
  Event->Branch("Energy", &Energy, "Energy");
  Event->Branch("Depth", &Depth, "Depth");
  Event->Branch("X", &X, "X");
  Event->Branch("Y", &Y, "Y");

  
  //....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
  
  // ... Real work starts here!! ... //
  while (InFile.get() != EOF)
    {
      InFile >> X >> Y >> Depth >> Energy; 

      h1000->Fill(Depth);
      h1010->Fill(Energy);
      h2000->Fill(X+5.0,Y+5.0);
     
      Event->Fill();
      
      NumberOfLines++;
      
      if (!(NumberOfLines % 10000)) 
	{
      cout << "Read " << NumberOfLines << " lines!" << endl;			
	}
    } 
  
  cout << "\n" << NumberOfLines << " events read in total!" << endl;
  
  InFile.close();
  // Write histograms to a root file. 
  File->Write();
  
}

  
  
