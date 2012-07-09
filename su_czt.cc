
#include "G4RunManager.hh"
#include "G4UImanager.hh"
#include "G4UIterminal.hh"
//#include "G4UIWin32.hh"
#include "G4UItcsh.hh"

#include "G4VisExecutive.hh"

#include "ExDetectorConstruction.hh"
#include "ExPhysicsList.hh"
#include "ExPrimaryGeneratorAction.hh"
#include "ExRunAction.hh"
#include "SteppingAction.hh"

//#include "G4UIQt.hh"
//#include "G4Qt.hh"

#include "ExGlobalParameters.hh"
#include <sys/stat.h>
#include <sys/types.h>

int main(int argc, char ** argv)
{
  CGlobalParameters::LoadParameters();
  CGlobalParameters::print();
  system("rm -r output/imaging_output.txt");
  system("rm -r output/electron_clouds.txt");
  system("rm -r output/stepping.txt");
  system("rm -r output/electron_clouds.bin");
  system("rm -r output/pix_energy.txt");

  // Construct the default run manager
  G4RunManager* runManager = new G4RunManager;

  // set mandatory initialization classes
  // detector
  ExDetectorConstruction * detectorConstruction = new ExDetectorConstruction;
  runManager->SetUserInitialization(detectorConstruction);

  // physics
  ExPhysicsList * physicsList = new ExPhysicsList;
  runManager->SetUserInitialization(physicsList);

  // set mandatory user action class
  // source generator
  ExPrimaryGeneratorAction* gen_action = new ExPrimaryGeneratorAction();
  runManager->SetUserAction(gen_action);

  // Arbitrary user action 
  // set run action
  ExRunAction * run_action = new ExRunAction;
  runManager->SetUserAction(run_action);

  //##################################
  SteppingAction * step_action = new SteppingAction;
  runManager->SetUserAction(step_action);

 // Initialize G4 kernel
  runManager->Initialize();
  
  if ( argc > 1 )
  {
    // UIsession intialization
     G4UIsession* session = new G4UIterminal();
        
    // visualization manager
    G4VisManager* visManager = new G4VisExecutive;
    visManager->Initialize();

    // get the pointer to the UI manager
    
    // User can set the verbosities from macro file.

    G4UImanager* UI = G4UImanager::GetUIpointer();
    //UI->ApplyCommand("/run/verbose 0");
    //UI->ApplyCommand("hit/verbose 0");
    //UI->ApplyCommand("/event/verbose 0");
    //UI->ApplyCommand("/process/verbose 0");
    //UI->ApplyCommand("/tracking/verbose 0");
    UI->ApplyCommand("/control/execute circlesource.mac");
    
    session->SessionStart();  
    UI->ApplyCommand("/run/dumpCouples");

    delete session;
    delete visManager;
  }
  else
  {
    //### start a run
    long numberOfEvent = CGlobalParameters::numberOfEvent;
    system("mkdir -p output");
    system("rm -f output/imaging_output.txt");
    system("rm -f output/imaging_output.bin");
    system("rm -f output/electron_clouds*.dat");
    system("rm -f output/init_parameters.txt");
    system("rm -f output/init_parameters.bin");

    runManager->BeamOn(numberOfEvent);


    // job termination
  }
 
  delete runManager;
  return 0;
}
