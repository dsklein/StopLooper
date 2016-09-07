#include "systematic.h"
#include <iostream>


// Function definitions for "systematic" class


// Constructor
systematic::systematic( TString systName, direction whichDir, double (*func)() )
  : name(systName),
	var_dir(whichDir),
	reweight_func(func)
	// weightVar(NULL)
{
}



// The other handy functions
systematic::direction systematic::GetDir() { return var_dir; }
TString systematic::GetName()  { return name; }
bool systematic::IsUp()        { return (var_dir==kUp || var_dir==kSkipUp); }
bool systematic::IsDown()      { return (var_dir==kDown || var_dir==kSkipDown); }
bool systematic::IsVariation() { return (var_dir == kVariation); }
bool systematic::IsSkip()      { return (var_dir==kSkipUp || var_dir==kSkipDown); }

double systematic::GetWeight() { return reweight_func(); }

TString systematic::GetNameLong() {
  TString longname = name;
  if( IsUp() ) longname += "up";
  else if( IsDown() ) longname += "down";
  return longname;
}
