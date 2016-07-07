#include "sigRegion.h"

//////////////////////////////////////////////////////
// Function definitions for generic "selection" class


// Use this constructor if you want to impose a minimum and/or maximum value on your variable
template <class T> selection<T>::selection( const T& (*func)(), T minval, T maxval ) {
  cms3Function = func;
  minimum = minval;
  maximum = maxval;
  equal_val = std::numeric_limits<T>::max(); // Set this to an unlikely value for safety
}

// Use this constructor if you want to require your variable == some value
template <class T> selection<T>::selection( const T& (*func)(), T eqval ) {
  cms3Function = func;
  equal_val = eqval;
  minimum = std::numeric_limits<T>::max(); // Set incompatible min/max values for safety
  maximum = std::numeric_limits<T>::lowest();
}

template <class T> bool selection<T>::Pass() {
  if( cms3Function() > minimum &&
	  cms3Function() <= maximum ) return true;
  if( cms3Function() == equal_val ) return true;
  return false;
}


// Specialized "pass" function for the <int> type

template <> bool selection<int>::Pass() {
  if( cms3Function() >= minimum &&            // Notice that the lower bound is now inclusive
	  cms3Function() <= maximum ) return true;
  if( cms3Function() == equal_val ) return true;
  return false;
}



///////////////////////////////////////////////
// Function definitions for "sigRegion" class

sigRegion::sigRegion( std::string mylabel, std::string niceName )
  : label(mylabel),
	name_table(niceName),
	name_root(niceName)
{}

sigRegion::sigRegion( std::string mylabel, std::string tabName, std::string rootName )
  : label(mylabel),
	name_table(tabName),
	name_root(rootName)
{}

void sigRegion::addSelection( selectionBase* mySelection ) { selections.push_back(mySelection); }

TString sigRegion::GetLabel()     { return static_cast<TString>(label);      }
TString sigRegion::GetTableName() { return static_cast<TString>(name_table); }
TString sigRegion::GetRootName()  { return static_cast<TString>(name_root);  }

bool sigRegion::PassAllCuts() {
  for( selectionBase* sel : selections ) if( !sel->Pass() ) return false;
  return true;
}


