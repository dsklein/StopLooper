#include "sigRegion.h"

//////////////////////////////////////////////////////
// Function definitions for generic "selection" class


// Constructor:  minimum < tas::someFunction() < maximum
template <class T> selection<T>::selection( const T& (*func)(), T minval, T maxval ) {
  cms3Function = func;
  cutVar = NULL;
  minimum = minval;
  maximum = maxval;
  equal_val = std::numeric_limits<T>::max();
}

// Constructor:  tas::someFunction() == specificValue
template <class T> selection<T>::selection( const T& (*func)(), T eqval ) {
  cms3Function = func;
  cutVar = NULL;
  equal_val = eqval;
  minimum = std::numeric_limits<T>::max();
  maximum = std::numeric_limits<T>::lowest();
}

// Constructor:  minimum < someGlobalVariable < maximum
template <class T> selection<T>::selection( T* myCutVar, T minval, T maxval) {
  cutVar = myCutVar;
  cms3Function = NULL;
  minimum = minval;
  maximum = maxval;
  equal_val = std::numeric_limits<T>::max();
}

// Constructor:  someGlobalVariable == specificValue
template <class T> selection<T>::selection( T* myCutVar, T eqval ) {
  cutVar = myCutVar;
  cms3Function = NULL;
  equal_val = eqval;
  minimum = std::numeric_limits<T>::max();
  maximum = std::numeric_limits<T>::lowest();
}

// Generic template function to evaluate if the current event passes this selection
template <class T> bool selection<T>::Pass() {
  if( cms3Function != NULL ) {
	if( cms3Function() > minimum &&
		cms3Function() <= maximum ) return true;
	if( cms3Function() == equal_val ) return true;
	return false;
  }
  else if( cutVar != NULL ) {
	if( *cutVar > minimum &&
		*cutVar <= maximum ) return true;
	if( *cutVar == equal_val ) return true;
	return false;
  }
  else {
	std::cout << "Error in sigRegion.cc: Function pointer and variable pointer are both null!" << std::endl;
	throw(5);
  }
  return false;
}

///////////////////////////////////////////////////////////////
// Specialized constructor and "pass" functions for <bool> type

template<> selection<bool>::selection( const bool& (*func)(), bool eqval ) {
  cms3Function = func;
  cutVar = NULL;
  equal_val = eqval;
}

template<> selection<bool>::selection( bool* myCutVar, bool eqval ) {
  cms3Function = NULL;
  cutVar = myCutVar;
  equal_val = eqval;
}

template<> bool selection<bool>::Pass() {
  if( cms3Function != NULL ) {
	if( cms3Function() == equal_val ) return true;
	return false;
  }
  else if( cutVar != NULL ) {
	if( *cutVar == equal_val ) return true;
	return false;
  }
  else {
	std::cout << "Error in sigRegion.cc: Function pointer and variable pointer are both null!" << std::endl;
	throw(5);
  }
  return false;
}

////////////////////////////////////////////////////
// Specialized "pass" function for the <int> type
// Notice that the lower bound is now inclusive

template <> bool selection<int>::Pass() {
  if( cms3Function != NULL ) {
	if( cms3Function() >= minimum &&
		cms3Function() <= maximum ) return true;
	if( cms3Function() == equal_val ) return true;
	return false;
  }
  else if( cutVar != NULL ) {
	if( *cutVar >= minimum &&
		*cutVar <= maximum ) return true;
	if( *cutVar == equal_val ) return true;
	return false;
  }
  else {
	std::cout << "Error in sigRegion.cc: Function pointer and variable pointer are both null!" << std::endl;
	throw(5);
  }
  return false;
}


//---------------------------------------------------------------------------------------------------------------//


///////////////////////////////////////////////
// Function definitions for "sigRegion" class

// Constructors
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

// Other functions
void sigRegion::addSelection( selectionBase* mySelection ) { selections.push_back(mySelection); }

void sigRegion::addSelections( int count, ... ) {
  va_list args;
  va_start( args, count );
  for( int i=0; i<count; i++ ) {
	selectionBase* mySel = va_arg( args, selectionBase* );
	selections.push_back( mySel );
  }
  va_end( args );
}

TString sigRegion::GetLabel()     { return static_cast<TString>(label);      }
TString sigRegion::GetTableName() { return static_cast<TString>(name_table); }
TString sigRegion::GetRootName()  { return static_cast<TString>(name_root);  }

bool sigRegion::PassAllCuts() {
  for( selectionBase* thisSelection : selections ) if( !thisSelection->Pass() ) return false;
  return true;
}


