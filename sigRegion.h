#ifndef SIGREGION_H
#define SIGREGION_H

#include <iostream>
#include <limits>

#include "TString.h"



// Selection class:
// A simple class to hold a pointer to a CMS3 function or global variable, and min/max/equal cut values


// Pure virtual base class
// (allows us to create a vector of generic selections)
class selectionBase{
public:
  virtual bool Pass()=0;
};

// Derived template classes
template <class T> class selection : public selectionBase {

public:
  selection( const T& (*func)(), T minval, T maxval);
  selection( const T& (*func)(), T eqval );
  selection( T* myCutVar, T minval, T maxval);
  selection( T* myCutVar, T eqval );
  bool Pass();

protected:
  const T& (*cms3Function)();
  T* cutVar;
  T minimum;
  T maximum;
  T equal_val;

};

//---------------------------------------------------------------------------------------------------------------//


// Signal region class:
// Essentially a name, a vector of selections, and a function that checks if we pass all those selections

class sigRegion{

public:
  sigRegion( std::string myLabel, std::string niceName );
  sigRegion( std::string myLabel, std::string tabName, std::string rootName );
  void AddSelection( selectionBase* mySelection );
  void AddSelections( std::vector<selectionBase*> mySelections );
  TString GetLabel();
  TString GetTableName();
  TString GetRootName();
  bool PassAllCuts();


private:
  std::string label;
  std::string name_table;
  std::string name_root;
  std::vector<selectionBase*> selections;

};


#endif
