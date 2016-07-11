#ifndef SIGREGION_H
#define SIGREGION_H

#include <limits>

#include "TString.h"



// Selection class: a simple class to hold a pointer to a CMS3 function, and min/max/equal cut values


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
  selection();
  bool Pass();

protected:
  const T& (*cms3Function)();
  T minimum;
  T maximum;
  T equal_val;

};




// Signal region class: Essentially a container of cuts, with identifiers, and a function to evaluate the logical AND of all cuts

class sigRegion{

public:
  sigRegion( std::string myLabel, std::string niceName );
  sigRegion( std::string myLabel, std::string tabName, std::string rootName );
  void addSelection( selectionBase* mySelection );
  void addSelections( int count, ... );
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
