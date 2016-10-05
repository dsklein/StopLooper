#include "sample.h"


// Constructors

sample::sample(std::string myLabel, std::string niceName)
	: storage_label(myLabel),
	  name_table(niceName),
	  name_legend(niceName),
	  myType(kBackground),
	  hist_color(kBlack),
	  chain("t")
{
}

sample::sample(std::string myLabel, std::string niceName, short int color, sampleType type = kBackground)
	: storage_label(myLabel),
	  name_table(niceName),
	  name_legend(niceName),
	  myType(type),
	  hist_color(color),
	  chain("t")
{
}

sample::sample(std::string myLabel, std::string tabName, std::string legName)
	: storage_label(myLabel),
	  name_table(tabName),
	  name_legend(legName),
	  myType(kBackground),
	  hist_color(kBlack),
	  chain("t")
{
}

sample::sample(std::string myLabel, std::string tabName, std::string legName, short int color, sampleType type = kBackground)
	: storage_label(myLabel),
	  name_table(tabName),
	  name_legend(legName),
	  myType(type),
	  hist_color(color),
	  chain("t")
{
}

// Other member functions
// Get or set various properties

void    sample::AddFile( TString filename ) { chain.Add(filename); }
TString sample::GetLabel()     { return static_cast<TString>(storage_label); }
TString sample::GetTableName() { return static_cast<TString>(name_table);   }
TString sample::GetLegName()   { return static_cast<TString>(name_legend);  }
bool    sample::IsData()       { return (myType==kData);       }
bool    sample::IsSignal()     { return (myType==kSignal);     }
bool    sample::IsBkg()        { return (myType==kBackground); }
short int sample::GetColor()   { return hist_color;   }
TChain* sample::GetChain()     { return &chain; }

void sample::SetNiceName(std::string name) { name_table=name; name_legend=name; }
void sample::SetColor(short int color)     { hist_color=color; }
void sample::SetSampleType(sampleType type = kBackground)  { myType=type; }
