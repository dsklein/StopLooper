#include "sample.h"


// Constructors

sample::sample(std::string stName, std::string niceName)
  : name_storage(stName),
	name_table(niceName),
	name_legend(niceName),
	myType(kBackground),
	hist_color(kBlack)
{
}

sample::sample(std::string stName, std::string niceName, short int color, sampleType type = kBackground)
  : name_storage(stName),
	name_table(niceName),
	name_legend(niceName),
	myType(type),
	hist_color(color)
{
}

sample::sample(std::string stName, std::string tabName, std::string legName)
  : name_storage(stName),
	name_table(tabName),
	name_legend(legName),
	myType(kBackground),
	hist_color(kBlack)
{
}

sample::sample(std::string stName, std::string tabName, std::string legName, short int color, sampleType type = kBackground)
  : name_storage(stName),
	name_table(tabName),
	name_legend(legName),
	myType(type),
	hist_color(color)
{
}

// Other member functions
// Get or set various properties

TString sample::GetIntName()   { return static_cast<TString>(name_storage); }
TString sample::GetTableName() { return static_cast<TString>(name_table);   }
TString sample::GetLegName()   { return static_cast<TString>(name_legend);  }
bool    sample::IsData()       { return (myType==kData);       }
bool    sample::IsSignal()     { return (myType==kSignal);     }
bool    sample::IsBkg()        { return (myType==kBackground); }
short int sample::GetColor()   { return hist_color;   }

void sample::SetNiceName(std::string name) { name_table=name; name_legend=name; }
void sample::SetColor(short int color)     { hist_color=color; }
void sample::SetSampleType(sampleType type = kBackground)  { myType=type; }
