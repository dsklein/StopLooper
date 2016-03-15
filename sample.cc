#include "sample.h"


// Constructors

sample::sample(std::string stName, std::string niceName)
  : name_storage(stName),
	name_table(niceName),
	name_legend(niceName),
	is_Data(false),
	is_Signal(false),
	hist_color(kBlack)
{
}

sample::sample(std::string stName, std::string niceName, short int color, bool isdata=false, bool issignal=false)
  : name_storage(stName),
	name_table(niceName),
	name_legend(niceName),
	is_Data(isdata),
	is_Signal(issignal),
	hist_color(color)
{
}

sample::sample(std::string stName, std::string tabName, std::string legName)
  : name_storage(stName),
	name_table(tabName),
	name_legend(legName),
	is_Data(false),
	is_Signal(false),
	hist_color(kBlack)
{
}

sample::sample(std::string stName, std::string tabName, std::string legName, short int color, bool isdata=false, bool issignal=false)
  : name_storage(stName),
	name_table(tabName),
	name_legend(legName),
	is_Data(isdata),
	is_Signal(issignal),
	hist_color(color)
{
}

// Other member functions
// Get or set various properties

TString sample::GetIntName()   { return static_cast<TString>(name_storage); }
TString sample::GetTableName() { return static_cast<TString>(name_table);   }
TString sample::GetLegName()   { return static_cast<TString>(name_legend);  }
bool    sample::IsData()       { return is_Data;      }
bool    sample::IsSignal()     { return is_Signal;    }
short int sample::GetColor()   { return hist_color;   }

void sample::SetNiceName(std::string name) { name_table=name; name_legend=name; }
void sample::SetColor(short int color)     { hist_color=color; }
void sample::SetDataFlag(bool isdata = true)      { is_Data=isdata; }
void sample::SetSignalFlag(bool issignal = true)  { is_Signal=issignal; }
