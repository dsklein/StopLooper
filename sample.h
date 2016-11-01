#ifndef SAMPLE_H
#define SAMPLE_H


#include <string>

#include "TString.h"
#include "TChain.h"
#include "sigRegion.h"


class sample{
public:
	enum sampleType { kData, kSignal, kBackground };

	sample(std::string myLabel, std::string niceName);
	sample(std::string myLabel, std::string niceName, short int color, sampleType type);
	sample(std::string myLabel, std::string tabName, std::string legName);
	sample(std::string myLabel, std::string tabName, std::string legName, short int color, sampleType type);

	void    AddFile(TString filename);
	void    AddSelections( std::vector<selectionBase*> newselections );
	TString GetLabel();
	TString GetTableName();
	TString GetLegName();
	bool    IsData();
	bool    IsSignal();
	bool    IsBkg();
	short int GetColor();
	TChain* GetChain();
	bool PassSelections();

	void SetNiceName(std::string name);
	void SetColor(short int color);
	void SetSampleType(sampleType type);

private:
	std::string storage_label;
	std::string name_table;
	std::string name_legend;
	sampleType myType;
	short int hist_color;
	TChain chain;
	std::vector<selectionBase*> selections;

};

#endif
