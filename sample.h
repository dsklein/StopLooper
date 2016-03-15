#ifndef SAMPLE_H
#define SAMPLE_H


#include <string>

#include "TString.h"


class sample{
public:
  sample(std::string stName, std::string niceName);
  sample(std::string stName, std::string niceName, short int color, bool isdata, bool issignal);
  sample(std::string stName, std::string tabName, std::string legName);
  sample(std::string stName, std::string tabName, std::string legName, short int color, bool isdata, bool issignal);

  TString GetIntName();
  TString GetTableName();
  TString GetLegName();
  bool    IsData();
  bool    IsSignal();
  short int GetColor();

  void SetNiceName(std::string name);
  void SetColor(short int color);
  void SetDataFlag(bool isdata);
  void SetSignalFlag(bool issignal);

private:
  std::string name_storage;
  std::string name_table;
  std::string name_legend;
  bool is_Data;
  bool is_Signal;
  short int hist_color;
  // TChain ch;

};

#endif
