#pragma once
#include <vector>
#include <sstream>
#include <string>

using namespace std;

template<typename T>
auto read_header(T& in) {
  string header; getline(in,header);
  vector<string> col_names;
  istringstream iss_header(header);
  string tmp_tok; getline(iss_header,tmp_tok,',');
  for(string tok; getline(iss_header,tok,',');col_names.push_back(tok));
  return col_names;
}
