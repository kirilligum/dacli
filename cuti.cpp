#include <regex>
#include "headers.hpp"

using namespace std;

int main(int argc, char* argv[]){
  vector<regex> rowsin;
  istringstream rowsinp(argv[1]);
  for(string t; getline(rowsinp,t,',');rowsin.emplace_back(t));
  vector<string> rows;
  vector<size_t> rowsi;

  vector<regex> colsin;
  istringstream colsinp(argv[2]);
  for(string t; getline(colsinp,t,',');colsin.emplace_back(t));
  vector<string> cols;
  vector<size_t> colsi;
  string colss;
  getline(cin,colss);
  istringstream colsiss(colss);
  size_t icol=1;
  for(string c;getline(colsiss,c,',');) {
    if(any_of(begin(colsin),end(colsin),[=](regex r){return regex_match(c,r);})){
      cols.push_back(c);
      colsi.push_back(icol);
    }
    ++icol;
  }

  vector<vector<string>> data;
  size_t irow=1;
  for(string line; getline(cin,line);){
    vector<string> cells;
    istringstream issl(line);
    string first_cell;
    getline(issl,first_cell,',');
    if(any_of(begin(rowsin),end(rowsin),[=](regex r){return regex_match(first_cell,r);})){
      icol=1;
      if(find(begin(colsi),end(colsi),icol)!=end(colsi)){
        cells.push_back(first_cell);
      }
      ++icol;
      rows.push_back(first_cell);
      rowsi.push_back(irow);
      for(string cell; getline(issl,cell,',');){
        if(find(begin(colsi),end(colsi),icol)!=end(colsi)){
          cells.push_back(cell);
        }
        ++icol;
      }
      data.push_back(cells);
      ++irow;
    }
  }

  for(size_t i=0; i<cols.size(); ++i) {
    cout << cols[i];
    if(i!=cols.size()-1)
      cout << ",";
  }
  cout << endl;
  for(auto i:data){
    for(size_t j=0; j<i.size(); ++j) {
      cout << i[j];
      if(j!=i.size()-1)
        cout << ",";
    }
    cout << endl;
  }

  return 0;
}
