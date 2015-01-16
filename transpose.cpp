#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>

using namespace std;

int main(){
  cin.sync_with_stdio(false);
  vector<vector<string>> m;
  for(string line; getline(cin,line);){
    vector<string> cells;
    istringstream issl(line);
    for(string cell; getline(issl,cell,',');cells.push_back(cell));
    m.push_back(cells);
  }

  for(size_t c=0;c<m.front().size();++c){
    for(size_t r=0;r<m.size();++r){
      cout << m[r][c] ;
      if(r<m.size()-1) cout << ",";
    }
    cout << endl;
  }

  return 0;
}
