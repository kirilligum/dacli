#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>

using namespace std;

int main(int argc, char* argv[]){
  fstream iss(argv[1]);
  vector<vector<string>> m;
  for(string line; getline(iss,line);){
    vector<string> cells;
    istringstream issl(line);
    for(string cell; getline(issl,cell,',');cells.push_back(cell));
    m.push_back(cells);
  }

  for(size_t c=0;c<m.front().size();++c){
    for(size_t r=0;r<m.size();++r){
      cout << m[r][c] << ",";
    }
    cout << endl;
  }

  return 0;
}
