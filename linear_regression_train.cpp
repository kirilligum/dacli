#include <iostream>
#include <iterator>
#include <vector>
#include <deque>
#include <sstream>
#include <string>
#include <cmath>
#include <cfloat>
#include <stdexcept>
#include <algorithm>
#include <numeric>
#include <unistd.h>
#include "read_header.hpp"
#include "counter.hpp"

using namespace std;

auto h_sample(const vector<double>& theta, const vector<double>& x){
  double h=theta[0];
  for(size_t i=0, n=x.size();i<n;++i){
    h+=theta[i+1]*x[i];
  }
  return 0;
}

int main() {
  vector<vector<double>> data;
  vector<string> rows;
  auto col_names = read_header(cin);
  typedef counter<unsigned __int128>  count_type;
  vector<count_type> missing(col_names.size(),count_type());
  vector<count_type> infs(col_names.size(),count_type());
  vector<count_type> nans(col_names.size(),count_type());
  vector<count_type> count(col_names.size(),count_type());
  for(string line;getline(cin,line);){
    istringstream issl(line);
    string idx; getline(issl,idx,',');///> skip the first field since it's id
    size_t icol =0;
    vector<double> line_values;
    for(string tok;getline(issl,tok,','); ){
      double x;
      if(!tok.empty()){
        try {
          x = stod(tok);
        } catch (const std::invalid_argument&) {
          cerr << "Error: stod convertion of the cell is invalid.\n";
          throw;
        } catch (const std::out_of_range&) {
          cerr << "Error: stod convertion of the cell is out of range.\n";
          throw;
        }
        if(std::isfinite(x)){
          line_values.push_back(x);
        }
        else if(std::isinf(x))
          infs[icol].increment();
        else if(std::isnan(x))
          nans[icol].increment();
        else
          cerr << "error: caticorization of input cell didn't work.\n";
      }else {
        missing[icol].increment();
      }
      ++icol;
    }
    data.push_back(line_values);
  }
  cout << data.size() << endl;
  for(auto l: data) {
    for(auto i:l) cout << i << " " ;
    cout << endl;
  }
  vector<double> theta(data.front().size()+1,0.0);
  vector<double> h;
  for(auto i: data) h.push_back(h_sample(theta,i));
  for(auto i : h) cout << i << " "; cout << endl;
  double cost=0.0;
  return 0;
}
