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
  return h;
}

template<typename T = std::istream>
struct read_data_lines {
  vector<vector<double>> data;
  typedef counter<unsigned __int128>  count_type;
  vector<string> col_names;
  vector<count_type> missing;
  vector<count_type> infs;
  vector<count_type> nans;
  vector<count_type> count;
  vector<string> rows;
  read_data_lines(T& in)
    : missing(col_names.size(),count_type()),
    infs(col_names.size(),count_type()),
    nans(col_names.size(),count_type()),
    count(col_names.size(),count_type()) {
      col_names = read_header(cin);
      for(string line;getline(in,line);){
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
    }
};

template< typename T = vector<double> >
auto cost(const T& h, const T& target){
  auto rows = h.size();
  double j=0.0;
  for(size_t i=0; i<rows; ++i)
    j+=pow(h[i]-target[i],2);
  j/=2*rows;
  return j;
}

auto split_target_train(vector<vector<double>>& data){
  vector<double> target;
  for(auto& l: data) {
    target.push_back(l.back());
    l.pop_back();
  }
  return target;
}

int main() {
  read_data_lines<> rdl(cin);
  auto data = std::move(rdl.data);
  auto target = split_target_train(data);
  //for(auto l: data) { for(auto i:l) cout << i << " " ; cout << endl; }
  //std::copy(begin(target),end(target),std::ostream_iterator<double>(cout," "));cout << endl;
  vector<double> theta(data.front().size()+1,0.0);
  //vector<double> h;
  //for(auto i: data) h.push_back(h_sample(theta,i));
  //for(auto i : h) cout << i << " "; cout << endl;
  //cout << " cost " << cost(h,target) << endl;
  double am = 1.0e-2/target.size();
  double lambda1 =1e-4;
  double lambda2 =1e-2;
  double old_cost=0.0;
  double abserr=1e-4;
  for(size_t i=0;i<1e4;++i){
    for(size_t j=0;j<theta.size();++j){///> th_j=th_j-a/m*Sum(hi-yi)*xij
      double dcost=0;
      vector<double> h;
      for(auto i: data) h.push_back(h_sample(theta,i));
      for(size_t irow=0, nrow=h.size();irow<nrow;++irow){
        if(j==0)
          dcost+=(h[irow]-target[irow]);
        else
          dcost+=(h[irow]-target[irow])*data[irow][j-1];
      }
      theta[j]-=am*dcost+abs(lambda1*theta[j])+pow(lambda2*theta[j],2);
    }
    vector<double> h;
    for(auto i: data) h.push_back(h_sample(theta,i));
    //for(auto i : h) cout << i << " "; cout << endl;
    double current_cost = cost(h,target);
    cout << " cost " << current_cost<< endl;
    std::copy(begin(theta),end(theta),ostream_iterator<double>(cout," ")); cout << endl;
    if(abs(old_cost-current_cost)<abserr) break;
    old_cost =current_cost;
  }
  return 0;
}
