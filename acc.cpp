//#include "headers.hpp"
//#include "strtk.hpp"
#include <iostream>
#include <vector>
#include <sstream>
#include <string>
#include <cmath>
#include <cfloat>
#include <stdexcept>
#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics/stats.hpp>
#include <boost/accumulators/statistics/mean.hpp>
#include <boost/accumulators/statistics/moment.hpp>
#include <boost/accumulators/statistics/variance.hpp>

using namespace std;
using namespace boost::accumulators;

int main()
//int main(int argc, char *argv[])
{
  //// by default, display (based on pandas describe) count, mean, std, min, 25, 50, 75th precentile, and max for each column excluding NaN. output number of nans and missing values.

  string header; getline(cin,header);
  vector<string> col_names;
  istringstream iss_header(header);
  string tmp_tok; getline(iss_header,tmp_tok,',');
  for(string tok; getline(iss_header,tok,',');col_names.push_back(tok));

  vector<accumulator_set<double, features<
    tag::count,
    tag::mean,
    tag::variance
    >>> accs(col_names.size());
  vector<size_t> missing(col_names.size(),0);
  vector<size_t> infs(col_names.size(),0);
  vector<size_t> nans(col_names.size(),0);
  for(string line;getline(cin,line);){
    istringstream issl(line);
    string idx; getline(issl,idx,',');
    size_t icol =0;
    for(string tok;getline(issl,tok,','); ){
      double x;
      if(tok.empty())
        missing[icol]++;
      else{
        try {
          x = stod(tok);
        } catch (const std::invalid_argument&) {
          cerr << "Error: stod convertion of the cell is invalid.\n";
          throw;
        } catch (const std::out_of_range&) {
          cerr << "Error: stod convertion of the cell is out of range.\n";
          throw;
        }
        if(std::isfinite(x))
          ref(accs[icol])(x);
        else if(std::isinf(x))
          infs[icol]++;
        else if(std::isnan(x))
          nans[icol]++;
        else
          cerr << "error: caticorization of input cell didn't work.\n";
      }
      ++icol;
    }
  }

  for(size_t i=0; i<col_names.size();++i){
    cout << col_names[i];
    if(i!=col_names.size()-1) cout << ',';
    else cout << " ,# acc_name";
  }
  cout << endl;
  for(size_t i=0; i<accs.size();++i){
    cout << count(accs[i]);
    if(i!=col_names.size()-1) cout << ',';
    else cout << " ,# count";
  }
  cout << endl;
  for(size_t i=0; i<accs.size();++i){
    cout << mean(accs[i]);
    if(i!=col_names.size()-1) cout << ',';
    else cout << " ,# mean";
  }
  cout << endl;
  for(size_t i=0; i<accs.size();++i){
    cout << sqrt(variance(accs[i]));
    if(i!=col_names.size()-1) cout << ',';
    else cout << " ,# stdev";
  }
  cout << endl;

  return 0;
}

  //std::cout << "Moment: " << accumulators::moment<2>(acc) << std::endl;
