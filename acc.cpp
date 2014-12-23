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
#include <boost/accumulators/statistics/min.hpp>
#include <boost/accumulators/statistics/extended_p_square.hpp>
#include <boost/accumulators/statistics/max.hpp>

using namespace std;
using namespace boost::accumulators;

int main()
//int main(int argc, char *argv[])
{

  string header; getline(cin,header);
  vector<string> col_names;
  istringstream iss_header(header);
  string tmp_tok; getline(iss_header,tmp_tok,',');
  for(string tok; getline(iss_header,tok,',');col_names.push_back(tok));

  //// by default, display (based on pandas describe) count, mean, std, min, 25, 50, 75th precentile, and max for each column excluding NaN. output number of nans and missing values.
  vector<double> probs {0.25,0.5,0.75};
  typedef accumulator_set<double, features<
    tag::count
    , tag::mean
    , tag::variance
    , tag::min
    , tag::extended_p_square
    , tag::max
    >> acc_type;
  vector<acc_type> accs(col_names.size(), acc_type( tag::extended_p_square::probabilities = probs));
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

  cout << "acc_name,"; for(size_t i=0; i<col_names.size();++i){
    cout << col_names[i];
    if(i!=col_names.size()-1) cout << ',';
  } cout << endl;
  cout << "count,"; for(size_t i=0; i<accs.size();++i){
    cout << extract_result<tag::count>(accs[i]);
    if(i!=col_names.size()-1) cout << ',';
  } cout << endl;
  cout << "mean,"; for(size_t i=0; i<accs.size();++i){
    cout << extract_result<tag::mean>(accs[i]);
    if(i!=col_names.size()-1) cout << ',';
  } cout << endl;
  cout << "stdev,"; for(size_t i=0; i<accs.size();++i){
    cout << sqrt(variance(accs[i]));
    if(i!=col_names.size()-1) cout << ',';
  } cout << endl;
  cout << "min,"; for(size_t i=0; i<accs.size();++i){
    cout << extract_result<tag::min>(accs[i]);
    if(i!=col_names.size()-1) cout << ',';
  } cout << endl;
  cout << "25%,"; for(size_t i=0; i<accs.size();++i){
    cout << extract_result<tag::extended_p_square>(accs[i])[0];
    if(i!=col_names.size()-1) cout << ',';
  } cout << endl;
  cout << "50%,"; for(size_t i=0; i<accs.size();++i){
    cout << extract_result<tag::extended_p_square>(accs[i])[1];
    if(i!=col_names.size()-1) cout << ',';
  } cout << endl;
  cout << "75%,"; for(size_t i=0; i<accs.size();++i){
    cout << extract_result<tag::extended_p_square>(accs[i])[2];
    if(i!=col_names.size()-1) cout << ',';
  } cout << endl;
  cout << "max,"; for(size_t i=0; i<accs.size();++i){
    cout << extract_result<tag::max>(accs[i]);
    if(i!=col_names.size()-1) cout << ',';
  } cout << endl;

  return 0;
}

  //std::cout << "Moment: " << accumulators::moment<2>(acc) << std::endl;
