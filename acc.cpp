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
#include <boost/accumulators/statistics/sum_kahan.hpp>
#include <boost/accumulators/statistics/moment.hpp>
#include <boost/accumulators/statistics/variance.hpp>
#include <boost/accumulators/statistics/skewness.hpp>
#include <boost/accumulators/statistics/kurtosis.hpp>
#include <boost/accumulators/statistics/min.hpp>
#include <boost/accumulators/statistics/extended_p_square.hpp>
#include <boost/accumulators/statistics/p_square_quantile.hpp>
#include <boost/accumulators/statistics/max.hpp>
#include <boost/accumulators/statistics/density.hpp>

using namespace std;
using namespace boost::accumulators;

int main()
//int main(int argc, char *argv[])
{
  cin.sync_with_stdio(false);

  string header; getline(cin,header);
  vector<string> col_names;
  istringstream iss_header(header);
  string tmp_tok; getline(iss_header,tmp_tok,',');
  for(string tok; getline(iss_header,tok,',');col_names.push_back(tok));

  //// by default, display (based on pandas describe) count, mean, std, min, 25, 50, 75th precentile, and max for each column excluding NaN. output number of nans and missing values.
  vector<double> probs {0.5};
  //vector<double> probs {0.25,0.5,0.75};
  typedef accumulator_set<double, features<
    tag::count
    //, tag::mean
    //, tag::moment<2>
    //, tag::moment<3>
    //, tag::moment<4>
    //, tag::sum_kahan
    //, tag::variance
    //, tag::skewness
    //, tag::kurtosis
    //, tag::min
    //, tag::extended_p_square
    , tag::p_square_quantile
    //, tag::max
    //, tag::density
    >> acc_type;
  vector<acc_type> accs(col_names.size(), acc_type( quantile_probability=0.5, tag::extended_p_square::probabilities = probs, tag::density::num_bins=8,tag::density::cache_size=40));
  vector<size_t> missing(col_names.size(),0);
  vector<size_t> infs(col_names.size(),0);
  vector<size_t> nans(col_names.size(),0);
  for(string line;getline(cin,line);){
    istringstream issl(line);
    string idx; getline(issl,idx,',');
    size_t icol =0;
    for(auto& i :missing) i++;
    for(string tok;getline(issl,tok,','); ){
      double x;
      if(!tok.empty()){
        missing[icol]--;
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
  //cout << "missing,"; for(size_t i=0; i<missing.size();++i){
    //cout << missing[i];
    //if(i!=col_names.size()-1) cout << ',';
  //} cout << endl;
  //cout << "infs,"; for(size_t i=0; i<missing.size();++i){
    //cout << infs[i];
    //if(i!=col_names.size()-1) cout << ',';
  //} cout << endl;
  //cout << "nans,"; for(size_t i=0; i<missing.size();++i){
    //cout << nans[i];
    //if(i!=col_names.size()-1) cout << ',';
  //} cout << endl;
  //int show_h=0;
  cout << "count,"; for(size_t i=0; i<accs.size();++i){
    int tmp_count = extract_result<tag::count>(accs[i]);
    //if(tmp_count>40) ++show_h;
    cout << tmp_count;
    //cout << extract_result<tag::count>(accs[i]);
    if(i!=col_names.size()-1) cout << ',';
  } cout << endl;
  //cout << "mean,"; for(size_t i=0; i<accs.size();++i){
    //cout << extract_result<tag::mean>(accs[i]);
    //if(i!=col_names.size()-1) cout << ',';
  //} cout << endl;
  //cout << "moment<2>,"; for(size_t i=0; i<accs.size();++i){
    //cout << extract_result<tag::moment<2>>(accs[i]);
    //if(i!=col_names.size()-1) cout << ',';
  //} cout << endl;
  //cout << "moment<3>,"; for(size_t i=0; i<accs.size();++i){
    //cout << extract_result<tag::moment<3>>(accs[i]);
    //if(i!=col_names.size()-1) cout << ',';
  //} cout << endl;
  //cout << "sum,"; for(size_t i=0; i<accs.size();++i){
    //cout << extract_result<tag::sum_kahan>(accs[i]);
    //if(i!=col_names.size()-1) cout << ',';
  //} cout << endl;
  //cout << "stdev,"; for(size_t i=0; i<accs.size();++i){
    //cout << sqrt(variance(accs[i]));
    //if(i!=col_names.size()-1) cout << ',';
  //} cout << endl;
  //cout << "variance,"; for(size_t i=0; i<accs.size();++i){
    //cout << variance(accs[i]);
    //if(i!=col_names.size()-1) cout << ',';
  //} cout << endl;
  //cout << "skewness,"; for(size_t i=0; i<accs.size();++i){
    //cout << skewness(accs[i]);
    //if(i!=col_names.size()-1) cout << ',';
  //} cout << endl;
  //cout << "kurtosis,"; for(size_t i=0; i<accs.size();++i){
    //cout << kurtosis(accs[i]);
    //if(i!=col_names.size()-1) cout << ',';
  //} cout << endl;
  //cout << "min,"; for(size_t i=0; i<accs.size();++i){
    //cout << extract_result<tag::min>(accs[i]);
    //if(i!=col_names.size()-1) cout << ',';
  //} cout << endl;
  //cout << "25%,"; for(size_t i=0; i<accs.size();++i){
    //cout << extract_result<tag::p_square_quantile>(accs[i])[0];
    //if(i!=col_names.size()-1) cout << ',';
  //} cout << endl;
  cout << "50%,"; for(size_t i=0; i<accs.size();++i){
    cout << extract_result<tag::p_square_quantile>(accs[i]);
    if(i!=col_names.size()-1) cout << ',';
  } cout << endl;
  //cout << "50%,"; for(size_t i=0; i<accs.size();++i){
    //cout << extract_result<tag::extended_p_square>(accs[i])[0];
    //if(i!=col_names.size()-1) cout << ',';
  //} cout << endl;
  //cout << "75%,"; for(size_t i=0; i<accs.size();++i){
    //cout << extract_result<tag::p_square_quantile>(accs[i])[2];
    //if(i!=col_names.size()-1) cout << ',';
  //} cout << endl;
  //cout << "max,"; for(size_t i=0; i<accs.size();++i){
    //cout << extract_result<tag::max>(accs[i]);
    //if(i!=col_names.size()-1) cout << ',';
  //} cout << endl;
  //if(show_h){
    //vector<vector<double>> hm;
    //vector<double> h_max;
    //for(auto i: accs) {
      //vector<double> hv;
      //for(int j=0; j<extract_result<tag::density>(i).size();++j){
        ////cout << extract_result<tag::density>(i)[j].second << " ";
        //hv.push_back(extract_result<tag::density>(i)[j].second);
      //}
      ////h_max.push_back(10.0);
      //h_max.push_back(*std::max_element(hv.begin(),hv.end()));
      //hm.push_back(hv);
    //}
    //for(int ih=0; ih <10; ++ih){
      //cout << "hist_"<<ih<<","; for(size_t i=0; i<accs.size();++i){
        //double level = hm[i][ih]/h_max[i]*8;
        ////cout << level;
        //for(unsigned i=0; i<level;++i) cout << '#';
        //if(i!=col_names.size()-1) cout << ',';
      //} cout << endl;
    //}
  //}
  return 0;
}

  //std::cout << "Moment: " << accumulators::moment<2>(acc) << std::endl;
