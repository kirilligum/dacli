//#include "headers.hpp"
//#include "strtk.hpp"
#include <iostream>
#include <vector>
#include <sstream>
#include <string>
#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics/stats.hpp>
#include <boost/accumulators/statistics/mean.hpp>
#include <boost/accumulators/statistics/moment.hpp>

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

  vector<accumulator_set<double, stats<tag::mean, tag::moment<2> > >> accs(col_names.size());
  for(string line;getline(cin,line);){
    istringstream issl(line);
    string idx; getline(issl,idx,',');
    size_t icol =0;
    for(string tok;getline(issl,tok,','); ){
      accs[icol++](stod(tok));
    }
  }

  for(size_t i=0; i<col_names.size();++i){
    cout << col_names[i];
    if(i!=col_names.size()-1) cout << ',';
    else cout << " # acc_name";
  }
  cout << endl;
  for(size_t i=0; i<accs.size();++i){
    cout << mean(accs[i]);
    if(i!=col_names.size()-1) cout << ',';
    else cout << " # mean";
  }
  cout << endl;

  return 0;
}

  //std::cout << "Moment: " << accumulators::moment<2>(acc) << std::endl;
