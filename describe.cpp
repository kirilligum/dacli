//#include "headers.hpp"
//#include "strtk.hpp"
#include <iostream>
#include <vector>
#include <sstream>
#include <string>
#include <cmath>
#include <cfloat>
#include <stdexcept>

using namespace std;

template<typename T>
auto read_header(T& in) {
  string header; getline(cin,header);
  vector<string> col_names;
  istringstream iss_header(header);
  string tmp_tok; getline(iss_header,tmp_tok,',');
  for(string tok; getline(iss_header,tok,',');col_names.push_back(tok));
  return col_names;
}

int main()
//int main(int argc, char *argv[])
{
  cin.sync_with_stdio(false);
  auto col_names = read_header(cin);

  //// by default, display (based on pandas describe) count, mean, std, min, 25, 50, 75th precentile, and max for each column excluding NaN. output number of nans and missing values.
  vector<size_t> missing(col_names.size(),0);
  vector<size_t> infs(col_names.size(),0);
  vector<size_t> nans(col_names.size(),0);
  vector<size_t> count(col_names.size(),0);
  vector<double> mean(col_names.size(),0);
  vector<double> m2(col_names.size(),0);
  vector<vector<double>> xv(col_names.size());
  vector<double> var2p(col_names.size(),0);
  vector<double> sum(col_names.size(),0);
  vector<int> min(col_names.size(),0);
  vector<int> max(col_names.size(),0);
  vector<size_t> q(col_names.size(),0);
  vector<size_t> q_old(col_names.size(),0);
  for(string line;getline(cin,line);){
    istringstream issl(line);
    string idx; getline(issl,idx,',');///> skip the first field since it's id
    size_t icol =0;
    //for(auto& i :missing) i++;
    for(string tok;getline(issl,tok,','); ){
      //count[icol]++;
      double x;
      if(!tok.empty()){
        //missing[icol]--;
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
          count[icol]++;
          //sum[icol]+=x;
          double delta = x - mean[icol];
          mean[icol]=mean[icol]+delta/count[icol];
          m2[icol] += delta*(x-mean[icol]);
          //xv[icol].push_back(x);
        }
        else if(std::isinf(x))
          infs[icol]++;
        else if(std::isnan(x))
          nans[icol]++;
        else
          cerr << "error: caticorization of input cell didn't work.\n";
      }else {
        missing[icol]++;
      }
      ++icol;
    }
    ///
    /// second pass of two-pass . turned out to be too slow
    ///
    //size_t xvn = xv[0].size();
    //for(size_t icol=0; icol<col_names.size();++icol){
      //double m= sum[icol]/count[icol];
      //mean[icol]= m;
      //double sum2=0,sum3=0;
      //for(size_t i=0;i<xvn;++i){
        //sum2+=pow((xv[icol][i]-m),2);
        //sum3+=(xv[icol][i]-m);
      //}
      //var2p[icol]=(sum2-pow(sum3,2)/count[icol])/(count[icol]-1);
    //}
  }

  cout << "acc_name,"; for(size_t i=0; i<col_names.size();++i){
    cout << col_names[i];
    if(i!=col_names.size()-1) cout << ',';
  } cout << endl;
  cout << "missing,"; for(size_t i=0; i<missing.size();++i){
    cout << missing[i];
    if(i!=col_names.size()-1) cout << ',';
  } cout << endl;
  cout << "infs,"; for(size_t i=0; i<missing.size();++i){
    cout << infs[i];
    if(i!=col_names.size()-1) cout << ',';
  } cout << endl;
  cout << "nans,"; for(size_t i=0; i<missing.size();++i){
    cout << nans[i];
    if(i!=col_names.size()-1) cout << ',';
  } cout << endl;
  cout << "count,"; for(size_t i=0; i<missing.size();++i){
    cout << count[i];
    if(i!=col_names.size()-1) cout << ',';
  } cout << endl;
  cout << "mean,"; for(size_t i=0; i<missing.size();++i){
    cout << mean[i];
    if(i!=col_names.size()-1) cout << ',';
  } cout << endl;
  cout << "var,"; for(size_t i=0; i<missing.size();++i){
    //cout << var2p[i];
    cout << m2[i]/count[i];
    if(i!=col_names.size()-1) cout << ',';
  } cout << endl;

  //cout << "mean,"; for(size_t i=0; i<accs.size();++i){
    //cout << extract_result<tag::mean>(accs[i]);
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
    //cout << extract_result<tag::extended_p_square>(accs[i])[0];
    //if(i!=col_names.size()-1) cout << ',';
  //} cout << endl;
  //cout << "50%,"; for(size_t i=0; i<accs.size();++i){
    //cout << extract_result<tag::extended_p_square>(accs[i])[1];
    //if(i!=col_names.size()-1) cout << ',';
  //} cout << endl;
  //cout << "75%,"; for(size_t i=0; i<accs.size();++i){
    //cout << extract_result<tag::extended_p_square>(accs[i])[2];
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

