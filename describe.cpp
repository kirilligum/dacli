//#include "headers.hpp"
//#include "strtk.hpp"
#include <iostream>
#include <vector>
#include <sstream>
#include <string>
#include <cmath>
#include <cfloat>
#include <stdexcept>
#include <algorithm>

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
  vector<double> std(col_names.size(),0);
  vector<double> var(col_names.size(),0);
  vector<double> ske(col_names.size(),0);
  vector<double> kur(col_names.size(),0);
  vector<double> sum(col_names.size(),0);
  vector<double> sum2(col_names.size(),0);
  vector<double> sum3(col_names.size(),0);
  vector<double> sum4(col_names.size(),0);
  vector<double> min(col_names.size());
  vector<double> max(col_names.size());
  vector<size_t> q(col_names.size(),0);
  vector<size_t> q_old(col_names.size(),0);
  bool first_line_not_read =1;
  for(string line;getline(cin,line);){
    istringstream issl(line);
    string idx; getline(issl,idx,',');///> skip the first field since it's id
    size_t icol =0;
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
          if(first_line_not_read){
            min[icol]=x;
            max[icol]=x;
            first_line_not_read=0;
          }else{
            if(x<min[icol]) min[icol]=x;
            if(x>max[icol]) max[icol]=x;
          }
          count[icol]++;
          sum[icol]+=x;
          double x2=x*x;
          sum2[icol]+=x2;
          sum3[icol]+=x2*x;
          sum4[icol]+=x2*x2;
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
  }
  for(size_t i=0; i<missing.size();++i){
    double n = count[i];
    double k = sum[i]/n,k2=k*k,k3=k2*k,k4=k2*k2;
    mean[i] = k;
    double sum1s = sum[i]-n*k;
    double sum2s = sum2[i] -2.0*k*sum[i] +n*k2;
    double sum3s = sum3[i] -3*k*sum2[i] +3*k2*sum[i] -n*k3;
    double sum4s = sum4[i] -4*k*sum3[i] +6*k2*sum2[i] -4*k3*sum[i] +n*k4;
    //double sum1s = sum[i];
    //double sum2s = sum2[i];
    //double sum3s = sum3[i];
    //double sum4s = sum4[i];
    double m1 = sum1s/n, m12=m1*m1, m13=m12*m1,m14=m12*m12;
    double m2 = sum2s/n;
    double m3 = sum3s/n;
    double m4 = sum4s/n;
    double variance = m2-m12;
    //double variance = m2-k2;
    //double sample_variance = variance*n/(n-1);
    var[i] = variance;
    double stdev = sqrt(variance);
    //double sample_stdev = sqrt((sum2s - sum1s*k)/(n-1));
    std[i] = stdev;
    double skewness = (m3 -3*m1*m2 +2*m13)/(variance*stdev);
    //double skewness = (sum3s/n -3*k*sum2s/n +2*k3)/((sum2s/n-k2)*sqrt(sum2s/n-k2));
    //double skewness = (sum3s/n -3*k*variance -k3)/(variance*stdev);
    //double skew_bias = sqrt(pow((n-1)/n,3));
    //double sample_skewness = skewness*skew_bias;
    ske[i] = skewness;
    double kurtosis = (m4 -4*m1*m3 +6*m12*m2 -3*m14)/(variance*variance)-3;
    //double kurtosis = (sum4s/n -4*k*variance*stdev +6*k2*variance -4*k3*stdev +k4)/(variance*variance)-3;
    kur[i] = kurtosis;
  }

  cout << "acc_name,"; for(size_t i=0; i<col_names.size();++i){
    cout << col_names[i];
    if(i!=col_names.size()-1) cout << ',';
  } cout << endl;
  if(none_of(begin(missing),end(missing),[](auto m){return m==0?1:0;})){
    cout << "missing,"; for(size_t i=0; i<missing.size();++i){
      cout << missing[i];
      if(i!=col_names.size()-1) cout << ',';
    } cout << endl;
  }
  if(none_of(begin(infs),end(infs),[](auto m){return m==0?1:0;})){
    cout << "infs,"; for(size_t i=0; i<infs.size();++i){
      cout << infs[i];
      if(i!=col_names.size()-1) cout << ',';
    } cout << endl;
  }
  if(none_of(begin(nans),end(nans),[](auto m){return m==0?1:0;})){
    cout << "nans,"; for(size_t i=0; i<nans.size();++i){
      cout << nans[i];
      if(i!=col_names.size()-1) cout << ',';
    } cout << endl;
  }
  cout << "count,"; for(size_t i=0; i<missing.size();++i){
    cout << count[i];
    if(i!=col_names.size()-1) cout << ',';
  } cout << endl;
  cout << "mean,"; for(size_t i=0; i<missing.size();++i){
    cout << mean[i];
    if(i!=col_names.size()-1) cout << ',';
  } cout << endl;
  cout << "m2,"; for(size_t i=0; i<missing.size();++i){
    cout << sum2[i]/count[i];
    if(i!=col_names.size()-1) cout << ',';
  } cout << endl;
  cout << "m3,"; for(size_t i=0; i<missing.size();++i){
    cout << sum3[i]/count[i];
    if(i!=col_names.size()-1) cout << ',';
  } cout << endl;
  cout << "std,"; for(size_t i=0; i<missing.size();++i){
    cout << std[i];
    if(i!=col_names.size()-1) cout << ',';
  } cout << endl;
  cout << "var,"; for(size_t i=0; i<missing.size();++i){
    cout << var[i];
    if(i!=col_names.size()-1) cout << ',';
  } cout << endl;
  cout << "ske,"; for(size_t i=0; i<missing.size();++i){
    cout << ske[i];
    if(i!=col_names.size()-1) cout << ',';
  } cout << endl;
  cout << "kur,"; for(size_t i=0; i<missing.size();++i){
    cout << kur[i];
    if(i!=col_names.size()-1) cout << ',';
  } cout << endl;
  cout << "min,"; for(size_t i=0; i<missing.size();++i){
    cout << min[i];
    if(i!=col_names.size()-1) cout << ',';
  } cout << endl;
  cout << "max,"; for(size_t i=0; i<missing.size();++i){
    cout << max[i];
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

