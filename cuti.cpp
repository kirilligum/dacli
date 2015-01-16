#include <regex>
#include <unordered_set>
#include "headersnb.hpp"

using namespace std;

template <typename T>
auto parse_header(T hin){
  int empty_begin = 0;
  int empty_end = 0;
  vector<unsigned> empty = {0,0};
  vector<regex> names_to_match,beginning_of_range_to_match,ending_of_range_to_match;
  istringstream colsinp(hin);
  for(string t; getline(colsinp,t,',');){
    if(find(begin(t),end(t),'-')!=end(t)){
      istringstream ist(t);
      string range_begin,range_end;
      getline(ist,range_begin,'-');
      if(range_begin.empty()) empty_begin =1;
      beginning_of_range_to_match.emplace_back(range_begin);
      getline(ist,range_end,'-');
      if(range_end.empty()) empty_end =1;
      ending_of_range_to_match.emplace_back(range_end);
      names_to_match.emplace_back(range_end);///> puts the end range into a list of single columns to match. after the ranged is closed the columns that match the closing of the range still match.could be a problem if the columns matching the end of the range come before the begining #BUG
      string check_if_more_than_two;
      if(getline(ist,check_if_more_than_two,'-'))
        cout << "error: correct range syntax should be: beginning-end. you have another -\n";
      //cout <<range_begin<< "  " << range_end<< "  " << check_if_more_than_two<< ".\n";
    }else{
      names_to_match.emplace_back(t);
    }
  }
  empty[0]=empty_begin;
  empty[1]=empty_end;
  return make_tuple(empty,names_to_match,beginning_of_range_to_match,ending_of_range_to_match);
  //return make_tuple(empty_begin,empty_end,names_to_match,beginning_of_range_to_match,ending_of_range_to_match);
}

template<typename E, typename B, typename Counter,typename Op>
void compare_with_range(E element,B& bounds,Counter& in_range,Op op){
  auto comp_rb = find_if(begin(bounds),end(bounds),[=](regex b){return regex_match(element,b);});
  while(comp_rb!=end(bounds)){
    bounds.erase(comp_rb);
    in_range = op(in_range,1); //++in_range or --
    comp_rb = find_if(begin(bounds),end(bounds),[=](regex b){return regex_match(element,b);});
  }
}

template <typename T>
auto select_cols(T header_line){
  vector<regex> mcols,mcols_rbegins,mcols_rends;////note set,unordered_multiset don't work with regex due to absense of comparision or hashing
  vector<unsigned> mcols_empty;
  tie(mcols_empty,mcols,mcols_rbegins,mcols_rends) = parse_header(header_line);
  vector<string> cols;
  vector<size_t> colsi;
  string colss;
  getline(cin,colss);
  istringstream colsiss(colss);
  size_t icol=1;
  unsigned in_col_range=0;
  if(mcols_empty[0]) in_col_range=1;
  string tmp_c;getline(colsiss,tmp_c,',');
  for(string c;getline(colsiss,c,',');) {
    compare_with_range(c,mcols_rbegins,in_col_range,plus<unsigned>());
    if(any_of(begin(mcols),end(mcols),[=](regex r){return regex_match(c,r);})||in_col_range){
      cols.push_back(c);
      colsi.push_back(icol);
    }
    ++icol;
    if(in_col_range) compare_with_range(c,mcols_rends,in_col_range,minus<unsigned>());
  }
  return make_tuple(cols,colsi);
}

template <typename CS,typename T>
auto select_rows_get_data(CS colsi, T header_line){
  vector<vector<string>> data;
  vector<regex> mrows,mrows_rbegins,mrows_rends;
  vector<unsigned> mrows_empty;
  tie(mrows_empty,mrows,mrows_rbegins,mrows_rends) = parse_header(header_line);
  unsigned in_row_range=0;
  if(mrows_empty[0]) in_row_range=1;
  for(string line; getline(cin,line);){
    vector<string> cells;
    istringstream issl(line);
    string first_cell;
    getline(issl,first_cell,',');
    cells.push_back(first_cell);
    compare_with_range(first_cell,mrows_rbegins,in_row_range,plus<unsigned>());
    if(any_of(begin(mrows),end(mrows),[=](regex r){return regex_match(first_cell,r);})||in_row_range){
      size_t icol=1;
      for(string cell; getline(issl,cell,',');){
        if(find(begin(colsi),end(colsi),icol)!=end(colsi)){
          cells.push_back(cell);
        }
        ++icol;
      }
      data.push_back(cells);
    }
    if(in_row_range) compare_with_range(first_cell,mrows_rends,in_row_range,minus<unsigned>());
  }
  return data;
}

int main(int argc, char* argv[]){
  cin.sync_with_stdio(false);///> makes reading stdin a lot faster

  vector<string> cols;
  vector<size_t> colsi;
  tie(cols,colsi) = select_cols(argv[2]);///> select cols listed in argv[2] from the first line of stdin


  vector<vector<string>> data;
  data = select_rows_get_data(colsi,argv[1]);///> checks if the row is listed in argv[1] and reads its data while acounting for selected cols

  ///
  /// output of the dataframe
  ///
  cout << "id,";
  for(size_t i=0; i<cols.size(); ++i) {
    cout << cols[i];
    if(i!=cols.size()-1) cout << ",";
  }
  cout << endl;
  for(size_t i=0; i<data.size();++i){
    cout << data[i][0]<<",";
    for(size_t j=1; j<data[i].size(); ++j) {
      cout << data[i][j];
      if(j!=data[i].size()-1)
        cout << ",";
    }
    cout << endl;
  }

  return 0;
}
