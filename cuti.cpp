#include <regex>
#include <unordered_set>
#include "headers.hpp"

using namespace std;

template <typename T>
auto parse_header_in(T hin){
  vector<regex> names_to_match,beginning_of_range_to_match,ending_of_range_to_match;
  istringstream colsinp(hin);
  for(string t; getline(colsinp,t,',');){
    if(find(begin(t),end(t),'-')!=end(t)){
      istringstream ist(t);
      string range_begin,range_end;
      getline(ist,range_begin,'-');
      beginning_of_range_to_match.emplace_back(range_begin);
      getline(ist,range_end,'-');
      ending_of_range_to_match.emplace_back(range_end);
      string check_if_more_than_two;
      if(getline(ist,check_if_more_than_two,'-'))
        cout << "error: correct range syntax should be: beginning-end. there is another -\n";
      //cout <<range_begin<< "  " << range_end<< "  " << check_if_more_than_two<< ".\n";
    }else{
      names_to_match.emplace_back(t);
    }
  }
  return make_tuple(names_to_match,beginning_of_range_to_match,ending_of_range_to_match);
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

int main(int argc, char* argv[]){
  vector<regex> mrows,mrows_rbegins,mrows_rends;////note set,unordered_multiset don't work with regex due to absense of comparision or hashing
  tie(mrows,mrows_rbegins,mrows_rends) = parse_header_in(argv[1]);
  unsigned in_row_range=0;
  vector<string> rows;
  vector<size_t> rowsi;

  ///// parese column header
  vector<regex> mcols,mcols_rbegins,mcols_rends;
  tie(mcols,mcols_rbegins,mcols_rends) = parse_header_in(argv[2]);
  vector<string> cols;
  vector<size_t> colsi;
  string colss;
  getline(cin,colss);
  istringstream colsiss(colss);
  size_t icol=1;
  unsigned in_col_range=0;
  for(string c;getline(colsiss,c,',');) {
    compare_with_range(c,mcols_rbegins,in_col_range,plus<unsigned>());
    if(any_of(begin(mcols),end(mcols),[=](regex r){return regex_match(c,r);})||in_col_range){
      cols.push_back(c);
      colsi.push_back(icol);
    }
    ++icol;
    compare_with_range(c,mcols_rends,in_col_range,minus<unsigned>());
  }

  vector<vector<string>> data;
  size_t irow=1;
  for(string line; getline(cin,line);){
    vector<string> cells;
    istringstream issl(line);
    //// parse  row  header (row's numbers)
    string first_cell;
    getline(issl,first_cell,',');
    cells.push_back(first_cell);
    compare_with_range(first_cell,mrows_rbegins,in_row_range,plus<unsigned>());
    if(any_of(begin(mrows),end(mrows),[=](regex r){return regex_match(first_cell,r);})||in_row_range){
      icol=1;
      if(find(begin(colsi),end(colsi),icol)!=end(colsi)){
        cells.push_back(first_cell);
      }
      ++icol;
      rows.push_back(first_cell);
      rowsi.push_back(irow);
      for(string cell; getline(issl,cell,',');){
        if(find(begin(colsi),end(colsi),icol)!=end(colsi)){
          cells.push_back(cell);
        }
        ++icol;
      }
      data.push_back(cells);
      ++irow;
    }
    compare_with_range(first_cell,mrows_rends,in_row_range,minus<unsigned>());
  }

  cout << "_,";
  for(size_t i=0; i<cols.size(); ++i) {
    cout << cols[i];
    if(i!=cols.size()-1)
      cout << ",";
  }
  cout << endl;
  //for(auto i:data){
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
