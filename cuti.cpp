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
  return make_tuple(names_to_match,mcols_rbegins,mcols_rends);
}


int main(int argc, char* argv[]){
  vector<regex> mrows,mrows_rbegins,mrows_rends;////note set,unordered_multiset don't work with regex due to absense of comparision or hashing
  tie(mrows,mrows_rbegins,mrows_rends) = parse_header_in(argv[1]);
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
  unsigned in_range=0;
  for(string c;getline(colsiss,c,',');) {
    auto comp_rb = find_if(begin(mcols_rbegins),end(mcols_rbegins),[=](regex r){return regex_match(c,r);});
    while(comp_rb!=end(mcols_rbegins)){
      mcols_rbegins.erase(comp_rb);
      ++in_range;
      comp_rb = find_if(begin(mcols_rbegins),end(mcols_rbegins),[=](regex r){return regex_match(c,r);});
    }
    if(any_of(begin(mcols),end(mcols),[=](regex r){return regex_match(c,r);})||in_range){
      cols.push_back(c);
      colsi.push_back(icol);
    }
    ++icol;
    auto comp_re = find_if(begin(mcols_rends),end(mcols_rends),[=](regex r){return regex_match(c,r);});
    while(comp_re!=end(mcols_rends)){
      mcols_rends.erase(comp_re);
      --in_range;
      comp_re = find_if(begin(mcols_rends),end(mcols_rends),[=](regex r){return regex_match(c,r);});
    }
  }

  vector<vector<string>> data;
  size_t irow=1;
  for(string line; getline(cin,line);){
    vector<string> cells;
    istringstream issl(line);
    //// parse  row  header (row's numbers)
    string first_cell;
    getline(issl,first_cell,',');
    if(any_of(begin(mrows),end(mrows),[=](regex r){return regex_match(first_cell,r);})){
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
  }

  for(size_t i=0; i<cols.size(); ++i) {
    cout << cols[i];
    if(i!=cols.size()-1)
      cout << ",";
  }
  cout << endl;
  for(auto i:data){
    for(size_t j=0; j<i.size(); ++j) {
      cout << i[j];
      if(j!=i.size()-1)
        cout << ",";
    }
    cout << endl;
  }

  return 0;
}
