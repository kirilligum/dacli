#include <iostream>
//#include <iterator>
//#include <vector>
//#include <deque>
//#include <sstream>
//#include <string>
//#include <cmath>
//#include <cfloat>
//#include <stdexcept>
//#include <algorithm>
#include <numeric>
#include <unistd.h>
#include <cstring>
#include <thread>
#include <future>
#include <vector>

using namespace std;

//void count_lines(int& count){
//int count_lines(char* buf, int len){
void count_lines(char* buf, int len){
//void count_lines(int& count, char* buf, int len){
  int count=0;
  //for( char buf[buf_size+1]; int len = in.read(buf,buf_size).gcount();)
    for(char *p=buf;(p=(char*)memchr(p,'\n',(buf+len)-p));++p)
      ++count;
    cout << count << endl;
  //return count;
}

int count_lines_serial(char* buf, int len){
//void count_lines(int& count, char* buf, int len){
  int count=0;
  //for( char buf[buf_size+1]; int len = in.read(buf,buf_size).gcount();)
    for(char *p=buf;(p=(char*)memchr(p,'\n',(buf+len)-p));++p)
      ++count;
  return count;
}

int main() {
  cin.sync_with_stdio(false);
  const int buf_size = 1024*16*16*16*64;
  int count = 0;
  int max_threads =  std::thread::hardware_concurrency();
  cout << "number of threads = " << max_threads << endl;
  //int ith =0;
  //vector<std::future<int>> vt;
  vector<thread> vt;
  //vector<std::future<int>> vt(max_threads);
  vector<int> vc(max_threads,0);
  for( char buf[buf_size+1]; int len = cin.read(buf,buf_size).gcount();){
    //int cnt =0;
    //count += count_lines_serial(buf,len);
    vt.emplace_back(count_lines,buf,len);
    //vt.emplace_back(async(count_lines,buf,len));
    //std::future<int> f = std::async(count_lines,buf,len);
    //vc[ith] += f.get();
    //vt[ith] = std::async(count_lines,buf,len);
    //vc[ith] += vt[ith].get();
    //++ith;
    //ith = ith%max_threads;
    //if(vt[ith].joinable())
      //vt[ith].join();
    //count+=vc[ith];
    //vt.back().join();
    //vc.push_back(0);
  }
  for(auto&t:vt) t.join();
  //for(auto&f:vt) count+=f.get();
  //cout << std::accumulate(begin(vc),end(vc),0) << endl;
  //cout << count<< endl;
  return 0;
}
