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

inline int count_lines(char* buf, int start, int len){
  int count=0;
    for(char *p=buf+start;(p=(char*)memchr(p,'\n',(buf+start+len)-p));++p)
      ++count;
  return count;
}

int main() {
  cin.sync_with_stdio(false);
  //const int buf_size = 1024*16;
  const int buf_size = 1024*16*16*64*64;
  int count = 0;
  //int max_threads =  1;
  int max_threads =  std::thread::hardware_concurrency()-1;
  //int max_threads =  std::thread::hardware_concurrency();
  cout << "number of threads = " << max_threads << endl;
  vector<std::future<int>> vt(max_threads);
  vector<int> vc(max_threads,0);
  //int bytes_read;
  //char buf[buf_size+1];
  //int fd = STDIN_FILENO;
  //while((bytes_read=read(fd,buf,buf_size))>0){
    //char *p = buf;
    //while((p=(char*)memchr(p,'\n',(buf+bytes_read)-p))){
      //++p;
      //++count;
    //}
  //}
  for( char buf[buf_size+1]; int len = cin.read(buf,buf_size).gcount();){
    count += count_lines(buf,0,len);
    //int inc = len / max_threads;
    //for(int ith=0;ith<max_threads;++ith){
      //vt[ith] = std::async(count_lines,buf,ith*inc,inc);
    //}
    //for(int ith=0;ith<max_threads;++ith){
      //count += vt[ith].get();
    //}
  }
  cout << count<< endl;
  return 0;
}
