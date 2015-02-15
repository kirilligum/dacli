#pragma once
#include <iostream>
#include <vector>
#include <cmath>

using namespace std;

template<typename C>
class counter{
  public:
    vector<C> c;
    counter():c{0},max{(pow(2,8*sizeof(c[0])))}{}
    void increment(){
      c[0]++;
      for(size_t i=1;c[i-1]==0;++i){
        c[i]++;
        if(i==c.size())
          c.push_back(1);
      }
    }
    vector<C> value(){
      return c;
    }
    double double_value(){
      double x=0.0;
      for(size_t i=0;i<c.size();++i){
        x+=static_cast<double>(c[i])*pow(max,i);
      }
      return x;
    }
    double max;
};

