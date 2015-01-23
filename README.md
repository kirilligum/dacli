dacli
======

dacli (Data Analysis Command Line Interface) is a set of high performance tools to analyse data in the simplest way -- by piping commands to each other in a command line. 

# Why:
- **no learning curve**  -- use linux instead of learning new environments and languages, such as python/pandas, R, Octave, Matlab, ...
- **tired of waiting?** profiler-optimizaed cache-aware C++ and automatic switching between the latest algorithms; even avoiding expensive operations, such as division. Faster than numpy, matlab, R, ...
- **really big data?** --  one-pass online algorithms that minimize I/O, distribution of reading and computing the data, and convinient for clusters' scheduling
- **prototype quickly** = streams and component programming by using pipes `|` for intuitive work flow + less typing + leverage cli tools that you already know


# Tools:

**cuti** -- similar to `cut` but cuts the table based on its header and index into a subtable not just by integer value of index but by regex and ranges

    > cat test.csv|column -ts,
    id   aa    b    dd   eed    e    f
    1    0.3   1.2  3    123    1.2  43
    2    -21   1.3  0    43     7    4
    3    0.98  8.1  -34  .0923  23   1
    4    1.2   433  -89  232.3  12   0.01
    545  98    112  43   65     73   23


    > cat test.csv| ./cuti "1-3,.4.\*" "-.\*d,f"|column -ts,
    id   aa    b    dd   eed    f
    1    0.3   1.2  3    123    43
    2    -21   1.3  0    43     4
    3    0.98  8.1  -34  .0923  1
    545  98    112  43   65     23

**transpose** -- transposes the table. it simply swaps rows and columns

    > cat test.csv| ./transpose |column -ts,
    id   1    2    3      4      545
    aa   0.3  -21  0.98   1.2    98
    b    1.2  1.3  8.1    433    112
    dd   3    0    -34    -89    43
    eed  123  43   .0923  232.3  65
    e    1.2  7    23     12     73
    f    43   4    1      0.01   23


# Intallation:
compile using C++14 compiler. gcc 4.9.1 is fine. Ex: `g++ -std=c++14 transpose.cpp -o transpose`

