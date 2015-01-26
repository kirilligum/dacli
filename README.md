dacli
======

dacli (Data Analysis Command Line Interface) is a set of high performance tools to analyse data in the simplest way -- by piping commands within a command line. 

# Why:
- **no learning curve**  -- use linux instead of learning new environments and languages.
- **tired of waiting?** profiler-optimizaed cache-aware C++ and automatic switching between the latest algorithms. Faster than numpy, matlab, R, ...
- **really big data?** --  one-pass online algorithms that minimize I/O and don't crush due to lack of memory, parallel reading and computing, and convenient to schedule on a clusters.
- **prototype quickly** = intuitive work flow + less typing + leverage cli tools that you already know. Simply pipe '|' once command into another and don't worry about compatibility or control-flow.

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

**describe** -- mean, variance, quantiles, histogram, and other descriptive statistics using big-data-friendly algorithms  that read data only once, scale linear and use constant memory.

    > cat test.csv| ./describe |column -ts,

**transpose** -- transposes the table. it simply swaps rows and columns

    > cat test.csv| ./transpose |column -ts,
    id   1    2    3      4      545
    aa   0.3  -21  0.98   1.2    98
    b    1.2  1.3  8.1    433    112
    dd   3    0    -34    -89    43
    eed  123  43   .0923  232.3  65
    e    1.2  7    23     12     73
    f    43   4    1      0.01   23
**describe** -- read the data once and outputs descriptive statistics. Note, quantiles and histograms are approximate.
# Intallation:
compile using C++14 compiler. gcc 4.9.1 is fine. Ex: `g++ -std=c++14 transpose.cpp -o transpose`

