dacli
======

dacli (Data Analysis Command Line Interface) is a set of high performance tools to analyse data in the simplest way -- by piping commands within a command line. 

# Why:
- **no learning curve**  -- use Linux instead of learning new environments and languages.
- **fastest** -- profiler-optimized cache-aware C++ and automatic switching between the best up-to-date algorithms. Faster than NumPy, Matlab, R, and others.
- **really big data** --  one-pass online algorithms that minimize I/O and don't crush when data can't fit into the memory, parallel reading and computing, and convenient to schedule on a clusters.
- **prototype quickly** -- intuitive work flow + less typing + leverage CLI tools that you already know. Simply pipe '|' one command into another and don't worry about compatibility or control-flow.

# Tools:

**describe** -- mean, variance, quantiles, histogram, and other descriptive statistics using big-data-friendly algorithms that read data only once, scale linear, and use constant memory.

    > cat test2.csv| ./describe |column -ts,

    acc_name                          var1       var3        var7       var8        var9          
    count                             452061     452061      452061     452061      452061        
    mean                              2.80827    3.31845     3.73274    2.64713     0.105195      
    m2                                8.17182    11.6423     16.3645    8.58346     0.097163      
    m3                                24.5554    36.8962     81.3101    32.0809     0.101295      
    m4                                76.0191    241.388     444.087    142.307     0.225859      
    std                               0.534279   0.793856    1.5592     1.25544     0.293423      
    var                               0.285455   0.630207    2.43109    1.57614     0.086097      
    ske                               0.0224654  -11.8349    0.547793   0.512863    2.88804       
    kur                               0.418838   392.495     -0.342701  3.51789     22.5399       
    quant_0                           0.015398   -29.4863    -1.54412   -1.56047    -1.45811      
    quant_0.25                        2.46498    3           2.97866    1.99999     -6.34337e-13  
    quant_0.5                         2.86482    3.31815     3.15344    2.66765     1.48035e-13   
    quant_0.75                        3.06161    3.69929     4.99991    4           0.0216103     
    quant_1                           5.3155     6.35004     8.67222    43.3548     10.0039       
    hist_0                            54         10          4          20          12455         
    hist_1                            2411       19          17         440992      404686        
    hist_2                            3145       67          16232      11042       34838         
    hist_3                            72047      25          88280      4           51            
    hist_4                            136913     17          121862     1           8             
    hist_5                            183702     31          101598     0           1             
    hist_6                            30357      252         59113      0           0             
    hist_7                            23328      3344        29078      1           3             
    hist_8                            101        448293      34315      1           13            
    hist_9                            3          3           1518       0           6             
    plot_hist_0                       +          +           +          +           +             
    plot_hist_1                       +          +           +          ########+   ########+     
    plot_hist_2                       +          +           #+         +           #-            
    plot_hist_3                       ###+       +           ######+    +           +             
    plot_hist_4                       ######+    +           ########+  +           +             
    plot_hist_5                       ########+  +           #######-   +           +             
    plot_hist_6                       ##:        +           ####+      +           +             
    plot_hist_7                       #+         +           ##+        +           +             
    plot_hist_8                       +          ########+   ###:       +           +             
    plot_hist_9                       +          +           +          +           +             

Notes: (1) `quant_0` is `min`, `quant_0.5` is `median`, and `quant_1` is `max`. (2) Histograms are interpolated and approximate just to give you an idea of how the data is distributed.

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

**transpose** -- transposes the table. It simply swaps rows and columns

    > cat test.csv| ./transpose |column -ts,
    id   1    2    3      4      545
    aa   0.3  -21  0.98   1.2    98
    b    1.2  1.3  8.1    433    112
    dd   3    0    -34    -89    43
    eed  123  43   .0923  232.3  65
    e    1.2  7    23     12     73
    f    43   4    1      0.01   23
# Installation:
compile using C++14 compiler. gcc 4.9.1 is fine. Ex: `g++ -std=c++14 transpose.cpp -o transpose`

