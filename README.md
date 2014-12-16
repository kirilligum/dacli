duckly
======

Command Line Interface toolkit for Data Analysis



Duckly (dacli == Data Analysis Command Line Interface) is a set of high performance tools to analyse data in the simplest way -- by piping commands to eachother in a command line. I follow linux philosophy, where I use a tool that integrates well with other toolS, for example Grep,cut, sort, paste, libsvm, gnuplot, and others that take text or csv file as an input including remote communication. Since all programs and packages have an option of csv file as an input and output, this toolkit allows component based programming, a style where components scan be easily replaced, and avoids the need of using api. 

The toolkit is made to do simple data manipulation quickly without use of variables, loading packages and so on. It is written in C++ to be high performance. The interface is intuitive due to resamblance with command line tools, such as cut, paste, grep. The command tools on their own are not that quick to use. they are not deisgned to work with headers, for example to sort based on collumn "age", you would need to find the index of that column, copy the header into a new file "final.csv", remove the header, sort and append to "final.xcsv". the operation is so common, that it should be a simple command, which is in our case is "cuti" (cut based on index or header). Also, there is no easy was to transpose a csv file, for example if you want to grep to find a column that contains a certain number. And what about arithmetic transformations or filtering based on inequalities. treatment of missing data? there needs to be a simple high performance way. I also plan to make a complementary C++, and may be C and D library with a similar interfac, so the algorithms can be easily rewritten as a code.


Limitations or why should you not use it:

This toolkit is simply data in, algorithm, data out, and combine data. for some complicated manipulations A language such as python might be better suited due to large collection of optimized libraries, that are still slow due to python's and cython's limitations; the same story is with java and scala. 
If you are comfortable making system calls from python interpriters, this kit is also not for you unless you want to oquicly manipulate data and get statistics without waiting for interpreter and importing pandas.
Another limitation is that conversion and copying takes time but also keep in mind that linux system keeps the files in ram so there is no difference of keeping the txt file loaded in interpreter or not. 

Intallation:
compile using C++14 compiler. C++11 should work too most of the time. gcc 4.9.1 is fine. `g++ -std=c++14 transpose.cpp -o transpose`

tools:

**transpose** -- transposes the table. it simply swaps rows and columns [done]

**cuti**-- similar to cut but cuts the matrix based on header and index into a submatrix not just by integer value of index but by regex and arithmetic expressions.

**cutd** -- similar to head,tail, sed -n '7p', grep with added arithmetics support. It builds a submatrix that contains matched elements with unmatched cells as missing values. Both headers, row header and column index, are kept.

**mergei** --merges columns or rows based on headers. As a spot, you can define the location where the mergin matrix would be placed.

**replace** -- replace values that are missing from another table or replace non-missing values. the shapes can be different

**sorti** -- sort based on a given header. basically it’s the same as sort but works on tables with headers.

**calc** -- arithmetic expresion (regex can be done with sed)

**stats** -- accumulaters. mean, std dev, median, higher moments, t-test, z-test, anova, and so on.

**gen** -- generators.  generate data from distribution, names headers can be supplied as a string or file.

**miss** -- missing value replacement using, value, mean median, regression, stochastic regression, expectation minimization (GMM), multiple imputation (bootstrap+regression), and may be also machine learning (in practice not useful due to long computational time)

**transform** -- transforms the data using normalization, one hot, indexing of catigorical, and so on.


**ml** -- machine learning random/extreame forest, gmm, svm, linear methods, gmm, pca, and other methods for classification, regression, density, clustering, variable importance transformation and so on. 

Examples of cli using core utils that are useful for data science:

===== Suppose we need to make a regression on data that has too many 0 as targets. We need a subset to test things and we want it to contain as many non-zero targets as zero targets. we simply need to sort the set on the target and head -n2000 (given that there are 1000 non-zero targets):

    cat -n transpose.cpp | tr , \\n | grep target | cut -f1

you’ll get the index of the target column let it be $x

    head -n1 train.csv > train_sorted.csv
    tail -n+2 train.csv|sort -t, -gk $x >> train_sorted.csv

as you can see, the problem is that for sort to be easy to work with it needs to deal well with headers. 

using duckly:
    sorti -gk”target” train.csv > train_sorted.csv


===== remove highly correlated varibales

the following line builds a correlation matrix with “stats corr”. “-n” displaces the names of correlated variables on each line with a threshold “0.96” starting from the variable at the diagonal of the correalation matrix. “-r” removes the first varibale in the line and lines with only one variables leaving only correlated variables; in addition, unique ran.

    cat train.csv| stats corr -rn “0.96” > correlated_names.csv


===== random forest regression:

    cat test.csv| ml -rf train.csv > predictions.csv
