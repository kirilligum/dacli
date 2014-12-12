duckly
======

Command Line Interface toolkit for Data Analysis



Duckly (dacli == Data Analysis Command Line Interface) is a set of high performance tools to analyse data in the simplest way -- by piping commands to eachother in a command line. I follow linux philosophy, where I use a tool that integrates well with other toolS, for example Grep,cut, sort, paste, libsvm, gnuplot, and others that take text or csv file as an input including remote communication. Since all programs and packages have an option of csv file as an input and output, this toolkit allows component based programming, a style where components scan be easily replaced, and avoids the need of using api. 

The toolkit is made to do simple data manipulation quickly without use of variables, loading packages and so on. It is written in C++ to be high performance. The interface is intuitive due to resamblance with command line tools, such as cut, paste, grep. The command tools on their own are not that quick to use. they are not deisgned to work with headers, for example to sort based on collumn "age", you would need to find the index of that column, copy the header into a new file "final.csv", remove the header, sort and append to "final.xcsv". the operation is so common, that it should be a simple command, which is in our case is "cuti" (cut based on index or header). Also, there is no easy was to transpose a csv file, for example if you want to grep to find a column that contains a certain number. And what about arithmetic transformations or filtering based on inequalities. treatment of missing data? there needs to be a simple high performance way. 

Limitations or why should you not use it:
This toolkit is simply data in, algorithm, data out, and combine data. for some complicated manipulations A language such as python might be better suited due to large collection of optimized libraries, that are still slow due to python's and cython's limitations; the same story is with java and scala. 
If you are comfortable making system calls from python interpriters, this kit is also not for you unless you want to oquicly manipulate data and get statistics without waiting for interpreter and importing pandas.
Another limitation is that conversion and copying takes time but also keep in mind that linux system keeps the files in ram so there is no difference of keeping the txt file loaded in interpreter or not. 

tools:

**transpose** -- transposes the table. it simply swaps rows and columns [done]

**cuti**-- similar to cut but cuts the matrix based on header and index into a submatrix not just by integer value of index but by regex and arithmetic expressions.

**cutd** -- similar to head,tail, sed -n '7p', grep with added arithmetics support. It builds a submatrix that contains matched elements with unmatched cells as missing values. Both headers, row header and column index, are kept.

**merge** --merges columns or rows based on headers. As a spot, you can define the location where the mergin matrix would be placed.

**replace** -- replace values that are missing from another table or replace non-missing values. the shapes can be different

**calc** -- arithmetic expresion (regex can be done with sed)

**accu** -- accumulaters. mean, std dev, median, higher moments, t-test, z-test, anova, and so on.

**gen** -- generators.  generate data from distribution, names headers can be supplied as a string or file.

**miss** -- missing value replacement using, value, mean median, regression, stochastic regression, expectation minimization (GMM), multiple imputation (bootstrap+regression), and may be also machine learning (in practice not useful due to long computational time)

**ml** -- machine learning random/extreame forest, gmm, svm and other methods for classification, regression, density, clustering and so on. 


