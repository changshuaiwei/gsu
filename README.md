# gsu

DESCRIPTION: Generalized Similarity U test

AUTHOR: Changshuai Wei

CONTACT: weichangshuai@gmail.com

YEAR: 2016

LICENSE: Released under GNU General Public License, v2 (see
COPYING.txt)

DOCUMENTATION: http://changshuaiwei.github.io/software.html


COMPILATION: You will need a standard C++ compiler such as GNU g++
(version 3) and additional Eigen and Boost library.

Using GSU:

./gsu --gtu --bfile segment --setc gene_set.txt --cov cov.PC.txt --phs multi_phenotype.txt --phs-wt weight.txt --LK --pjG --out assoc.multi.wt.rst


./gsu --gtu --bfile segment --setc gene_set.txt cov.PC.txt --phs uni_phenotype.txt --LK --pjG --out assoc.uni.wt.rst
