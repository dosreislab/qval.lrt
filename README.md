# qval.lrt: Calculate q-values from LRT statistics

The package calculates q-values under the Benjamini and Hochberg (BH) and Storey
and Tibshirani (ST) false discovery rate (FDR) methods, using likelihood-ratio
test (LRT) statistics as input values. For one-sided LRTs, when the null model
is always true, the LRT statistic will be asymptotically distributed as a 50:50
mixture of a $\chi^2$ distribution and point mass zero (Self and Liang, 1987).
The resulting p-values will thus be distributed as a 50:50 mixture of a uniform
distribution betwen 0 and 1 and point mass at 1. The qval.lrt package
accommodates this distribution for correct calculation under the ST and BH
methods. Furthermore, under certain model misspecifications, the proportion of
LRTs equal to zero can be larger than 50%. The package implements a robust
version of the ST method (STR) to deal with this case. The package was originally 
developed for calculation of FDR tests in analysis of positive selection in DNA
sequences under the branch-site model (Yang and dos Reis, 2011).

# References

* Benjamini, Y. and Hochberg, Y. 1995, Controlling the false discovery rate: 
a practical and powerful approach to multiple testing. Journal of the Royal 
Statistical Society: Series B (Methodological) 57:289–300.

* Self, S. G. and Liang, K.-Y. 1987, Asymptotic properties of maximum likelihood 
estimators and likelihood ratio tests under nonstandard conditions. Journal of 
the American Statistical Association 82:605–610.

* Storey, J. D. and Tibshirani, R. 2003, Statistical significance for genomewide 
studies. Proceedings of the National Academy of Sciences 100:9440–9445.

* Yang, Z. and dos Reis, M. 2011. Statistical properties of the branch-site test 
of positive selection. Molecular Biology and Evolution, 28:1217–1228. 







