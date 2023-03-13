---
title: 'Libscientific: A Powerful C Library for Multivariate Analysis'
tags:
  - C
  - python
  - chemometrics
  - multivariate analysis
authors:
  - name: Giuseppe Marco Randazzo
    orcid: 0000-0003-1585-0019
    affiliation: "1"
affiliations:
 - name: independent researcher
   index: 1
date: 17 March 2023
bibliography: paper.bib

---

# Summary

Multivariate analysis is a powerful technique that allows researchers to analyze and interpret data with multiple variables.
In today's data-driven world, multivariate analysis has become essential for the exploration of complex data sets.
Libscientific is a powerful library written in C that provides a comprehensive set of multivariate analysis tools.
The library includes several multivariate analysis algorithms, such as principal component analysis (PCA), partial least squares regression (PLS), consensus principal component analysis (CPCA), multiblock principal component analysis, and multiblock partial least squares (UPLS). Libscientific also includes several other tools to analyze data, such as cluster analysis using KMeans Hierarchical clustering and other methods to run linear algebra calculations. The library also provides a python foreign function to be used inside python scripts.



# Statement of need

The library is designed to be easy to use and can be integrated into any C or C++ project. Additionally, libscientific comes with a foreign function Python bindings, making it accessible within Python scripts and easier to perform data analysis tasks.
One of the main advantages of libscientific is its performance. Because the library is written in C, it is highly optimized for performance. This means that large data sets can be analyzed quickly and efficiently, making it an ideal choice for applications where speed is critical.
The library depends only on lapack for SVD and eigenvalues decomposition and can be easily integrated into embedded systems.
The current library version is 1.4.1, and here is a list of the current library features:

* Principal Component Analysis (PCA)
* Consensus Principal Component Analysis (CPCA)
* Partial Least Squares (PLS)
* Multiple Linear Regression (MLR)
* Unfold Principal Component Analysis (UPCA)
* Unfold Partial Least Squares (UPLS)
* Fisher Linear Discriminant Analysis (LDA)
* Kmeans++ Clustering
* Hierarchical Clustering
* Sample selection algorithms: Most Descriptive Compound (MDC), Most Dissimilar Compound (MaxDis)
* Statistical measures: R2, MSE, MAE, RMSE, Sensitivity, PPV
* Yates Analysis
* Receiver Operating Characteristic curve anaysis (ROC)
* Precision-Recal curve analysis
* Matrix-matrix Euclidean, Manhattan, Cosine and Mahalanobis distances
* Numerical integration
* Natural cubic spline interpolation and prediction
* Linear algebra (Eigenvector/value and SVD operated by Lapack library)
* Ordinary Least Squares solver
* Linear equation Solver
* Nelder-Mead Simplex Optimization
* Cross validation methods: Bootstrap k-fold, Leave-One-Out, Y-Scrambling

Libscientific was designed to analyze multivariate chemical data in the general cheminformatic and -omic fields.
It has already been used indirectly by researchers in scientific publications [@Randazzo16;@Randazzo171;@Randazzo172;@Randazzo20;@Kwon21;@Kwon22], and it is also used in student courses on general multivariate analysis.


# Multivariate analyisis algorithms specs

Principal component analysis (PCA) is one of the most commonly used methods for multivariate analysis. PCA is an unsupervised method that compresses data into low-dimensional representations that capture the dominant variation in the data. Libscientific provides a robust implementation of PCA using the NIPALS algorithm described in Geladi paper [@Geladi86]. libscientific implementation can handle data sets with many variables, few instances, and missing values.

Partial least squares (PLS) is another commonly used method for multivariate analysis. PLS is a supervised method that captures the dominant covariation between the data matrix and the target/response. Libscientific provides one version of PLS described by Geladi et al [@Geladi86]. This implementation works with single-task and multi-task regression problems.

In addition to PCA and PLS, libscientific provides implementations of Consensus PCA (CPCA) to analyze time series and multi-block data, algorithm described in the Westerhuis paper [@Westerhuis98], and other multi-block methods such as Unfold Principal Component Analysis (UPCA) and Unfold Partial Least Squares (UPLS) both implementation from Wold et al[@SWold87].

All multivariate algorithms admit missing values since the core linear algebra functions are coded to skip missing values, according to Martens et al (page 381)[@Martens2001].

# Other algorithms

The library also provides compound selection algorithms such as Most Descriptive Compounds [@Hudson96] or Most Dissimilar Compound [@Holliday1996] selections, allowing one to analyze scores plots or original data matrices and select samples based on the object/sample diversity.
Moreover, multi-thread cross-validation methodologies such as "Bootstrap k-fold" Leave-One-Out (LOO), and Y-Scrambling tests are implemented to facilitate the scientist in testing model prediction abilities.

# Conclusions

Libscientific is a powerful library that provides a comprehensive set of multivariate analysis tools for researchers and analysts. Whether a scientists work on research or data analytics, libscientific can help gain deeper insights into the data. Its C-based implementation and Python bindings offer high performance and ease of use, making it an ideal choice for data-driven applications.


# Acknowledgements

Libcientific was born as an open-source project from the Ph.D. thesis of the author Giuseppe Marco Randazzo.
The author acknowledges the support from the University of Perugia, the valuable code review made by the people from Freaknet Medialab, and the bug reports from the whole open-source community using this library.


# References
