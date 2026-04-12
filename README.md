# libscientific

[![DOI](https://joss.theoj.org/papers/10.21105/joss.05420/status.svg)](https://doi.org/10.21105/joss.05420)
[![Licence: GPL v3](https://img.shields.io/github/license/gmrandazzo/libscientific)](https://github.com/gmrandazzo/libscientific/blob/master/LICENSE)
[![Pylint](https://github.com/gmrandazzo/libscientific/actions/workflows/pylint.yml/badge.svg)](https://github.com/gmrandazzo/libscientific/actions/workflows/pylint.yml)
[![Pytest](https://github.com/gmrandazzo/libscientific/actions/workflows/pytest.yml/badge.svg)](https://github.com/gmrandazzo/libscientific/actions/workflows/pytest.yml)
[![Coverage](https://sonarqube.gmrandazzo.com/api/project_badges/measure?project=libscientific&metric=coverage&token=sqb_426f9a683b14ac981c8e32a1782672fc11c1a789)](https://sonarqube.gmrandazzo.com/dashboard?id=libscientific)
[![Maintainability Rating](https://sonarqube.gmrandazzo.com/api/project_badges/measure?project=libscientific&metric=sqale_rating&token=sqb_426f9a683b14ac981c8e32a1782672fc11c1a789)](https://sonarqube.gmrandazzo.com/dashboard?id=libscientific)

`libscientific` is a lightweight C framework for multivariate statistical analysis and numerical computing. Designed for high performance and portability, it is independent of most standard calculus libraries, making it ideal for both research environments and resource-constrained embedded systems.

## Key Features

- **Multivariate Analysis:** PCA (NIPALS), PLS, Consensus PCA (CPCA), Multiple Linear Regression (MLR), Fisher LDA, UPCA, and UPLS.
- **Clustering:** K-means++ and Hierarchical Clustering.
- **Instance Selection:** Most Descriptive Compounds (MDC) and Most Dissimilar Compounds (DIS).
- **Numerical Computing:** SVD, QR decomposition, cubic spline interpolation, and Nelder-Mead optimization.
- **Statistical Tools:** R², MSE, MAE, ROC curves, Precision-Recall, and cross-validation (Leave-One-Out, Bootstrap, Y-Scrambling).
- **Persistence:** Seamlessly save and load models using **SQLite** for robust, cross-platform storage.
- **Python Bindings:** Full-featured Python interface via `ctypes` for modern data science workflows.

## Quick Start (Python)

```python
from libscientific.pca import PCA
import libscientific.matrix as mx

# 1. Create dummy data (3 samples, 2 features)
data = [[1.0, 2.0], [2.0, 4.0], [3.0, 6.0]]

# 2. Fit a PCA model (1 component)
pca = PCA(scaling=1, npc=1)
pca.fit(data)

# 3. Get results
print(f"Explained Variance: {pca.get_exp_variance()}")
print(f"Scores: {pca.get_scores()}")
print(f"Loadings: {pca.get_loadings()}")

# 4. Save the model for later use
pca.save("my_pca_model.db")
```

## Quick Start (C)

```c
#include <scientific.h>

int main() {
    matrix *m;
    PCAMODEL *model;

    // 1. Initialize data (3 samples, 2 features)
    NewMatrix(&m, 3, 2);
    m->data[0][0] = 1.0; m->data[0][1] = 2.0;
    m->data[1][0] = 2.0; m->data[1][1] = 4.0;
    m->data[2][0] = 3.0; m->data[2][1] = 6.0;

    // 2. Fit PCA model (scaling type 1: unit variance, 1 component)
    NewPCAModel(&model);
    PCA(m, 1, 1, model, NULL);

    // 3. Print results and cleanup
    PrintPCA(model);
    
    DelPCAModel(&model);
    DelMatrix(&m);
    return 0;
}
```

**Compiling with C:**
```bash
gcc example.c -o example -lscientific -lm
```

## Installation

To use `libscientific`, you must first compile and install the C library on your system, then install the Python bindings.

### Prerequisites
- LAPACK/BLAS (or a Fortran compiler)
- SQLite3
- C Compiler (GCC or Clang)
- CMake

### Step 1: Build and Install the C Library
```bash
git clone https://github.com/gmrandazzo/libscientific.git
cd libscientific
mkdir build && cd build
cmake -DCMAKE_INSTALL_PREFIX=/usr/local ..
make -j
sudo make install
```

### Step 2: Install the Python Bindings
```bash
cd ../src/python_bindings
pip install .
```

### Alternative: macOS (Homebrew)
```bash
brew tap gmrandazzo/homebrew-gmr
brew install --HEAD libscientific
```

## Documentation & Resources
- **Full API Documentation:** [gmrandazzo.github.io/libscientific/](http://gmrandazzo.github.io/libscientific/)
- **Interactive Examples (Google Colab):**
  - [Drug Dataset Sampling](https://colab.research.google.com/drive/1gp8ppAsGlUbC4qGT-1Frc9ru1PiDvBl1) [![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)]
  - [Solubility Dataset PLS](https://colab.research.google.com/drive/1eQxLoZOrDMnTkxSkjuyTS_SuAyV3BcSF) [![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)]

## Development & Contributing
The library uses a specific data structure for every model, stored on the heap to support dynamic allocation. Every object (matrix, vector, tensor, model) must be manually allocated and deallocated using the `NewSOMETHING` and `DelSOMETHING` constructs.

**How to contribute:**
1. Fork the repository and update to the latest version.
2. Ensure your code is memory-leak free (validate with Valgrind).
3. Document all functions with parameters, attributes, and references.
4. Add a unit test in `src/tests` and a Python test in `src/python_bindings/tests`.
## Citation
If you use `libscientific` in your research, please cite the following paper:

> G.M. Randazzo, *libscientific: A C framework for multivariate and statistical analysis*, Journal of Open Source Software, 8(90), 5420, 2023. [https://doi.org/10.21105/joss.05420](https://doi.org/10.21105/joss.05420)

## Documentation & Resources
...
## References
1. P. Geladi, B.R. Kowalski, *Partial least-squares regression: a tutorial*, Analytica Chimica Acta, 185, 1986.
2. S. Wold et al., *Multi-way principal components and PLS-analysis*, Journal of Chemometrics, 1, 1987.
3. T. Kanungo et al., *An efficient k-means clustering algorithm*, IEEE Transactions on Pattern Analysis and Machine Intelligence, 2002.
4. B.D. Hudson et al., *Parameter Based Methods for Compound Selection*, Quantitative Structure-Activity Relationships, 15, 1996.
5. J. Holliday, P. Willett, *Definitions of "Dissimilarity" for Dissimilarity-Based Compound Selection*, Journal of Biomolecular Screening, 1, 1996.
6. R.D. Clark, P.C. Fox, *Statistical variation in progressive scrambling*, J Comput Aided Mol Des, 18, 2004.
7. J. A. Westerhuis et al., *Analysis of multiblock and hierarchical PCA and PLS models*, Journal of Chemometrics, 12, 1998.
## License
`libscientific` is distributed under the **GPLv3 License**. See the `LICENSE` file or [gnu.org/licenses/gpl-3.0.html](http://www.gnu.org/licenses/gpl-3.0.en.html) for details.
