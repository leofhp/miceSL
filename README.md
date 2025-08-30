# miceSL: Super Learning Methods for `mice`

`miceSL` provides a unified super learning method for multiple imputation by chained equations (MICE).

It includes the local imputation approach of Laqueur et al., 2022 (*superMICE*) and the matching approach of Carpenito & Manjourides, 2022 (*MISL*) as special cases.

## Installation

You can install the development version from GitHub:

```r
# install.packages("devtools")
devtools::install_github("leofhp/miceSL")
```

## Usage

The `strategy` and `bootstrap` arguments specify the sub-method to use:

- **`bootstrap`**:  
  Boolean. Determines whether the super learner is trained on the full data (`FALSE`) or on a bootstrap sample in each iteration (`TRUE`).  

- **`strategy`**:  
  Character string specifying how imputed values for continuous variables are generated from the trained super learner model:
  - `"localkernel"` uses the local normal approximation of Laqueur et al. with a kernel-estimated variance. 
  - `"matching"` uses the predictive mean matching approach of Carpenito & Manjourides.
  - `"dirichlet"` perturbs predictions by sampling new super learner weights from a Dirichlet distribution centered on the optimal weights. 

Special cases:
- `strategy = "localkernel", bootstrap = FALSE` → **superMICE**  
- `strategy = "matching", bootstrap = TRUE` → **MISL**

---

## Imputation for Different Variable Types

- `strategy` only affects **continuous variables**.  
- Missing values in **binary** or **categorical** variables are always imputed by random draws from binomial or categorical distributions centered around the super learner predictions.  
- The `outcome_type` argument controls this and can be `"continuous"`, `"binary"`, or `"categorical"`.  

Defaults:
- Numeric 0/1 or factor with two levels → `"binary"`  
- Factor with more than two levels → `"categorical"`  
- Otherwise → `"continuous"`

---

## Hyperparameters

### Local kernel approach (`strategy = "localkernel"`)
- **`kernel`**: Choice of `"gaussian"` (default), `"uniform"`, or `"triangular"`.  
- **`bandwidth`**: Should match the data scale.  
  - Laqueur et al. suggest setting it so that 0.7–2% of observed values lie within one standard deviation of the prediction under the Gaussian kernel. 
  - Alternatively, they propose a jackknife selection procedure, which is called by passing a numeric vector of candidate bandwidths.  
- Note: the `bw.update` argument from Laqueur et al. is *not* included due to the error-prone environment handling required for maintaining bandwidth configurations across iterations.

### Dirichlet approach (`strategy = "dirichlet"`)
- **`alpha`**: Scales the variance of the Dirichlet distribution used to perturb weights.  
  - Larger `alpha` → lower variance of imputations.  
  - Smaller `alpha` → higher variance.  

---

## Super Learner Setup

The base learner set is specified via:

- `continuous_learners`  
- `binary_learners`  
- `categorical_learners`  

Recommendations:
- Choose a more **diverse learner set** than the defaults.  
- Any algorithm from **`SuperLearner`** or **`sl3`** can be used, depending on `sl_engine`.  

Further arguments:
- **`sl_engine`**:
    - Defaults to `sl_engine = "SuperLearner"`, as it tends to have a faster runtime on a single core, which is the relevant performance measure since parallelization is recommended at the `mice` level (e.g., via `mice::futuremice`) rather than within the super learner itself.
    - However, `outcome_type = "categorical"` is currently only supported with `sl3`.  
- **`cv_folds`**:
    - Number of cross-validation folds (default: `10`). Can be reduced (e.g., to 5) to save computation. The impact of this is likely small since predictions are perturbed during imputation, reducing the risk of overfitting.   

---

## Example

Here is a simple example using `miceSL`:

```r
# Bivariate data (x, y) with missing values in x based on an MAR mechanism.
N <- 50
p <- 0.2
x <- rnorm(n = N, mean = 3, sd = 1)
y <- 10 + 2 * sin(x^2) - 0.05 * exp(x) - 0.5 * x +
  rnorm(n = N, mean = 0, sd = 0.5)
x[(runif(N) < (2 * p * pnorm(y, mean = mean(y), sd = sd(y))))] <- NA
df <- data.frame(x, y)

imps <- mice::mice(
  df, m = 50, maxit = 1, method = "sl", outcome_type = "continuous",
  continuous_learners = c("SL.mean", "SL.glm"), strategy = "localkernel",
  bootstrap = FALSE, bandwidth = 0.5, cv_folds = 5
)
```

---

## References

- Laqueur HS, Shev AB, Kagawa RMC. SuperMICE: An Ensemble Machine Learning Approach to Multiple Imputation by Chained Equations. American Journal of Epidemiology. 2022;191(3):516-525.
- Carpenito T, Manjourides J. MISL: Multiple imputation by super learning. Statistical Methods in Medical Research. 2022;31(10):1904-1915.