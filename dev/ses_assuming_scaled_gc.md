Idea to estimate theta and sigma^2_v_plus_epsilon in the case where we assume that we know that g_c is scaled ($E[G_c]=0$, $\text{Var}(G_c)=1$).

## 1. Model & Assumptions

- **Model:** $G_{c,i} = \beta^T w_i + \epsilon_i$
- **Constraints:** $E[G_c]=0$ and $\text{Var}(G_c) = 1$.
- **Covariates:** $\Sigma_W = \text{Var}(w)$ is known.
- **Error:** $\epsilon \sim N(0, \sigma^2)$.
- **Implied identity:** The unit variance constraint implies $\sigma^2 = 1 - \beta^T \Sigma_W \beta$.

## 2. Point Estimation

We estimate the parameters $\theta = (\beta, \sigma^2)$ using a two-step "plug-in" approach:

1. **Estimate $\beta$ (OLS):**
   $$\hat{\beta} = (W^T W)^{-1} W^T G_c$$
   (Note: Since $\Sigma_W$ is known, asymptotically $(W^T W/n) \to \Sigma_W$, but standard OLS is used in finite samples.)

2. **Estimate $\sigma^2$ (Constrained):**
   $$\hat{\sigma}^2 = 1 - \hat{\beta}^T \Sigma_W \hat{\beta}$$

## 3. Joint Asymptotic Variance-Covariance Matrix

For the joint estimator $\hat{\theta} = (\hat{\beta}^T, \hat{\sigma}^2)^T$, the asymptotic distribution is:
$$\sqrt{n} \begin{pmatrix} \hat{\beta} - \beta \\ \hat{\sigma}^2 - \sigma^2 \end{pmatrix} \xrightarrow{d} N\left( 0, \mathbf{V} \right)$$

Where the covariance matrix $\mathbf{V}$ is:
$$\mathbf{V} = \begin{bmatrix}
\sigma^2 \Sigma_W^{-1} & -2\sigma^2 \beta \\
-2\sigma^2 \beta^T & 4\sigma^2 (1 - \sigma^2)
\end{bmatrix}$$

**Practical calculation:** To get the standard errors for your specific sample of size $n$, compute $\frac{1}{n} \mathbf{\hat{V}}$ by plugging in your estimates $\hat{\beta}$ and $\hat{\sigma}^2$ into the matrix above.
