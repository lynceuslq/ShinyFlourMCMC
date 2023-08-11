# ShinyFlourMCMC

Bayesian Inference for genetic gate parameterisation (sensor example)
-   Gate model $$\mu_i = k*(\alpha + \frac{{I}^{n_\text{1}}}{{K}^{n_\text{1}} + {I_i}^{n_1}})$$
-   Model output variation $$Y_i \sim N(\mu_i,\sigma) , \sigma = \frac{1}{\sqrt \tau} $$
-   Parameter priors  
    $$k \sim \mathit{N(k_{exp}, k_{sd})}$$
    $$\alpha \sim \mathit{N(\alpha_{exp}, \alpha_{sd})}$$
    $$n1 \sim \mathit{N(n1_{exp}, n1_{sd})}$$
    $$K1 \sim \mathit{N(\alpha_{exp}, K1_{sd})}$$
    $$\tau \sim \mathit{\Gamma (0.1,0.1)}$$

## Workflow
![]
