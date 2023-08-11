# ShinyFlourMCMC

## Example Work page: sensor parameterisation with ShinyFlourMCMC
![example work page](https://github.com/lynceuslq/ShinyFlourMCMC/blob/main/www/exampleworkpage.png)

## Example Work page: comparing metrics of genetic elements
![Compare radar](https://github.com/lynceuslq/ShinyFlourMCMC/blob/main/www/radar_compare.png)

****

## Implemented models for genetic gate parameterisation

-   Sensor
    $$f(I)=k*(\alpha + \frac{{I}^{n_\text{1}}}{{K}^{n_\text{1}} + {I}^{n_\text{1}}})$$
-   NOT Gate
    $$f(R_\text{3}) = k_\text{3}*(\alpha_\text{3} + \frac{{K_\text{3}}^{n_\text{3}}}{{K_\text{3}}^{n_\text{3}} + {R_\text{3}}^{n_\text{3}}})$$
-   AND Gate
    $$f(R,S) = G*\frac{{R/{K_\text{r}}}^{n_\text{r}}}{1+{R/{K_\text{r}}}^{n_\text{r}}}*\frac{{S/{K_\text{s}}}^{n_\text{s}}}{1+{S/{K_\text{s}}}^{n_\text{s}}}$$

## Bayesian Inference for genetic gate parameterisation (sensor example)
-   Gate model $$\mu_i = k*(\alpha + \frac{{I}^{n_\text{1}}}{{K}^{n_\text{1}} + {I_i}^{n_1}})$$
-   Model output variation $$Y_i \sim N(\mu_i,\sigma) , \sigma = \frac{1}{\sqrt \tau} $$
-   Parameter priors  
    $$k \sim \mathit{N(k_{exp}, k_{sd})}$$
    $$\alpha \sim \mathit{N(\alpha_{exp}, \alpha_{sd})}$$
    $$n1 \sim \mathit{N(n1_{exp}, n1_{sd})}$$
    $$K1 \sim \mathit{N(\alpha_{exp}, K1_{sd})}$$
    $$\tau \sim \mathit{\Gamma (0.1,0.1)}$$

## Workflow
![Workflow for sensor parameterisation](https://github.com/lynceuslq/ShinyFlourMCMC/blob/main/www/MCMC_workflow_single%20input.png)

****
## Hints for flourescent file uploading
![File input](https://github.com/lynceuslq/ShinyFlourMCMC/blob/main/www/fileinput.png)
