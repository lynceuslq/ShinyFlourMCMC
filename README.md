# ShinyFlourMCMC
ShinyFlourMCMC is a Shiny App developed to paramterise the metrics of genetic gates and sensors with inflourescence datasets. ShinyFlourMCMC is a updated version of ShinyFlour and can use Bayesian inference and Markov Chain Monte Carlo algorithms (more precisely, Hamilton Monte Carlo) to assess the performance of genetic elements, which is proved to be a more statistically robust method here. Here is also an online version (https://lynceuslq.shinyapps.io/shinyflour_new/) of ShinyFlourMCMC you can use, yet using Bayesian inference may cause timeout in connection, it is recommended to build the App on your server with source codes if you want to use Bayesian and MCMC methods.

## Example Work page: sensor parameterisation with ShinyFlourMCMC
![example work page](https://github.com/lynceuslq/ShinyFlourMCMC/blob/main/www/exampleworkpage.png)

## Example Work page: comparing performance metrics of genetic elements
![Compare radar](https://github.com/lynceuslq/ShinyFlourMCMC/blob/main/www/radar_compare.png)

****

## Implemented steady-state models for genetic gate parameterisation

-   Sensor
    $$f(I) \sim k*(\alpha + \frac{{I}^{n_\text{1}}}{{K}^{n_\text{1}} + {I}^{n_\text{1}}})$$
-   NOT Gate
    $$f(R_\text{3}) \sim k_\text{3}*(\alpha_\text{3} + \frac{{K_\text{3}}^{n_\text{3}}}{{K_\text{3}}^{n_\text{3}} + {R_\text{3}}^{n_\text{3}}})$$
-   AND Gate
    $$f(R,S) \sim G*\frac{{R/{K_\text{r}}}^{n_\text{r}}}{1+{R/{K_\text{r}}}^{n_\text{r}}}*\frac{{S/{K_\text{s}}}^{n_\text{s}}}{1+{S/{K_\text{s}}}^{n_\text{s}}}$$

## Bayesian Inference for genetic gate parameterisation (sensor example)
-   Gate model $$\mu_i = k*(\alpha + \frac{{I}^{n_\text{1}}}{{K}^{n_\text{1}} + {I_i}^{n_1}})$$
-   Model output variation $$Y_i \sim N(\mu_i,\sigma) , \sigma = \frac{1}{\sqrt \tau} $$
-   Parameter priors: priors are set from the results NLS optimisation
    $$k \sim \mathit{N(k_{exp}, k_{sd})}$$
    $$\alpha \sim \mathit{N(\alpha_{exp}, \alpha_{sd})}$$
    $$n1 \sim \mathit{N(n1_{exp}, n1_{sd})}$$
    $$K1 \sim \mathit{N(\alpha_{exp}, K1_{sd})}$$
    $$\tau \sim \mathit{\Gamma (0.1,0.1)}$$

## Workflow

### Single-input gate modeling workflow
![Workflow for sensor parameterisation](https://github.com/lynceuslq/ShinyFlourMCMC/blob/main/www/MCMC_workflow_single%20input.png)

### 2-input gate (AND gate) workflow
![2-input](https://github.com/lynceuslq/ShinyFlourMCMC/blob/main/www/2inputworkflow.png)

****
## Hints for flourescent file uploading
![File input](https://github.com/lynceuslq/ShinyFlourMCMC/blob/main/www/fileinput.png)
