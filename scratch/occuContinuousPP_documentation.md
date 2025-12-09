

# Likelihood of occupancy with a Poisson Process detection process

Let's note the following:

- $i\in{1,\dots, \mathcal{S}}$ sites
- $j \in 1, \ldots, \mathcal{J_i}$ deployments per site
- $k \in 0, 1, \ldots, \mathcal{K_ij}$ detection events per deployment. 
  - $t_{ijk}$ the time of the $k^{th}$ detection event in deployment $j$, site $i$
  - $\mathcal{K}_i = \sum_{j=1}^{\mathcal{J}_{i}} (\mathcal{K}_{ij})$ the number of detection events in site $i$
  - $\mathcal{S}_+$ the sites with at least one detection ($\mathcal{K}_i>0$) and  $\mathcal{S}_-$  the sites with no detection  ($\mathcal{K}_i=0$)

With the following process:

- Each site is associated with an unobserved occupancy state ($Z_i=0$ if unoccupied, $Z_i=1$ if occupied))
  - $Z_i\sim\mathrm{Bernoulli}(\psi_i)$
  - $\operatorname{logit}(\psi_i) = \alpha_0 + \alpha_1  \text{siteCov}_{1i} + \alpha_2  \text{siteCov}_{2i} + \ldots$
- Conditional on occupancy ($Z_i=1$), detection events arise from a Poisson process during each deployment.
  - Poisson process with intensity function $\lambda_j(t)$
  - $\log \left( \lambda_j(t) \right) = \beta_0 + \beta_1 obsCov_{1j}(t) + \beta_2 obsCov_{2j}(t) + \cdots$
  - We note $\Lambda_{ij}=\int_{T_0}^{T_1}\lambda(t) \, dt$ the integral of the intensity function during deployment $j$ in site $i$, with $T_0$ and $T_1$ the beginning and ending times of the deployment, respectively. It can also be thought of as the mean expected number of detection events if the site is occupied.

The probability of obtaining the data can be declined into two: for $\mathcal{S}_+$ and $\mathcal{S}_-$

> Note: the likelihood of the inhomogeneous Poisson Process **with ordered events** is:  $L(\lambda;data)=e^{-\Lambda}\prod_{i=1}^{n}\lambda(x_i)$
>
> - Drazek, Ludwik Czeslaw. 2013. Intensity Estimation for Poisson Processes. https://teaching.seehuhn.de/projects/2012-Poisson/Drazek.pdf. (*Chapter 5; Section 5.2. Maximum Likelihood Estimation of the intensity*). 
> - Ross, Sheldon M. 1996. Stochastic Processes. (p66-67) 
>   *if the sample is ordered then the likelihood of the ordered sample is the likelihood of the unordered sample above multiplied by $n!$*

- For a given site in $\mathcal{S}_+$, we know that the site is occupied

$$
\mathbb{P}(\text{this data}, \mathcal{K}_i>0) = 

\overbrace{
\overbrace{\psi_i}^{\text{Occupied}} 
\times
\overbrace{
    \prod_{j=1}^{\mathcal{J}_i} \left(
    	e^{-\Lambda_{ij}}
    	\prod_{k=1}^{\mathcal{K}_{ij}}
    	{\lambda_j(t_{ijk})}    \right)
}^{\text{Observed}}
}^{\text{Occupied and observed}}
$$



- For a given site in $\mathcal{S}_-$, we have two possibilities: either the site is occupied and we observed no detection event, or the site is not occupied

$$
\mathbb{P}(\text{this data}, \mathcal{K}_i=0) = 

\overbrace{
\overbrace{\psi_i}^{\text{Occupied}} 
\times
\overbrace{
    \prod_{j=1}^{\mathcal{J}_i} \left(
    	e^{-\Lambda_{ij}}
    \right)
}^{\text{Unobserved}}
}^{\text{Occupied and unobserved}}
+
\overbrace{(1-\psi_i)}^{\text{Not occupied}}
$$

## Expressed directly from observation covariates functions

We have: $\log \left( \lambda_j(t_{ijk}) \right) = \beta_0 + \beta_1 obsCov_{1ij}(t_{ijk}) + \beta_2 obsCov_{2ij}(t_{ijk}) + \cdots$

So, we have: 

$$
\lambda_j(t_{ijk}) = \exp \left( X_{ij}(t_{ijk})^\top \beta \right)
$$

- With $X_j(t)$ the covariates' vector (including a constant vector of 1s)
- And $\beta$ the coefficients' vector (including the intercept)

And, we have:

$$
\begin{align}
\Lambda_{ij}
&= \int_{T_0}^{T_1} \lambda(t) \, dt \\
&= \int_{T_0}^{T_1} \exp \left( X_{ij}(t)^\top \beta \right) \, dt
\end{align}
$$

































































