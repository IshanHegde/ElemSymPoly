# RaschModel
A high performence python package for Rasch model analysis

## Algorithim Discription 

PROX is used to get the initial values for JMLE
### PROX

Initially, the mean abilties of all persons and the mean difficulties of all items are assumed to be zero. Then it is updated as follows:

Let $\beta_n$ be the the updated ability estimate for person $n$:


$$\beta_n = \mu_n + \sqrt{1+ \frac{\sigma_n^2}{2.9}}\log{\left(\frac{R_n}{N_n-R_n}\right)}$$

Where $\mu_n$ and $\sigma_n^2$ are the mean and variance of item difficulties faced by person $n$. $\R_n$ is the observed raw score of person $n$ and $\N_n$ is the maximum possible combined score of all items for person $n$.

Let $\delta_i$ be the updated difficulty estimate of item $i$:

$$\delta_i = \mu_i - \sqrt{1+\frac{\sigma_i^2}{2.9}}\log{\left(\frac{R_i}{N_i-R_i}\right)}$$ 

Where $\mu_i$ and $\sigma_i^2$ are the mean and variance of person abilities who attempted item $i$. Likewise, $\R_i$ is the the observed raw score of item $i$ and $\N_i$ is the maximum possible combined score of all persons for item $i$.

## Credits:

Linacre, Winsteps®, Iterations - PROX & JMLE, https://www.winsteps.com/winman/iterations.htm 

Wright, B. D. and Stone, M. H. (1979) Best test design. MESA Press: Chicago, IL
https://research.acer.edu.au/measurement/1/ Chapter 3, Section 6, pg.61
