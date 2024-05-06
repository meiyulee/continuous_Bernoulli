# Bivariate continuous trinomial distribution

## Model setting

$f(x,y)=C \times λ_{1}^x \times λ_{2}^y \times (1-λ_{1}-λ_{2})^{1-x-y}$

where $0<x<1, 0<y<1, 0<x+y<1$.

This distribution has two parameters, $λ_{1}, λ_{2}$, and the parameter space are
$0<λ_{1}<1, 0<λ_{2}<1$.

Let $f(x)=\int f(x,y)dy$. and $X$ **is not** Continuous Bernoulli distribution( $λ_{1}$ ).

Let $f(y)=\int f(x,y)dx$ , and $Y$ **is not** Continuous Bernoulli distribution( $λ_{2}$ ).

The marginal probability density function of $X$ and $Y$ have two parameters, $λ_{1}, λ_{2}$ and $\int f(x)dx= \int f(y)dy=1$.

$X$ and $Y$ **are not** independent random variables.

$X+Y$ is not Continuous Bernoulli distribution ( $λ_{1}, λ_{2}$ ).

The $(x,y,f(x,y))$ dynamic diagrams is affected by $λ_{1}$ and setting $λ_{1}+λ_{2}=c$.

This displayed method can understand $f(x,y)$ diagram changed when the $λ_{1}$ different value in simply.

1. $X$ ~ Continuous Bernoulli distribution ( $λ$ )


    $f(x)=C \times λ^x \times (1-λ)^{1-x}, 0<x<1$, 
    $\int f(x)dx=1$.
    Let $λ$ =0.01 to 0.99 and step = 0.01.

    Video
    
    https://github.com/meiyulee/continuous_Bernoulli/blob/master/_videos/Continuous_Bernoulli.mp4
