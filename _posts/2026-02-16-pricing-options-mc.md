---
title: 'Monte Carlo Pricing of Barrier Options '
date: 2026-02-12
permalink: /posts/2026/02/pricing_dki_mc/
tags:
  - technical
---





# Monte Carlo Pricing of Barrier Options  

From Constant Volatility to Local Volatility, Stochastic Volatility, and Jumps

This write-up explains **what we are doing** and **why we are doing it** when pricing barrier options using Monte Carlo (MC).

---

# A. Foundations

## A1. Why Monte Carlo pricing? When do analytical models fail?

In risk-neutral pricing, the value of a derivative is the discounted expectation of its payoff under the **risk-neutral measure** $Q$:

$$
V_0 = e^{-rT}\,\mathbb{E}^Q[\text{Payoff}]
$$

Analytical (closed-form) solutions exist only when the model and payoff are simple enough that we can evaluate the expectation exactly. Examples:

- **Black–Scholes** gives closed forms for European calls/puts under constant volatility.
- **Heston** has semi-closed forms for European options via Fourier methods.
- **PDE** methods can price many products in low dimensions (typically 1–2 state variables).

Analytical methods start to fail (or become unwieldy) when:

### 1) Payoffs are path-dependent
Barrier options depend on the entire path:

- Down-and-in put payoff depends on whether $S_t$ ever crossed a barrier.

$$
\mathbf{1}\left\{\min_{0\le t\le T} S_t \le B\right\}
$$

This indicator makes the pricing expectation more complex than standard European payoff $f(S_T)$.

### 2) The model contains additional state variables
Examples:
- stochastic volatility adds $v_t$
- local volatility makes $\sigma$ a function of $(t,S)$
- jumps make the process discontinuous

Each added feature increases the complexity of either:
- solving a PDE, or
- finding a closed-form transform

### 3) High dimensionality
Baskets, multi-asset barriers, worst-of, correlation products, etc. quickly become high-dimensional. PDE methods suffer from the curse of dimensionality.

### Why MC works in all these cases
MC does not require closed forms. It only requires the ability to simulate sample paths under the model:

$$
V_0 \approx e^{-rT}\,\frac{1}{N}\sum_{i=1}^N \text{Payoff}^{(i)}
$$

Its power is generality: if you can simulate the model, you can price the payoff.

---

## A2. Assumptions in Monte Carlo pricing

### (1) We assume a *risk-neutral* model is specified
We must choose dynamics for $S_t$ under $Q$. For example:

- Constant vol (GBM)
- Local vol: $\sigma_{\text{loc}}(t,S)$
- Heston: stochastic variance $v_t$
- Merton: jumps + diffusion

The MC estimate is only as meaningful as the model.

### (2) We assume we can discretize time
We simulate on a grid:

$$
0 = t_0 < t_1 < \cdots < t_n = T,\quad \Delta t = T/n
$$

Even if the true process is continuous-time, MC runs in discrete time.

Some processes (GBM, Merton log-price over a step) have **exact** step simulation. Others require numerical schemes (Euler, Milstein, QE).

### (3) We assume the law of large numbers and CLT
As the number of paths $N \to \infty$:

$$
\hat V_N \to V_0
$$

and asymptotically the statistical error scales like:

$$
\text{SE}(\hat V_N) \propto \frac{1}{\sqrt{N}}
$$

This is the key limitation: doubling accuracy requires roughly **4×** more paths.

---

# B. Constant Volatility Case (GBM)

We start with the baseline Black–Scholes / GBM model under $Q$:

$$
dS_t = (r-q)S_t\,dt + \sigma S_t\,dW_t
$$

Barrier product of interest (Down-and-In Put):

$$
\text{Payoff} = \mathbf{1}\left\{\min_{0\le t\le T} S_t \le B\right\}\,(K-S_T)^+
$$

---

## B1. Discretization and simulation

### Log transform
Define $X_t = \ln S_t$. Using Itô’s lemma:

$$
dX_t = \left(r-q-\frac12\sigma^2\right)dt + \sigma\,dW_t
$$

Over a step $\Delta t$, this has an **exact** normal increment:

$$
X_{t+\Delta t} = X_t + \left(r-q-\frac12\sigma^2\right)\Delta t + \sigma\sqrt{\Delta t}\,Z,\quad Z\sim N(0,1)
$$

So for GBM, the simulation of endpoints is exact—there is no discretization bias in $S$ itself.

### Where discretization bias still appears: barrier monitoring
Even if $S_t$ endpoints are exact, barriers depend on the *minimum over continuous time*. If we only check $S_{t_i}$ on the grid, we may miss barrier hits between times.

This becomes the major source of bias for barrier options.

---

## B2. Variance reduction (why and how)

The MC estimator is:

$$
\hat V = e^{-rT}\frac{1}{N}\sum_{i=1}^N X_i
$$

with variance:

$$
\text{Var}(\hat V) = \frac{\text{Var}(X)}{N}
$$

Variance reduction aims to reduce $\text{Var}(X)$ without changing the expectation.

### 1) Antithetic variates
For each path using normals $Z$, we also simulate a path using $-Z$. For many payoffs, this induces negative correlation between paired payoffs, lowering variance.

Intuition:
- If $Z$ creates an unusually high $S_T$, then $-Z$ tends to create an unusually low $S_T$.
- Averaging the two reduces random fluctuations.

---

## B3. Brownian bridge correction (continuous barrier monitoring)

### Problem: discrete monitoring bias
If we monitor only at grid points:

$$
\min_i S_{t_i} \le B
$$

we miss events where the true continuous path crosses $B$ between $t_i$ and $t_{i+1}$.

This typically underestimates barrier hits and biases barrier option prices.

### Brownian bridge idea
Work in log space $X_t = \ln S_t$. Consider one time step $[t,t+\Delta t]$.

We know the endpoints:

$$
X_t=x,\quad X_{t+\Delta t}=y
$$

Conditioned on endpoints, the path behaves like a **Brownian bridge**. Under constant volatility, the probability the bridge crosses a lower barrier $b=\ln B$ is known:

- If $\min(x,y)\le b$, barrier hit is certain: $p_{\text{hit}}=1$
- Otherwise:

$$
p_{\text{hit}}
=
\exp\left(
-\frac{2(x-b)(y-b)}{\sigma^2\Delta t}
\right)
$$

### How we use it in MC
For each step, compute $p_{\text{hit}}$. Sample $U\sim \text{Uniform}(0,1)$. If $U<p_{\text{hit}}$, mark the path as having hit the barrier in that interval.

This approximates **continuous monitoring** without shrinking $\Delta t$.

Key insight:
- We don’t simulate the exact crossing time.
- We simulate whether a crossing occurred based on the correct conditional probability.

---



# C. Local Volatility Model

## C1. What is local volatility?

Local vol assumes the underlying follows:

$$
dS_t = (r-q)S_t\,dt + \sigma_{\text{loc}}(t,S_t)\,S_t\,dW_t
$$

Here $\sigma_{\text{loc}}$ depends on both time and level.

This is a **diffusion model** (continuous paths), but with state-dependent volatility.

---

## C2. Why constant volatility “sucks” in practice

Black–Scholes implies a flat implied volatility surface:
- implied vol does not depend on strike or maturity

But markets show:
- **skew** (vol varies with strike)
- **smile** (curvature)
- **term structure** (vol varies with maturity)

A constant $\sigma$ model cannot match observed option prices across strikes/maturities simultaneously.

---

## C3. How do we get local vol from market quotes?

The market provides implied vol quotes on a grid:

$$
\sigma_{\text{imp}}(T,K)
$$

But local vol is not equal to implied vol. Implied vol is the BS volatility that reproduces today’s option price for a specific $(T,K)$.

Pipeline:

1. Collect market implied vol quotes $\sigma_{\text{imp}}(T_i,K_j)$
2. Convert each quote into a call price using BS:

$$
C(T_i,K_j) = BSCall(S_0,K_j,T_i,r,q,\sigma_{\text{imp}}(T_i,K_j))
$$

3. Smooth the call price surface (we must ensure no-arbitrage constraints)
4. Apply **Dupire** to compute $\sigma_{\text{loc}}(T,K)$

---

## C4. Dupire formula (and what it means)

Dupire provides local variance from the call surface:

$$
\sigma_{\text{loc}}^2(T,K)
=
\frac{
\frac{\partial C}{\partial T}
+
(r-q)K\frac{\partial C}{\partial K}
+
qC
}{
\frac12 K^2\frac{\partial^2 C}{\partial K^2}
}
$$

### Why we need strike *and* maturity dimensions
- $\partial C/\partial T$ requires variation across maturities ⇒ term structure matters
- $\partial^2C/\partial K^2$ requires curvature across strikes ⇒ smile/skew matters

### Key assumptions behind Dupire
- The underlying is a 1D diffusion (no jumps)
- The option surface is smooth enough to differentiate
- The surface is arbitrage-free (monotonicity + convexity in strike)

In practice, raw market quotes are noisy, so desks fit a smooth implied vol surface (e.g., SVI) before applying Dupire.

---

## C5. How interpolation works (why it’s needed)

Dupire gives $\sigma_{\text{loc}}(T,K)$ on a discrete grid.

But in simulation, we need $\sigma_{\text{loc}}(t,S_t)$ for arbitrary $(t,S_t)$.

So we build an interpolator:

$$
\sigma_{\text{loc}}(t,S) \approx \text{Interp}\big(\{(T_i,K_j,\sigma_{ij})\},\,t,\,S\big)
$$

Common choices:
- bilinear interpolation (simple, stable)
- bicubic/splines (smoother but can overshoot)

In MC, a common mapping is:
- evaluate Dupire surface at $K=S$ (or $K \approx S$) along the path.

Note that we do bilinear interpolation in our code

---

## C6. Practical note — why we assume an implied volatility surface

In a real production environment, constructing a local volatility model starts from **actual market option data**. Traders observe option prices or implied volatilities for many strikes and maturities:

$$
\sigma_{\text{imp}}(T,K)
$$

These quotes form a discrete surface, and the usual pipeline is:

$$
\text{Market Quotes}
\;\longrightarrow\;
\text{Call Prices}
\;\longrightarrow\;
\text{Smooth Surface}
\;\longrightarrow\;
\text{Dupire}
\;\longrightarrow\;
\sigma_{\text{loc}}(T,K)
$$

However, for this blog, collecting real option data and cleaning it introduces significant complexity:

- market quotes are noisy  
- arbitrage violations must be removed  
- surfaces must be smoothed  
- interpolation must be stable  

All of these steps are important in practice, but they distract from the main goal here: understanding the **mechanics of local volatility Monte Carlo pricing**.

---

### Synthetic implied volatility surface

Instead, we assume a smooth functional form for the implied volatility surface:

$$
\sigma_{\text{imp}}(K,T)
=
0.20
+
0.10\left(\frac{K}{F(T)}-1\right)^2
+
0.05\,\sqrt{T}
$$

where

$$
F(T)=S_0 e^{(r-q)T}
$$

is the forward price.

---

### Why this choice is realistic enough

This surface is deliberately simple but still captures the three key empirical features observed in real markets:

1. **Smile / curvature**

The quadratic term

$$
\left(\frac{K}{F(T)}-1\right)^2
$$

creates curvature across strikes.

---

2. **Term structure**

The linear maturity term

$$
0.05\,\sqrt{T}
$$

makes volatility vary with maturity.

---

3. **Forward moneyness scaling**

Using $K/F(T)$ instead of $K/S_0$ ensures the smile shape is consistent across maturities and not distorted by rates or dividends.

---

### Why this is useful pedagogically

Using a synthetic surface allows us to focus on the essential modeling pipeline:

$$
\text{Implied Vol Surface}
\;\rightarrow\;
\text{Call Prices}
\;\rightarrow\;
\text{Dupire Local Vol}
\;\rightarrow\;
\text{Monte Carlo Pricing}
$$

This keeps the implementation:

- clean
- reproducible
- stable numerically
- easy to understand

while still remaining realistic.

---

## C6. Pricing under local vol: naive vs bridge

### Simulation step (log-Euler)
Because $\sigma_{\text{loc}}(t,S)$ changes with $S$, exact transition is not known. We use log-Euler:

$$
\ln S_{t+\Delta t}
=
\ln S_t
+
\left(r-q-\frac12\sigma^2\right)\Delta t
+
\sigma\sqrt{\Delta t}Z
$$

where $\sigma = \sigma_{\text{loc}}(t,S_t)$ is frozen over the step.

### Barrier monitoring
- **Naive**: check barrier at time grid only
- **Bridge approximation**: freeze $\sigma$ within step and reuse GBM bridge probability

$$
p_{\text{hit}} \approx
\exp\left(
-\frac{2(x-b)(y-b)}{\sigma^2\Delta t}
\right)
$$

Assumption:
- local vol doesn’t change too much within a step (better for small $\Delta t$)

---

# D. Stochastic Volatility (Heston)

## D1. Why do we need stochastic volatility?

Local vol fits the implied surface today, but implies deterministic volatility dynamics. Real markets show:

- volatility clustering
- random variance
- leverage effect (returns and variance are negatively correlated)

Stochastic volatility models address these by modeling variance as its own random process.

---

## D2. Heston dynamics and parameter intuition

Heston under $Q$:

$$
\begin{aligned}
dS_t &= (r-q)S_t\,dt + \sqrt{v_t}\,S_t\,dW_1 \\
dv_t &= \kappa(\theta - v_t)\,dt + \xi\sqrt{v_t}\,dW_2
\end{aligned}
$$

Correlation:

$$
dW_1\,dW_2 = \rho\,dt
$$

Parameter meanings:

- $v_t$: instantaneous variance
- $\theta$: long-run variance (mean level)
- $\kappa$: mean reversion speed (how quickly variance reverts to $\theta$)
- $\xi$: vol of vol (variance randomness)
- $\rho$: correlation (negative in equities → skew)

---

## D3. Discretization (full truncation Euler)

Variance must remain nonnegative. A common scheme:

1) truncate variance for diffusion term:

$$
v_t^+ = \max(v_t,0)
$$

2) update:

$$
v_{t+\Delta t} =
v_t + \kappa(\theta-v_t^+)\Delta t + \xi\sqrt{v_t^+}\sqrt{\Delta t}\,Z_2
$$

then clip:

$$
v_{t+\Delta t} \leftarrow \max(v_{t+\Delta t},0)
$$

Stock update (log-Euler):

$$
\ln S_{t+\Delta t} =
\ln S_t + \left(r-q-\frac12v_t^+\right)\Delta t + \sqrt{v_t^+}\sqrt{\Delta t}\,Z_1
$$

with correlated normals:

$$
Z_2 = \rho Z_1 + \sqrt{1-\rho^2}\,\tilde Z
$$

---

## D4. Pricing naive and with bridge

### Naive
Check barrier only at endpoints.

### Bridge approximation under Heston
Exact bridge probability is unknown because volatility varies inside the step. A practical approximation is:

- freeze variance $v_t$ during the step
- treat it like constant vol with $\sigma^2=v_t$

$$
p_{\text{hit}} \approx \exp\left(-\frac{2(x-b)(y-b)}{v_t\,\Delta t}\right)
$$

Assumption:
- variance doesn’t move too much over one step (good when $\Delta t$ is small)

---

# E. Jump Diffusion (Merton)

## E1. Why do we need jumps?

Diffusions imply continuous paths. Markets exhibit:

- earnings gaps
- flash crashes
- macro shocks
- heavy tails

Jump diffusion adds discontinuities, improving tail behavior and skew.

---

## E2. Merton jump diffusion equation

Under $Q$:

$$
\frac{dS_t}{S_{t^-}}
=
(r-q-\lambda\kappa_J)\,dt
+
\sigma\,dW_t
+
(J-1)\,dN_t
$$

Where:
- $N_t$ is Poisson with intensity $\lambda$ (jumps per year)
- $J$ is the jump multiplier (random)
- $S_{t^-}$ is the pre-jump price

Interpretation:
- most of the time, $S_t$ evolves like GBM
- occasionally, a jump occurs and multiplies the price by $J$

---

## E3. How jump process is modeled (J, J−1, drift adjustment)

### Jump multiplier
A standard choice:

$$
J = e^Y,\quad Y\sim N(\mu_J,\sigma_J^2)
$$

So jump size is lognormal, ensuring $J>0$.

### Why $J-1$?
If a jump occurs:

$$
S_t = S_{t^-}\cdot J
$$

so relative jump return is:

$$
\frac{S_t - S_{t^-}}{S_{t^-}} = J-1
$$

That is why the SDE uses $(J-1)dN_t$.

### Drift adjustment (martingale condition)
Jumps add expected growth:

$$
\kappa_J = \mathbb{E}[J-1]
$$

For lognormal jumps:

$$
\kappa_J = e^{\mu_J+\frac12\sigma_J^2}-1
$$

Risk-neutral pricing requires:

$$
\mathbb{E}^Q[S_T] = S_0e^{(r-q)T}
$$

So we subtract $\lambda\kappa_J$ from drift:

$$
(r-q-\lambda\kappa_J)
$$

---

## E4. Simulation and naive DKI put pricing

Over one step, log-price evolves as:

$$
\ln S_{t+\Delta t} =
\ln S_t
+
\left(r-q-\lambda\kappa_J-\frac12\sigma^2\right)\Delta t
+
\sigma\sqrt{\Delta t}Z
+
\sum_{k=1}^{N}Y_k
$$

where:
- $N\sim \text{Poisson}(\lambda\Delta t)$
- $Y_k\sim N(\mu_J,\sigma_J^2)$

And efficiently:

$$
\sum_{k=1}^N Y_k \sim N(N\mu_J,\;N\sigma_J^2)
$$

Naive barrier pricing checks:

$$
\min_i S_{t_i}\le B
$$

Limitation:
- barrier might be crossed inside step due to diffusion
- barrier might be crossed instantly due to jumps

A better implementation combines:
- Brownian bridge when $N=0$
- substepping or jump-time simulation when $N>0$

---

# Summary (what matters in practice)

1. **Barrier monitoring is the biggest source of bias** in naive MC.
2. Brownian bridge corrections give “almost continuous” monitoring for diffusion models.
3. Local vol matches today’s surface; Heston adds realistic variance dynamics; Merton adds tail jumps.
4. The modeling choice changes barrier hit probabilities significantly, which is why these models matter for exotics.

---
