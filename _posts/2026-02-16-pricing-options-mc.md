---
title: 'Monte Carlo Pricing of Barrier Options'
date: 2026-02-16
permalink: /posts/2026/02/pricing-mc-options/
tags:
  - places
---

# Monte Carlo Pricing of Barrier Options  
### From Constant Volatility to Local Volatility, Stochastic Volatility, and Jumps

This write-up explains **what we are doing** and **why we are doing it** when pricing barrier options using Monte Carlo (MC). My goal is a technical but accescible article.

---

# A. Foundations

## A1. Why Monte Carlo pricing? When do analytical models fail?

In risk-neutral pricing, the value of a derivative is the discounted expectation of its payoff under the **risk-neutral measure** \(Q\):

$V_0 = e^{-rT}\mathbb{E}^Q[\text{Payoff}]$



Analytical (closed-form) solutions exist only when the model and payoff are simple enough that we can evaluate the expectation exactly. Examples:

- **Black–Scholes** gives closed forms for European calls/puts under constant volatility.
- **Heston** has semi-closed forms for European options via Fourier methods.
- **PDE** methods can price many products in low dimensions (typically 1–2 state variables).

Analytical methods start to fail (or become un-reliable/too complex to model) when:

### 1) Payoffs are path-dependent
Barrier options depend on the entire path:

- Down-and-in put payoff depends on whether \(S_t\) ever crossed a barrier.

\[
\mathbf{1}\left\{\min_{0\le t\le T} S_t \le B\right\}
\]

This indicator makes the pricing expectation more complex than standard European payoff \(f(S_T)\).

### 2) The model contains additional state variables
Examples:
- stochastic volatility adds \(v_t\)
- local volatility makes \(\sigma\) a function of \((t,S)\)
- jumps make the process discontinuous

Each added feature increases the complexity of either:
- solving a PDE, or
- finding a closed-form transform