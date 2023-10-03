from scipy.stats import gamma

# parameters
alpha = 10
beta = 1 / 0.06
x = 150

# compute the cumulative distribution function (CDF)
gamma_cdf = 1 - gamma.cdf(x, alpha, scale=beta)

print(gamma_cdf)
