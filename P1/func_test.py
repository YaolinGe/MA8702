import numpy as np
import random
import scipy.stats as st
import matplotlib.pyplot as plt


def normal(x, mu, sigma):
    numerator = np.exp(-1 * ((x - mu) ** 2) / (2 * sigma ** 2))
    denominator = sigma * np.sqrt(2 * np.pi)
    return numerator / denominator


def neg_log_prob(x, mu, sigma):
    return -1 * np.log(normal(x=x, mu=mu, sigma=sigma))


def HMC(mu=0.0, sigma=1.0, path_len=1.0, step_size=0.25, initial_position=0.0, epochs=1000):
    # setup
    steps = int(path_len / step_size)  # path_len and step_size are tricky parameters to tune...
    samples = [initial_position]
    momentum_dist = st.norm(0, 1)
    # generate samples
    for e in range(epochs):
        q0 = np.copy(samples[-1])
        q1 = np.copy(q0)
        p0 = momentum_dist.rvs()
        p1 = np.copy(p0)
        dVdQ = -1 * (q0 - mu) / (sigma ** 2)  # gradient of PDF wrt position (q0) aka potential energy wrt position

        # leapfrog integration begin
        for s in range(steps):
            p1 += step_size * dVdQ / 2  # as potential energy increases, kinetic energy decreases, half-step
            q1 += step_size * p1  # position increases as function of momentum
            p1 += step_size * dVdQ / 2  # second half-step "leapfrog" update to momentum
        # leapfrog integration end
        p1 = -1 * p1  # flip momentum for reversibility

        # metropolis acceptance
        q0_nlp = neg_log_prob(x=q0, mu=mu, sigma=sigma)
        q1_nlp = neg_log_prob(x=q1, mu=mu, sigma=sigma)

        p0_nlp = neg_log_prob(x=p0, mu=0, sigma=1)
        p1_nlp = neg_log_prob(x=p1, mu=0, sigma=1)

        # Account for negatives AND log(probabiltiies)...
        target = q0_nlp - q1_nlp  # P(q1)/P(q0)
        adjustment = p1_nlp - p0_nlp  # P(p1)/P(p0)
        acceptance = target + adjustment

        event = np.log(random.uniform(0, 1))
        if event <= acceptance:
            samples.append(q1)
        else:
            samples.append(q0)

    return samples

mu = 0
sigma = 1
trial = HMC(mu=mu, sigma=sigma, path_len=1.5, step_size=0.025)

lines = np.linspace(-6, 6, 10000)
normal_curve = [normal(x=l, mu=mu, sigma=sigma) for l in lines]

plt.plot(lines, normal_curve)
plt.hist(trial, density=True, bins=20)
plt.show()

#%%
import numpy as np
import matplotlib.pyplot as plt

def hmc(U, dU, x0 = None, D = None, EP = 1000, L = 5, delta = 0.05):
    burnin = int(L / 2)

    D = x0.size

    K = lambda p: np.dot(p, p) / 2
    dK = lambda p: p
    xStar = np.random.randn(D)
    pStar = np.random.randn(D)

    x = [xStar]
    p = [pStar]

    for i in range(EP):
        x1 = x[-1]
        U_x1 = U(x1)
        p1 = np.random.randn(D)
        K_p1 = K(p1)

        p2 = np.random.randn(D)

        K_p2 = K(p2)
        H0 = U_x1 + K_p2

        pStar = p2
        xStar = x1

        for j in range(L):
            xStar += delta * dK(pStar)
            pStar -= delta * dU(xStar)

        x2 = xStar
        p3 = pStar
        U_x2 = U(x2)
        K_p3 = K(p3)
        Hstar = U_x2 + K_p3








