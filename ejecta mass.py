#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar 24 15:22:41 2023

@author: neilcao
"""

import numpy as np
import random
import matplotlib.pyplot as plt


# Constants(?)
G = 6.6743e-20 # km^3 kg^-1 s^-2, gravitational constant
c = 299792.458 # km s^-1, speed of light
M_sun = 1.989e30 # kg, mass of the sun


# Fitted coefficients
a=-8.1324
b=1.4820
d=1.7784

e=-9.3335
f=114.17
g=-337.56
n=1.5465


# Input parameters
R1 = 12 # km, radius of primary component
R2 = 12 # km, radius of secondary component
alpha_dyn = 1 # Coefficient of uncertainty
f_loss = 0.5 # Fraction of disk mass ejected

def calc_compactness(M, R, G, c):
    C = G * M * M_sun / (R * c**2)
    return C


def calc_disk_mass(M1, M2, R1, R2, G, c, a, b, d):
    if M1 < M2:
        M_min = M1
        C_min = calc_compactness(M1, R1, G, c)
    else:
        M_min = M2
        C_min = calc_compactness(M2, R2, G, c)
            
            
    if a * C_min + b < 0:
        m_disk = M_min * 5e-4
    else:
        m_disk = M_min * max(5e-4, (a * C_min + b)**d)

    return m_disk


def calc_dynamical_mass(M1, M2, R1, R2, G, c, e, f, g, n):
    C1 = calc_compactness(M1, R1, G, c)
    C2 = calc_compactness(M2, R2, G, c)

    m_dyn = ((e / C1 + f * (M2**n) / (M1**n) + g * C1) * M1 + (e / C2 + f * (M2**n) / (M2**n) + g * C2) * M2) * 1e-3
    if m_dyn < 0:
        m_dyn = 0
    return m_dyn


def calc_single_merger_ejecta(M1, M2, R1, R2, G, c, alpha_dyn, f_loss):
    m_disk = calc_disk_mass(M1, M2, R1, R2, G, c, a, b, d)
    m_dyn = calc_dynamical_mass(M1, M2, R1, R2, G, c, e, f, g, n)

    m_ej = alpha_dyn * m_dyn + f_loss * m_disk
    
    return m_ej


def calc_total_ejecta(M1, M2, R1, R2, G, c, alpha_dyn, f_loss, r_merger):
    m_ej = calc_single_merger_ejecta(M1, M2, R1, R2, G, c, alpha_dyn, f_loss)

    m_total = m_ej * r_merger
    return m_total


a1, b1 = 1.1, 2.0
N = 10000
mu1, sigma1 = 440, 500


samples1 = np.random.normal(mu1, sigma1, N)

def sample_top_hat(a1, b1, N):
    samples = [random.uniform(a1, b1) for _ in range(N)]
    return samples


samples2 = sample_top_hat(a1, b1, N)
samples3 = sample_top_hat(a1, b1, N)

accepted_single_ejecta_samples = []
accepted_total_ejecta_samples = []


for M1, M2 in zip(samples2, samples3):
    m_ej = calc_single_merger_ejecta(M1, M2, R1, R2, G, c, alpha_dyn, f_loss)

    accepted_single_ejecta_samples.append(m_ej)


max_m_total = 0
max_M1 = 0
max_M2 = 0
max_r_merger = 0
accepted_total_ejecta_samples = []

for r_merger in samples1:
    for i, m_ej in enumerate(accepted_single_ejecta_samples):
        m_total = r_merger * m_ej
        if m_total > 0:
            accepted_total_ejecta_samples.append(m_total)
            if m_total > max_m_total:
                max_m_total = m_total
                max_M1 = samples2[i]
                max_M2 = samples3[i]
                max_r_merger = r_merger


C1 = calc_compactness(M1, R1, G, c)
C2 = calc_compactness(M2, R2, G, c)
m_disk = calc_disk_mass(M1, M2, R1, R2, G, c, a, b, d)
m_dyn = calc_dynamical_mass(M1, M2, R1, R2, G, c, e, f, g, n)
m_ej = calc_single_merger_ejecta(M1, M2, R1, R2, G, c, alpha_dyn, f_loss)


print("Total Ejecta Mass Mean:", np.mean(accepted_total_ejecta_samples))
print("Total Ejecta Mass Uncertainty:", np.std(accepted_total_ejecta_samples))

print("Maximum Total Ejecta Mass:", max_m_total)
print("M1 for Maximum Total Ejecta Mass:", max_M1)
print("M2 for Maximum Total Ejecta Mass:", max_M2)
print("r_merger for Maximum Total Ejecta Mass:", max_r_merger)

# Plotting the results
plt.hist(accepted_single_ejecta_samples, bins=100, density=True)
plt.xlabel('Total Ejecta Mass')
plt.ylabel('Density')
plt.title('Distribution of BNS Ejecta Mass')
plt.show()

plt.hist(accepted_total_ejecta_samples, bins=100, density=True)
plt.xlabel('Total Ejecta Mass')
plt.ylabel('Density')
plt.title('Distribution of BNS Ejecta Mass')
plt.show()
