import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.ticker
from mpl_toolkits import mplot3d
from astropy import units as u
import pdb

#pdb.set_trace()

'''
Equations and appendicies are adapted from https://ui.adsabs.harvard.edu/abs/1996ApOpt..35.1566C/

'''

# wavelengths in Angstrom as numpy array.
# In the source code, the wavelength comes from an ndarray called wave. It is 
# convereted to Angstrom using astropy.units:
#
# wave = wave.to(u.AA)
# wavelegnth = wave.value
wavelength = np.asarray([3000, 6000, 9000])

sigma = (1.e4/wavelength)
sigma_sq = sigma**2


## Variables to Change ##

# fractional humidity (%)
h = .199

# temperature [C]
t = 0
# temperature [K]
T = t + 273.15

# pressure [mbar]
p = 61783

# ppm of CO2
x_c = 450


### Refractive Index for standard air

# Constants for equation (1) from Appendix A
k0 = 238.0185
k1 = 5792105
k2 = 57.362
k3 = 167917

# The equation used for index of refraction for standard air is:
# 10^8(n_as - 1) = k1/(k0 - sigma_sq) + k3/(k2 - sigma_sq)

# Index of refraction for standard air. From equation (1)
n_as = ((k1 / (k0 - sigma_sq) + k3 / (k2 - sigma_sq)) / 10**8) + 1

# Index of refraction for standard air with x_c ppm CO2. From equation (2)
n_axs = ((n_as - 1) * (1 + 0.534 * 10**-6 * (x_c - 450))) + 1


### Refractive index for pure water vapor ###

# constants for equation
cf = 1.022 #cf is the correction factor
w0 = 295.235
w1 = 2.6422
w2 = -0.032380
w3 = 0.004028

# The equation used to determine the index of refraction for standard water vapor 
# (20 deg C, 1333 Pa) is:
#
# 10^8(n_ws - 1) = cf(w0 + w1(sigma^2) + w2(sigma^4) + w3(sigma^6))
# Adapted from equation (3)

n_ws_step1 = cf * (w0 + (w1 * sigma**2) + (w2 * sigma**4) + (w3 * sigma**6))
n_ws = n_ws_step1 / 10**8 + 1


def M_a_func(x_c):
    """

    Returns the molar mass of dry air with x_c ppm of CO2. Adapted from Section 3
    below Equation (4)

    Args:
        x_c (int): ppm of CO2 

    Returns:
        Molar mass of dry air with x_c ppm of CO2

    """
    
    return 10**-3 * (28.9635 + 12.011 * 10**-6 * (x_c - 400))

def x_w_func(T, p, h):
    """

    Returns the molar fraction of water vapor in moist air. Adapted from Section 3
    below Equation (4)

    Args:
        T (float): temperature in Kelvin [K]
        p (float): pressure in Pascals [Pa]
        h (float): fractional humidity 

    Returns:
        Molar fraction of water vapor in moist air

    """

    # temperature [C]
    t = T - 273.15

    # constants for f from Appendix A
    alpha = 1.00062
    beta = 3.14 * 10**-8
    gamma  = 5.6 * 10**-7

    # enchancement factor
    # found in Section 3 below equation (4)
    f = alpha + beta * p + gamma * t**2

    # constants for svp from Appendix A
    A = 1.2378847 * 10**-5
    B = -1.9121316 * 10**-2
    C = 33.93711047
    D = -6.3431645 * 10**3 

    # saturation vapor pressure (svp) of water vapor in air at temperature, T [Pa]
    # found in Section 3 below equation (4)
    svp = np.exp(A * T**2 + B * T + C + (D/T)) 

    # molar fraction of water vapor in moist air
    # found in Section 3 below equation (4)
    x_w = f * h * svp / p

    return x_w 
 
def Z_func(T, p, x_w):
    """

    Returns the compressability of moist air. Adapted from equation (12) in Appendix A.

    Args:
        T (float): temperature in Kelvin [K]
        p (float): pressure in Pascals [Pa]
        x_w (float): molar fraction of water vapor in moist air

    Returns:
        Compresability of the moist air

    """

    # t = temperature [C]
    t = T - 273.15

    # constants for Z (from Appendix A)
    a0 = 1.58123 * 10**-6
    a1 = -2.9331 * 10**-8
    a2 = 1.1043 * 10**-10
    b0 = 5.707 * 10**-6
    b1 = -2.051 * 10**-8
    c0 = 1.9898 * 10**-4
    c1 = -2.376 * 10**-6
    d = 1.83 * 10**-11
    e = -0.765 * 10**-8

    # compresibility of moist air
    # adapted from equation (12) in Appendix A
    Z = 1 - (p/T) * (a0 + a1*t + a2 * t**2 + (b0 + b1*t) * x_w + (c0 + c1*t) * x_w**2) + (p/T)**2 * (d + e * x_w**2)
    return Z

def density_func(T, p, h, x_c, x_w=None):
    """
    Return the density of moist air given. Adapted from equation (4)

    Args:
        T (float): temperature in Kelvin [K]
        p (float): pressure in Pascals [Pa]
        h (float): fractional humidity 
        x_c (int): ppm of CO2
        x_w (int): optional parameter for molar fraction of water vapor in moist air

    Returns:
        Density of the moist air

    """

    # molar mass of water vapor [kg/mol]
    M_w = 0.018015

    # molar mass of dry air with x_c ppm of CO2 [kg/mol]
    M_a = M_a_func(x_c)

    # gas constant [J/(mol*K)]
    R = 8.314510

    # molar fraction of water vapor in moist air
    if x_w is None:
        x_w = x_w_func(T, p, h)

    # compresibility of moist air
    Z = Z_func(T, p, x_w)

    # density of moist air
    # adapted from equation (4)
    rho = (p * M_a / Z * R * T) * (1 - x_w * (1 - M_w / M_a))
    
    return rho

#### Refractive Index of Moist Air ####

def n_prop_func(T, p, h, x_c, x_w=None):
    # molar mass of water vapor [kg/mol]
    M_w = 0.018015
    # molar mass of dry air with x_c ppm of CO2 [kg/mol]
    M_a = 10**-3 * (28.9635 + 12.011 * 10**-6 * (x_c - 400))
    # gas constant [J/(mol*K)]
    R = 8.314510

    # density of standard dry air (15 deg C, 101325 Pa, x_w = 0, 450 ppm CO2)
    rho_axs = density_func(288.15, 101325, 0, 450, 0)

    # density of standard water vapor (20 deg C, 1333 Pa, x_w = 1)
    rho_ws = density_func(293.15, 1333, 1, 0, 1)

    # molar fraction of water vapor in moist air
    if x_w is None:
        x_w = x_w_func(T, p, h)

    # compressability of moist air
    Z = Z_func(T, p, x_w)

    # density of dry air component of moist air with actual conditions
    # From step 8. in Appendix B
    rho_a = p * M_a * (1 - x_w) / Z * R * T
    # density of water vapor component of moist air with actual conditons 
    # From step 9. in Appendix B
    rho_w = p * M_w * x_w / Z * R * T

    # refractive index of moist air
    # From Equation (5)
    n_prop = 1 + ((rho_a/rho_axs) * (n_axs - 1) + (rho_w/rho_ws) * (n_ws - 1))

    return n_prop

# refractive index of moist air
# adapted from Equation (5)
n_prop = n_prop_func(T, p, h, x_c)

# only modify above 2000A
n_prop = n_prop*(wavelength >= 2000.) + 1.*(wavelength < 2000.)

# assign variables for the different methods
corrected_wavelengths = wavelength*n_prop

print(f'Refractive index for pressure at a very small number (can not be 0 because p is the divisor in the x_w equation): {n_prop_func(275.25, 0.00000001, .159, 429)}')

print(f'Refractive index for ppm of CO2 at 1 million (0 deg C, 1 atm, 0% humidity): {n_prop_func(273.15, 101325, 0, 1000000)}')
# at these same numbers, the index of refraction is 1.00044856, kinda seems like a meaningful difference

# h = molar fraction of water in air = xw = nw / (na + nw)
# na = mols of air
# nw = mols of water

'''
# test prints

print("Old refractive indicies: ", n_axs)
print("New refractive indicies: ", n_prop, "\n")

print("Old corrected wavelengths: ", old_corrected_wavelengths)
print("New corrected wavelengths: ", corrected_wavelengths)

v_old = ((old_corrected_wavelengths - wavelength)/ wavelength) * 3.0 * 10e8
v = ((corrected_wavelengths - wavelength)/ wavelength) * 3.0 * 10e8

print("v_old: ", v_old)
print("v_new: ", v)
print(v_old-v)
'''

# line styles for different wavelengths
styles = ['-', (0, (8, 3)), (0, (6, 3, 2, 3))]

# color labels
color_labels = ['Old Refractive Index', 'New Refractive Index', 'Delta Wavelength']


'''Varying Temperature Graph at Standard Air Conditions'''


# temp range
temp_range = np.arange(253, 298, .1)

# get refractive index and change in wavelength 
n_temp1 = []
n_temp2 = []
n_temp3 = []

delta_temp1 = []
delta_temp2 = []
delta_temp3 = []
for temperature in temp_range:
    # standard dry air has T=288.15[K], 101325 [Pa], 0% humidity, 450 ppm CO2
    n_prop = n_prop_func(temperature, 101325, 0, 450)
    n_temp1.append(n_prop[0])
    n_temp2.append(n_prop[1])
    n_temp3.append(n_prop[2])

    delta_wavelength = np.abs(wavelength * n_prop - wavelength)
    delta_temp1.append(delta_wavelength[0])
    delta_temp2.append(delta_wavelength[1])
    delta_temp3.append(delta_wavelength[2])

# setup plot
fig, ax1 = plt.subplots()
plt.ticklabel_format(axis='y', style='plain', useOffset=False)
plt.title("Varying Temperature [K] at Standard Dry Air Conditions")

# graph the original n vals
n_axs_temp = np.zeros(len(temp_range)) + n_axs[0]
ax1.plot(temp_range, n_axs_temp, color='g',label=wavelength[0], ls=styles[0], lw=.8)

n_axs_temp = np.zeros(len(temp_range)) + n_axs[1]
ax1.plot(temp_range, n_axs_temp, color='g', label=wavelength[1], ls=styles[1], lw=.8)

n_axs_temp = np.zeros(len(temp_range)) + n_axs[2]
ax1.plot(temp_range, n_axs_temp, color='g', label=wavelength[2], ls=styles[2], lw=.8)

# graph the Refracive Index vs Temperature for each wavelength
color = 'tab:orange'
ax1.set_xlabel('Temperature [K]')
ax1.set_ylabel('Refractive Index', color=color)
ax1.plot(temp_range, n_temp1, color=color, label=wavelength[0], ls=styles[0], lw=.8)
ax1.plot(temp_range, n_temp2, color=color, label=wavelength[1], ls=styles[1], lw=.8)
ax1.plot(temp_range, n_temp3, color=color, label=wavelength[2], ls=styles[2], lw=.8)
ax1.tick_params(axis='y', labelcolor=color)
ax1.set_ylim(top=1.00031)

ax2 = ax1.twinx()

# graph the change in wavelength vs Temperature for each wavelength
color = 'tab:blue'
ax2.set_ylabel('Change in Wavelength [Å]', color=color)
ax2.plot(temp_range, delta_temp1, color=color, label=f'{wavelength[0]} Å', ls=styles[0], lw=.8)
ax2.plot(temp_range, delta_temp2, color=color, label=f'{wavelength[1]} Å', ls=styles[1], lw=.8)
ax2.plot(temp_range, delta_temp3, color=color, label=f'{wavelength[2]} Å', ls=styles[2], lw=.8)
ax2.tick_params(axis='y', labelcolor=color)
ax2.set_ylim(top=2.85)

# Get handles and labels for all data sets
handles1, _ = ax1.get_legend_handles_labels()
handles2, labels2 = ax2.get_legend_handles_labels()

# get the solid color lines from the 1st and 4th handles of ax1 and the 1st handle of ax2
color_handles = [handles1[0], handles1[3], handles2[0]]

# make new lines for linestyle legend
ls_handle1 = plt.Line2D([0], [0], ls='-', color='black', markersize=10)
ls_handle2 = plt.Line2D([0], [0], ls='--', color='black', markersize=10)
ls_handle3 = plt.Line2D([0], [0], ls='-.', color='black', markersize=10)

# setup legends 
color_legend = ax1.legend(handles=color_handles, labels=color_labels, loc='upper center')

linestyle_legend = ax2.legend([ls_handle1, ls_handle2, ls_handle3], labels2, loc='upper left')

ax1.add_artist(color_legend)
ax2.add_artist(linestyle_legend)

# save plot as png
fig.tight_layout() 
plt.savefig('temp_standard_air.png')
plt.close()


'''Varying Temperature Graph at Standard Air Conditions vs Old'''


# temp range
temp_range = np.arange(253, 298, .1)

# get refractive index and change in wavelength 
n_temp1 = []
n_temp2 = []
n_temp3 = []

delta_temp1 = []
delta_temp2 = []
delta_temp3 = []
for temperature in temp_range:
    # standard dry air has T=288.15[K], 101325 [Pa], 0% humidity, 450 ppm CO2
    n_prop = n_prop_func(temperature, 101325, 0, 450)
    n_temp1.append(n_prop[0])
    n_temp2.append(n_prop[1])
    n_temp3.append(n_prop[2])

    delta_wavelength = np.abs(wavelength * n_prop - wavelength * n_axs)
    delta_temp1.append(delta_wavelength[0])
    delta_temp2.append(delta_wavelength[1])
    delta_temp3.append(delta_wavelength[2])

# setup plot
fig, ax1 = plt.subplots()
plt.ticklabel_format(axis='y', style='plain', useOffset=False)
plt.title("New vs Old Refractive Index w/ Varying Temperature at Std. Air")

# graph the original n vals
n_axs_temp = np.zeros(len(temp_range)) + n_axs[0]
ax1.plot(temp_range, n_axs_temp, color='g',label=wavelength[0], ls=styles[0], lw=.8)

n_axs_temp = np.zeros(len(temp_range)) + n_axs[1]
ax1.plot(temp_range, n_axs_temp, color='g', label=wavelength[1], ls=styles[1], lw=.8)

n_axs_temp = np.zeros(len(temp_range)) + n_axs[2]
ax1.plot(temp_range, n_axs_temp, color='g', label=wavelength[2], ls=styles[2], lw=.8)

# graph the Refracive Index vs Temperature for each wavelength
color = 'tab:orange'
ax1.set_xlabel('Temperature [K]')
ax1.set_ylabel('Refractive Index', color=color)
ax1.plot(temp_range, n_temp1, color=color, label=wavelength[0], ls=styles[0], lw=.8)
ax1.plot(temp_range, n_temp2, color=color, label=wavelength[1], ls=styles[1], lw=.8)
ax1.plot(temp_range, n_temp3, color=color, label=wavelength[2], ls=styles[2], lw=.8)
ax1.tick_params(axis='y', labelcolor=color)
ax1.set_ylim(top=1.00031)

ax2 = ax1.twinx()

# graph the change in wavelength vs Temperature for each wavelength
color = 'tab:blue'
ax2.set_ylabel('Old and New Model Wavelength Adjustments Differnece [Å]', fontsize=9, color=color)
ax2.plot(temp_range, delta_temp1, color=color, label=f'{wavelength[0]} Å', ls=styles[0], lw=.8)
ax2.plot(temp_range, delta_temp2, color=color, label=f'{wavelength[1]} Å', ls=styles[1], lw=.8)
ax2.plot(temp_range, delta_temp3, color=color, label=f'{wavelength[2]} Å', ls=styles[2], lw=.8)
ax2.tick_params(axis='y', labelcolor=color)

# Get handles and labels for all data sets
handles1, _ = ax1.get_legend_handles_labels()
handles2, labels2 = ax2.get_legend_handles_labels()

# get the solid color lines from the 1st and 4th handles of ax1 and the 1st handle of ax2
color_handles = [handles1[0], handles1[3], handles2[0]]

# make new lines for linestyle legend
ls_handle1 = plt.Line2D([0], [0], ls='-', color='black', markersize=10)
ls_handle2 = plt.Line2D([0], [0], ls='--', color='black', markersize=10)
ls_handle3 = plt.Line2D([0], [0], ls='-.', color='black', markersize=10)

# setup legends 
color_legend = ax1.legend(handles=color_handles, labels=color_labels, loc='upper center')

linestyle_legend = ax2.legend([ls_handle1, ls_handle2, ls_handle3], labels2, loc='upper left')

ax1.add_artist(color_legend)
ax2.add_artist(linestyle_legend)

# save plot as png
fig.tight_layout() 
plt.savefig('temp_standard_air_vs_old.png')
plt.close()


'''Varying Temperature Graph at Summit Air Conditions'''


# get refractive index and change in wavelength 
n_temp1 = []
n_temp2 = []
n_temp3 = []

delta_temp1 = []
delta_temp2 = []
delta_temp3 = []
for temperature in temp_range:
    # temperature, pressure, humidity from KR.20240505.22262.90_icubes.fits
    n_prop = n_prop_func(temperature, 61785, .154, 426.9)
    n_temp1.append(n_prop[0])
    n_temp2.append(n_prop[1])
    n_temp3.append(n_prop[2])

    delta_wavelength = np.abs(wavelength * n_prop - wavelength)
    delta_temp1.append(delta_wavelength[0])
    delta_temp2.append(delta_wavelength[1])
    delta_temp3.append(delta_wavelength[2])

# print slopes
slopes = []
slopes.append(np.polyfit(temp_range, delta_temp1, 1)[0])
slopes.append(np.polyfit(temp_range, delta_temp2, 1)[0])
slopes.append(np.polyfit(temp_range, delta_temp3, 1)[0])
print(f'Slopes of delta_temp/temp_range: {slopes}')

# setup plot
fig, ax1 = plt.subplots()
plt.ticklabel_format(axis='y', style='plain', useOffset=False)
plt.title("Varying Temperature [K] at Summit Air Conditions")

# plot average condition range on summit
plt.axvline(x=265.15, color='darkgrey')
plt.axvline(x=285.15, color='darkgrey')

# graph the original n vals
n_axs_temp = np.zeros(len(temp_range)) + n_axs[0]
ax1.plot(temp_range, n_axs_temp, color='g',label=wavelength[0], ls=styles[0], lw=.8)

n_axs_temp = np.zeros(len(temp_range)) + n_axs[1]
ax1.plot(temp_range, n_axs_temp, color='g', label=wavelength[1], ls=styles[1], lw=.8)

n_axs_temp = np.zeros(len(temp_range)) + n_axs[2]
ax1.plot(temp_range, n_axs_temp, color='g', label=wavelength[2], ls=styles[2], lw=.8)

# graph the Refracive Index vs Temperature for each wavelength
color = 'tab:orange'
ax1.set_xlabel('Temperature [K]')
ax1.set_ylabel('Refractive Index', color=color)
ax1.plot(temp_range, n_temp1, color=color, label=wavelength[0], ls=styles[0], lw=.8)
ax1.plot(temp_range, n_temp2, color=color, label=wavelength[1], ls=styles[1], lw=.8)
ax1.plot(temp_range, n_temp3, color=color, label=wavelength[2], ls=styles[2], lw=.8)
ax1.tick_params(axis='y', labelcolor=color)
ax1.set_ylim(top=1.00034)

ax2 = ax1.twinx()

# graph the change in wavelength vs Temperature for each wavelength
color = 'tab:blue'
ax2.set_ylabel('Change in Wavelength [Å]', color=color)
ax2.plot(temp_range, delta_temp1, color=color, label=f'{wavelength[0]} Å', ls=styles[0], lw=.8)
ax2.plot(temp_range, delta_temp2, color=color, label=f'{wavelength[1]} Å', ls=styles[1], lw=.8)
ax2.plot(temp_range, delta_temp3, color=color, label=f'{wavelength[2]} Å', ls=styles[2], lw=.8)
ax2.tick_params(axis='y', labelcolor=color)
ax2.set_ylim(top=1.8)

# Get handles and labels for all data sets
handles1, _ = ax1.get_legend_handles_labels()
handles2, labels2 = ax2.get_legend_handles_labels()

# get the solid color lines from the 1st and 4th handles of ax1 and the 1st handle of ax2
color_handles = [handles1[0], handles1[3], handles2[0]]

# make new lines for linestyle legend
ls_handle1 = plt.Line2D([0], [0], ls='-', color='black', markersize=10)
ls_handle2 = plt.Line2D([0], [0], ls='--', color='black', markersize=10)
ls_handle3 = plt.Line2D([0], [0], ls='-.', color='black', markersize=10)

# setup legends 
color_legend = ax1.legend(handles=color_handles, labels=color_labels, loc='upper center')

linestyle_legend = ax2.legend([ls_handle1, ls_handle2, ls_handle3], labels2, loc='upper left')

ax1.add_artist(color_legend)
ax2.add_artist(linestyle_legend)

# save plot as png
fig.tight_layout() 
plt.savefig('temp_summit_air.png')
plt.close()


'''Varying Temperature Graph at Summit Air Conditions vs Old Refractive Index'''


# get refractive index and change in wavelength 
n_temp1 = []
n_temp2 = []
n_temp3 = []

delta_temp1 = []
delta_temp2 = []
delta_temp3 = []
for temperature in temp_range:
    # temperature, pressure, humidity from KR.20240505.22262.90_icubes.fits
    n_prop = n_prop_func(temperature, 61785, .154, 426.9)
    n_temp1.append(n_prop[0])
    n_temp2.append(n_prop[1])
    n_temp3.append(n_prop[2])

    delta_wavelength = np.abs(wavelength * n_prop - wavelength * n_axs)
    delta_temp1.append(delta_wavelength[0])
    delta_temp2.append(delta_wavelength[1])
    delta_temp3.append(delta_wavelength[2])

# setup plot
fig, ax1 = plt.subplots()
plt.ticklabel_format(axis='y', style='plain', useOffset=False)
plt.title("New vs Old Refractive Index w/ Varying Temperature at Summit Air")

# plot average condition range on summit
plt.axvline(x=265.15, color='darkgrey')
plt.axvline(x=285.15, color='darkgrey')

# graph the original n vals
n_axs_temp = np.zeros(len(temp_range)) + n_axs[0]
ax1.plot(temp_range, n_axs_temp, color='g',label=wavelength[0], ls=styles[0], lw=.8)

n_axs_temp = np.zeros(len(temp_range)) + n_axs[1]
ax1.plot(temp_range, n_axs_temp, color='g', label=wavelength[1], ls=styles[1], lw=.8)

n_axs_temp = np.zeros(len(temp_range)) + n_axs[2]
ax1.plot(temp_range, n_axs_temp, color='g', label=wavelength[2], ls=styles[2], lw=.8)

# graph the Refracive Index vs Temperature for each wavelength
color = 'tab:orange'
ax1.set_xlabel('Temperature [K]')
ax1.set_ylabel('Refractive Index', color=color)
ax1.plot(temp_range, n_temp1, color=color, label=wavelength[0], ls=styles[0], lw=.8)
ax1.plot(temp_range, n_temp2, color=color, label=wavelength[1], ls=styles[1], lw=.8)
ax1.plot(temp_range, n_temp3, color=color, label=wavelength[2], ls=styles[2], lw=.8)
ax1.tick_params(axis='y', labelcolor=color)
ax1.set_ylim(top=1.00034)

ax2 = ax1.twinx()

# graph the change in wavelength vs Temperature for each wavelength
color = 'tab:blue'
ax2.set_ylabel('Old and New Model Wavelength Adjustments Differnece [Å]', fontsize=9, color=color)
ax2.plot(temp_range, delta_temp1, color=color, label=f'{wavelength[0]} Å', ls=styles[0], lw=.8)
ax2.plot(temp_range, delta_temp2, color=color, label=f'{wavelength[1]} Å', ls=styles[1], lw=.8)
ax2.plot(temp_range, delta_temp3, color=color, label=f'{wavelength[2]} Å', ls=styles[2], lw=.8)
ax2.tick_params(axis='y', labelcolor=color)
ax2.set_ylim(top=1.4)

# Get handles and labels for all data sets
handles1, _ = ax1.get_legend_handles_labels()
handles2, labels2 = ax2.get_legend_handles_labels()

# get the solid color lines from the 1st and 4th handles of ax1 and the 1st handle of ax2
color_handles = [handles1[0], handles1[3], handles2[0]]

# make new lines for linestyle legend
ls_handle1 = plt.Line2D([0], [0], ls='-', color='black', markersize=10)
ls_handle2 = plt.Line2D([0], [0], ls='--', color='black', markersize=10)
ls_handle3 = plt.Line2D([0], [0], ls='-.', color='black', markersize=10)

# setup legends 
color_legend = ax1.legend(handles=color_handles, labels=color_labels, loc='upper center')

linestyle_legend = ax2.legend([ls_handle1, ls_handle2, ls_handle3], labels2, loc='upper left')

ax1.add_artist(color_legend)
ax2.add_artist(linestyle_legend)

# save plot as png
fig.tight_layout() 
plt.savefig('temp_summit_air_vs_old.png')
plt.close()


""" Varying Pressure Graph at Standard Dry Air Conditions """

# pressure graph
pres_range = np.arange(40000, 130000, 100)

# get refractive index and change in wavelength 
n_pres1 = []
n_pres2 = []
n_pres3 = []

delta_pres1 = []
delta_pres2 = []
delta_pres3 = []
for pressure in pres_range:
    # standard dry air has T=288.15[K], 101325 [Pa], 0% humidity, 450 ppm CO2
    n_prop = n_prop_func(288.15, pressure, 0, 450)
    n_pres1.append(n_prop[0])
    n_pres2.append(n_prop[1])
    n_pres3.append(n_prop[2])

    delta_wavelength = np.abs(wavelength * n_prop - wavelength)
    delta_pres1.append(delta_wavelength[0])
    delta_pres2.append(delta_wavelength[1])
    delta_pres3.append(delta_wavelength[2])

# setup plot
fig, ax1 = plt.subplots()
plt.ticklabel_format(axis='y', style='plain', useOffset=False)
plt.title("Varying Pressure [Pa] at Standard Dry Air Conditions")

# graph the original n vals
n_axs_pres = np.zeros(len(pres_range)) + n_axs[0]
ax1.plot(pres_range, n_axs_pres, color='g',label=wavelength[0], ls=styles[0], lw=.8)

n_axs_pres = np.zeros(len(pres_range)) + n_axs[1]
ax1.plot(pres_range, n_axs_pres, color='g', label=wavelength[1], ls=styles[1], lw=.8)

n_axs_pres = np.zeros(len(pres_range)) + n_axs[2]
ax1.plot(pres_range, n_axs_pres, color='g', label=wavelength[2], ls=styles[2], lw=.8)

# graph the Refracive Index vs Pressure for each wavelength
color = 'tab:orange'
ax1.set_xlabel('Pressure [Pa]')
ax1.set_ylabel('Refractive Index', color=color)
ax1.plot(pres_range, n_pres1, color=color, label=wavelength[0], ls=styles[0], lw=.8)
ax1.plot(pres_range, n_pres2, color=color, label=wavelength[1], ls=styles[1], lw=.8)
ax1.plot(pres_range, n_pres3, color=color, label=wavelength[2], ls=styles[2], lw=.8)
ax1.tick_params(axis='y', labelcolor=color)
ax2 = ax1.twinx()

# graph the change in wavelength vs Pressure for each wavelength
color = 'tab:blue'
ax2.set_ylabel('Change in Wavelength [Å]', color=color)
ax2.plot(pres_range, delta_pres1, color=color, label=f'{wavelength[0]} Å', ls=styles[0], lw=.8)
ax2.plot(pres_range, delta_pres2, color=color, label=f'{wavelength[1]} Å', ls=styles[1], lw=.8)
ax2.plot(pres_range, delta_pres3, color=color, label=f'{wavelength[2]} Å', ls=styles[2], lw=.8)
ax2.tick_params(axis='y', labelcolor=color)


# Get handles and labels for all data sets
handles1, _ = ax1.get_legend_handles_labels()
handles2, labels2 = ax2.get_legend_handles_labels()

# get the solid color lines from the 1st and 4th handles of ax1 and the 1st handle of ax2
color_handles = [handles1[0], handles1[3], handles2[0]]

# make new lines for linestyle legend
ls_handle1 = plt.Line2D([0], [0], ls='-', color='black', markersize=10)
ls_handle2 = plt.Line2D([0], [0], ls='--', color='black', markersize=10)
ls_handle3 = plt.Line2D([0], [0], ls='-.', color='black', markersize=10)

# setup legends 
color_legend = ax1.legend(handles=color_handles, labels=color_labels, loc='upper center')

linestyle_legend = ax2.legend([ls_handle1, ls_handle2, ls_handle3], labels2, loc='upper left')

ax1.add_artist(color_legend)
ax2.add_artist(linestyle_legend)

# save plot as png
fig.tight_layout() 
plt.savefig('pres_standard_air.png')
plt.close()


""" Varying Pressure Graph at Standard Dry Air Conditions vs Old Refractive Index """

# pressure graph
pres_range = np.arange(40000, 130000, 100)

# get refractive index and change in wavelength 
n_pres1 = []
n_pres2 = []
n_pres3 = []

delta_pres1 = []
delta_pres2 = []
delta_pres3 = []
for pressure in pres_range:
    # standard dry air has T=288.15[K], 101325 [Pa], 0% humidity, 450 ppm CO2
    n_prop = n_prop_func(288.15, pressure, 0, 450)
    n_pres1.append(n_prop[0])
    n_pres2.append(n_prop[1])
    n_pres3.append(n_prop[2])

    delta_wavelength = np.abs(wavelength * n_prop - wavelength * n_axs)
    delta_pres1.append(delta_wavelength[0])
    delta_pres2.append(delta_wavelength[1])
    delta_pres3.append(delta_wavelength[2])

# setup plot
fig, ax1 = plt.subplots()
plt.ticklabel_format(axis='y', style='plain', useOffset=False)
plt.title("New vs Old Refractive Index w/ Varying Pressure at Std. Air")

# graph the original n vals
n_axs_pres = np.zeros(len(pres_range)) + n_axs[0]
ax1.plot(pres_range, n_axs_pres, color='g',label=wavelength[0], ls=styles[0], lw=.8)

n_axs_pres = np.zeros(len(pres_range)) + n_axs[1]
ax1.plot(pres_range, n_axs_pres, color='g', label=wavelength[1], ls=styles[1], lw=.8)

n_axs_pres = np.zeros(len(pres_range)) + n_axs[2]
ax1.plot(pres_range, n_axs_pres, color='g', label=wavelength[2], ls=styles[2], lw=.8)

# graph the Refracive Index vs Pressure for each wavelength
color = 'tab:orange'
ax1.set_xlabel('Pressure [Pa]')
ax1.set_ylabel('Refractive Index', color=color)
ax1.plot(pres_range, n_pres1, color=color, label=wavelength[0], ls=styles[0], lw=.8)
ax1.plot(pres_range, n_pres2, color=color, label=wavelength[1], ls=styles[1], lw=.8)
ax1.plot(pres_range, n_pres3, color=color, label=wavelength[2], ls=styles[2], lw=.8)
ax1.tick_params(axis='y', labelcolor=color)
ax2 = ax1.twinx()

# graph the change in wavelength vs Pressure for each wavelength
color = 'tab:blue'
ax2.set_ylabel('Old and New Model Wavelength Adjustments Differnece [Å]', fontsize=9, color=color)
ax2.plot(pres_range, delta_pres1, color=color, label=f'{wavelength[0]} Å', ls=styles[0], lw=.8)
ax2.plot(pres_range, delta_pres2, color=color, label=f'{wavelength[1]} Å', ls=styles[1], lw=.8)
ax2.plot(pres_range, delta_pres3, color=color, label=f'{wavelength[2]} Å', ls=styles[2], lw=.8)
ax2.tick_params(axis='y', labelcolor=color)


# Get handles and labels for all data sets
handles1, _ = ax1.get_legend_handles_labels()
handles2, labels2 = ax2.get_legend_handles_labels()

# get the solid color lines from the 1st and 4th handles of ax1 and the 1st handle of ax2
color_handles = [handles1[0], handles1[3], handles2[0]]

# make new lines for linestyle legend
ls_handle1 = plt.Line2D([0], [0], ls='-', color='black', markersize=10)
ls_handle2 = plt.Line2D([0], [0], ls='--', color='black', markersize=10)
ls_handle3 = plt.Line2D([0], [0], ls='-.', color='black', markersize=10)

# setup legends 
color_legend = ax1.legend(handles=color_handles, labels=color_labels, loc='upper center')

linestyle_legend = ax2.legend([ls_handle1, ls_handle2, ls_handle3], labels2, loc='upper right')

ax1.add_artist(color_legend)
ax2.add_artist(linestyle_legend)

# save plot as png
fig.tight_layout() 
plt.savefig('pres_standard_air_vs_old.png')
plt.close()


""" Varying Pressure Graph at Summit Air Conditions """

# get refractive index and change in wavelength 
n_pres1 = []
n_pres2 = []
n_pres3 = []

delta_pres1 = []
delta_pres2 = []
delta_pres3 = []
for pressure in pres_range:
    # temperature, pressure, humidity from KR.20240505.22262.90_icubes.fits
    n_prop = n_prop_func(275.25, pressure, .154, 426.9)
    n_pres1.append(n_prop[0])
    n_pres2.append(n_prop[1])
    n_pres3.append(n_prop[2])

    delta_wavelength = np.abs(wavelength * n_prop - wavelength)
    delta_pres1.append(delta_wavelength[0])
    delta_pres2.append(delta_wavelength[1])
    delta_pres3.append(delta_wavelength[2])

# print slopes
slopes = []
slopes.append(np.polyfit(pres_range, delta_pres1, 1)[0])
slopes.append(np.polyfit(pres_range, delta_pres2, 1)[0])
slopes.append(np.polyfit(pres_range, delta_pres3, 1)[0])
print(f'Slopes of delta_pres/pres_range: {slopes}')

# setup plot
fig, ax1 = plt.subplots()
plt.ticklabel_format(axis='y', style='plain', useOffset=False)
plt.title("Varying Pressure [Pa] at Summit Air Conditions")

# plot average condition range on summit
plt.axvline(x=50000, color='darkgrey')
plt.axvline(x=70000, color='darkgrey')

# graph the original n vals
n_axs_pres = np.zeros(len(pres_range)) + n_axs[0]
ax1.plot(pres_range, n_axs_pres, color='g',label=wavelength[0], ls=styles[0], lw=.8)

n_axs_pres = np.zeros(len(pres_range)) + n_axs[1]
ax1.plot(pres_range, n_axs_pres, color='g', label=wavelength[1], ls=styles[1], lw=.8)

n_axs_pres = np.zeros(len(pres_range)) + n_axs[2]
ax1.plot(pres_range, n_axs_pres, color='g', label=wavelength[2], ls=styles[2], lw=.8)

# graph the Refracive Index vs Pressure for each wavelength
color = 'tab:orange'
ax1.set_xlabel('Pressure [Pa]')
ax1.set_ylabel('Refractive Index', color=color)
ax1.plot(pres_range, n_pres1, color=color, label=wavelength[0], ls=styles[0], lw=.8)
ax1.plot(pres_range, n_pres2, color=color, label=wavelength[1], ls=styles[1], lw=.8)
ax1.plot(pres_range, n_pres3, color=color, label=wavelength[2], ls=styles[2], lw=.8)
ax1.tick_params(axis='y', labelcolor=color)
ax2 = ax1.twinx()

# graph the change in wavelength vs Pressure for each wavelength
color = 'tab:blue'
ax2.set_ylabel('Change in Wavelength [Å]', color=color)
ax2.plot(pres_range, delta_pres1, color=color, label=f'{wavelength[0]} Å', ls=styles[0], lw=.8)
ax2.plot(pres_range, delta_pres2, color=color, label=f'{wavelength[1]} Å', ls=styles[1], lw=.8)
ax2.plot(pres_range, delta_pres3, color=color, label=f'{wavelength[2]} Å', ls=styles[2], lw=.8)
ax2.tick_params(axis='y', labelcolor=color)


# Get handles and labels for all data sets
handles1, _ = ax1.get_legend_handles_labels()
handles2, labels2 = ax2.get_legend_handles_labels()

# get the solid color lines from the 1st and 4th handles of ax1 and the 1st handle of ax2
color_handles = [handles1[0], handles1[3], handles2[0]]

# make new lines for linestyle legend
ls_handle1 = plt.Line2D([0], [0], ls='-', color='black', markersize=10)
ls_handle2 = plt.Line2D([0], [0], ls='--', color='black', markersize=10)
ls_handle3 = plt.Line2D([0], [0], ls='-.', color='black', markersize=10)

# setup legends 
color_legend = ax1.legend(handles=color_handles, labels=color_labels, loc='upper center')

linestyle_legend = ax2.legend([ls_handle1, ls_handle2, ls_handle3], labels2, loc='upper left')

ax1.add_artist(color_legend)
ax2.add_artist(linestyle_legend)

# save plot as png
fig.tight_layout() 
plt.savefig('pres_summit_air.png')
plt.close()


""" Varying Pressure Graph at Summit Air Conditions vs Old"""

# get refractive index and change in wavelength 
n_pres1 = []
n_pres2 = []
n_pres3 = []

delta_pres1 = []
delta_pres2 = []
delta_pres3 = []
for pressure in pres_range:
    # temperature, pressure, humidity from KR.20240505.22262.90_icubes.fits
    n_prop = n_prop_func(275.25, pressure, .154, 426.9)
    n_pres1.append(n_prop[0])
    n_pres2.append(n_prop[1])
    n_pres3.append(n_prop[2])

    delta_wavelength = np.abs(wavelength * n_prop - wavelength * n_axs)
    delta_pres1.append(delta_wavelength[0])
    delta_pres2.append(delta_wavelength[1])
    delta_pres3.append(delta_wavelength[2])

# setup plot
fig, ax1 = plt.subplots()
plt.ticklabel_format(axis='y', style='plain', useOffset=False)
plt.title("New vs Old Refractive Index w/ Varying Pressure at Summit Air")

# plot average condition range on summit
plt.axvline(x=50000, color='darkgrey')
plt.axvline(x=70000, color='darkgrey')

# graph the original n vals
n_axs_pres = np.zeros(len(pres_range)) + n_axs[0]
ax1.plot(pres_range, n_axs_pres, color='g',label=wavelength[0], ls=styles[0], lw=.8)

n_axs_pres = np.zeros(len(pres_range)) + n_axs[1]
ax1.plot(pres_range, n_axs_pres, color='g', label=wavelength[1], ls=styles[1], lw=.8)

n_axs_pres = np.zeros(len(pres_range)) + n_axs[2]
ax1.plot(pres_range, n_axs_pres, color='g', label=wavelength[2], ls=styles[2], lw=.8)

# graph the Refracive Index vs Pressure for each wavelength
color = 'tab:orange'
ax1.set_xlabel('Pressure [Pa]')
ax1.set_ylabel('Refractive Index', color=color)
ax1.plot(pres_range, n_pres1, color=color, label=wavelength[0], ls=styles[0], lw=.8)
ax1.plot(pres_range, n_pres2, color=color, label=wavelength[1], ls=styles[1], lw=.8)
ax1.plot(pres_range, n_pres3, color=color, label=wavelength[2], ls=styles[2], lw=.8)
ax1.tick_params(axis='y', labelcolor=color)
ax2 = ax1.twinx()

# graph the change in wavelength vs Pressure for each wavelength
color = 'tab:blue'
ax2.set_ylabel('Old and New Model Wavelength Adjustments Differnece [Å]', fontsize=9, color=color)
ax2.plot(pres_range, delta_pres1, color=color, label=f'{wavelength[0]} Å', ls=styles[0], lw=.8)
ax2.plot(pres_range, delta_pres2, color=color, label=f'{wavelength[1]} Å', ls=styles[1], lw=.8)
ax2.plot(pres_range, delta_pres3, color=color, label=f'{wavelength[2]} Å', ls=styles[2], lw=.8)
ax2.tick_params(axis='y', labelcolor=color)


# Get handles and labels for all data sets
handles1, _ = ax1.get_legend_handles_labels()
handles2, labels2 = ax2.get_legend_handles_labels()

# get the solid color lines from the 1st and 4th handles of ax1 and the 1st handle of ax2
color_handles = [handles1[0], handles1[3], handles2[0]]

# make new lines for linestyle legend
ls_handle1 = plt.Line2D([0], [0], ls='-', color='black', markersize=10)
ls_handle2 = plt.Line2D([0], [0], ls='--', color='black', markersize=10)
ls_handle3 = plt.Line2D([0], [0], ls='-.', color='black', markersize=10)

# setup legends 
color_legend = ax1.legend(handles=color_handles, labels=color_labels, loc='upper center')

linestyle_legend = ax2.legend([ls_handle1, ls_handle2, ls_handle3], labels2, loc='upper right')

ax1.add_artist(color_legend)
ax2.add_artist(linestyle_legend)

# save plot as png
fig.tight_layout() 
plt.savefig('pres_summit_air_vs_old.png')
plt.close()


""" Varying Humidity Graph at Standard Dry Air Conditions """
hum_range = np.arange(0, 1, .01)

# get refractive index and change in wavelength 
n_hum1 = []
n_hum2 = []
n_hum3 = []

delta_hum1 = []
delta_hum2 = []
delta_hum3 = []
for humidity in hum_range:
    # standard dry air has T=288.15[K], 101325 [Pa], 0% humidity, 450 ppm CO2
    n_prop = n_prop_func(288.15, 101325, humidity, 450)
    n_hum1.append(n_prop[0])
    n_hum2.append(n_prop[1])
    n_hum3.append(n_prop[2])

    delta_wavelength = np.abs(wavelength * n_prop - wavelength)
    delta_hum1.append(delta_wavelength[0])
    delta_hum2.append(delta_wavelength[1])
    delta_hum3.append(delta_wavelength[2])

# setup plot
fig, ax1 = plt.subplots()
plt.ticklabel_format(axis='y', style='plain', useOffset=False)
plt.title("Varying Fractional Humidity at Standard Dry Air Conditions")

# graph the original n vals
n_axs_hum = np.zeros(len(hum_range)) + n_axs[0]
ax1.plot(hum_range, n_axs_hum, color='g',label=wavelength[0], ls=styles[0], lw=.8)

n_axs_hum = np.zeros(len(hum_range)) + n_axs[1]
ax1.plot(hum_range, n_axs_hum, color='g', label=wavelength[1], ls=styles[1], lw=.8)

n_axs_hum = np.zeros(len(hum_range)) + n_axs[2]
ax1.plot(hum_range, n_axs_hum, color='g', label=wavelength[2], ls=styles[2], lw=.8)

# graph the Refracive Index vs Humidity for each wavelength
color = 'tab:orange'
ax1.set_xlabel('Fractional Humidity')
ax1.set_ylabel('Refractive Index', color=color)
ax1.plot(hum_range, n_hum1, color=color, label=wavelength[0], ls=styles[0], lw=.8)
ax1.plot(hum_range, n_hum2, color=color, label=wavelength[1], ls=styles[1], lw=.8)
ax1.plot(hum_range, n_hum3, color=color, label=wavelength[2], ls=styles[2], lw=.8)
ax1.tick_params(axis='y', labelcolor=color)
ax1.set_ylim(top=1.0003)

ax2 = ax1.twinx()

# graph the change in wavelength vs humidity for each wavelength
color = 'tab:blue'
ax2.set_ylabel('Change in Wavelength [Å]', color=color)
ax2.plot(hum_range, delta_hum1, color=color, label=f'{wavelength[0]} Å', ls=styles[0], lw=.8)
ax2.plot(hum_range, delta_hum2, color=color, label=f'{wavelength[1]} Å', ls=styles[1], lw=.8)
ax2.plot(hum_range, delta_hum3, color=color, label=f'{wavelength[2]} Å', ls=styles[2], lw=.8)
ax2.tick_params(axis='y', labelcolor=color)
ax2.set_ylim(top=2.9)

# Get handles and labels for all data sets
handles1, _ = ax1.get_legend_handles_labels()
handles2, labels2 = ax2.get_legend_handles_labels()

# get the solid color lines from the 1st and 4th handles of ax1 and the 1st handle of ax2
color_handles = [handles1[0], handles1[3], handles2[0]]

# make new lines for linestyle legend
ls_handle1 = plt.Line2D([0], [0], ls='-', color='black', markersize=10)
ls_handle2 = plt.Line2D([0], [0], ls='--', color='black', markersize=10)
ls_handle3 = plt.Line2D([0], [0], ls='-.', color='black', markersize=10)

# setup legends 
color_legend = ax1.legend(handles=color_handles, labels=color_labels, loc='upper center')

linestyle_legend = ax2.legend([ls_handle1, ls_handle2, ls_handle3], labels2, loc='upper right')

ax1.add_artist(color_legend)
ax2.add_artist(linestyle_legend)

# save plot as png
fig.tight_layout() 
plt.savefig('hum_standard_air.png')
plt.close()


""" Varying Humidity Graph at Standard Dry Air Conditions vs Old Refractive Index """
hum_range = np.arange(0, 1, .01)

# get refractive index and change in wavelength 
n_hum1 = []
n_hum2 = []
n_hum3 = []

delta_hum1 = []
delta_hum2 = []
delta_hum3 = []
for humidity in hum_range:
    # standard dry air has T=288.15[K], 101325 [Pa], 0% humidity, 450 ppm CO2
    n_prop = n_prop_func(288.15, 101325, humidity, 450)
    n_hum1.append(n_prop[0])
    n_hum2.append(n_prop[1])
    n_hum3.append(n_prop[2])

    delta_wavelength = np.abs(wavelength * n_prop - wavelength * n_axs)
    delta_hum1.append(delta_wavelength[0])
    delta_hum2.append(delta_wavelength[1])
    delta_hum3.append(delta_wavelength[2])

# setup plot
fig, ax1 = plt.subplots()
plt.ticklabel_format(axis='y', style='plain', useOffset=False)
plt.title("New vs Old Refractive Index w/ Varying Humidity at Std. Air")

# graph the original n vals
n_axs_hum = np.zeros(len(hum_range)) + n_axs[0]
ax1.plot(hum_range, n_axs_hum, color='g',label=wavelength[0], ls=styles[0], lw=.8)

n_axs_hum = np.zeros(len(hum_range)) + n_axs[1]
ax1.plot(hum_range, n_axs_hum, color='g', label=wavelength[1], ls=styles[1], lw=.8)

n_axs_hum = np.zeros(len(hum_range)) + n_axs[2]
ax1.plot(hum_range, n_axs_hum, color='g', label=wavelength[2], ls=styles[2], lw=.8)

# graph the Refracive Index vs Humidity for each wavelength
color = 'tab:orange'
ax1.set_xlabel('Fractional Humidity')
ax1.set_ylabel('Refractive Index', color=color)
ax1.plot(hum_range, n_hum1, color=color, label=wavelength[0], ls=styles[0], lw=.8)
ax1.plot(hum_range, n_hum2, color=color, label=wavelength[1], ls=styles[1], lw=.8)
ax1.plot(hum_range, n_hum3, color=color, label=wavelength[2], ls=styles[2], lw=.8)
ax1.tick_params(axis='y', labelcolor=color)
ax1.set_ylim(top=1.0003)

ax2 = ax1.twinx()

# graph the change in wavelength vs humidity for each wavelength
color = 'tab:blue'
ax2.set_ylabel('Old and New Model Wavelength Adjustments Differnece [Å]', fontsize=9, color=color)
ax2.plot(hum_range, delta_hum1, color=color, label=f'{wavelength[0]} Å', ls=styles[0], lw=.8)
ax2.plot(hum_range, delta_hum2, color=color, label=f'{wavelength[1]} Å', ls=styles[1], lw=.8)
ax2.plot(hum_range, delta_hum3, color=color, label=f'{wavelength[2]} Å', ls=styles[2], lw=.8)
ax2.tick_params(axis='y', labelcolor=color)
ax2.set_ylim(top=.008)

# Get handles and labels for all data sets
handles1, _ = ax1.get_legend_handles_labels()
handles2, labels2 = ax2.get_legend_handles_labels()

# get the solid color lines from the 1st and 4th handles of ax1 and the 1st handle of ax2
color_handles = [handles1[0], handles1[3], handles2[0]]

# make new lines for linestyle legend
ls_handle1 = plt.Line2D([0], [0], ls='-', color='black', markersize=10)
ls_handle2 = plt.Line2D([0], [0], ls='--', color='black', markersize=10)
ls_handle3 = plt.Line2D([0], [0], ls='-.', color='black', markersize=10)

# setup legends 
color_legend = ax1.legend(handles=color_handles, labels=color_labels, loc='upper center')

linestyle_legend = ax2.legend([ls_handle1, ls_handle2, ls_handle3], labels2, loc='upper left')

ax1.add_artist(color_legend)
ax2.add_artist(linestyle_legend)

# save plot as png
fig.tight_layout() 
plt.savefig('hum_standard_air_vs_old.png')
plt.close()


""" Varying Humidity Graph at Summit Air Conditions vs Old Refractive Index"""

# get refractive index and change in wavelength 
n_hum1 = []
n_hum2 = []
n_hum3 = []

delta_hum1 = []
delta_hum2 = []
delta_hum3 = []
for humidity in hum_range:
    # temperature, pressure, humidity from KR.20240505.22262.90_icubes.fits
    n_prop = n_prop_func(275.25, 61785, humidity, 426.9)
    n_hum1.append(n_prop[0])
    n_hum2.append(n_prop[1])
    n_hum3.append(n_prop[2])

    delta_wavelength = np.abs(wavelength * n_prop - wavelength * n_axs)
    delta_hum1.append(delta_wavelength[0])
    delta_hum2.append(delta_wavelength[1])
    delta_hum3.append(delta_wavelength[2])

# setup plot
fig, ax1 = plt.subplots()
plt.ticklabel_format(axis='y', style='plain', useOffset=False)
plt.title("New vs Old Refractive Index w/ Varying Humidity at Summit Air")

# plot average condition range on summit
plt.axvline(x=.1, color='darkgrey')
plt.axvline(x=.3, color='darkgrey')

# graph the original n vals
n_axs_hum = np.zeros(len(hum_range)) + n_axs[0]
ax1.plot(hum_range, n_axs_hum, color='g',label=wavelength[0], ls=styles[0], lw=.8)

n_axs_hum = np.zeros(len(hum_range)) + n_axs[1]
ax1.plot(hum_range, n_axs_hum, color='g', label=wavelength[1], ls=styles[1], lw=.8)

n_axs_hum = np.zeros(len(hum_range)) + n_axs[2]
ax1.plot(hum_range, n_axs_hum, color='g', label=wavelength[2], ls=styles[2], lw=.8)

# graph the Refracive Index vs Humidity for each wavelength
color = 'tab:orange'
ax1.set_xlabel('Fractional Humidity')
ax1.set_ylabel('Refractive Index', color=color)
ax1.plot(hum_range, n_hum1, color=color, label=wavelength[0], ls=styles[0], lw=.8)
ax1.plot(hum_range, n_hum2, color=color, label=wavelength[1], ls=styles[1], lw=.8)
ax1.plot(hum_range, n_hum3, color=color, label=wavelength[2], ls=styles[2], lw=.8)
ax1.tick_params(axis='y', labelcolor=color)
ax1.set_ylim(top=1.00034)

ax2 = ax1.twinx()

# graph the change in wavelength vs Humidity for each wavelength
color = 'tab:blue'
ax2.set_ylabel('Old and New Model Wavelength Adjustments Differnece [Å]', fontsize=9, color=color)
ax2.plot(hum_range, delta_hum1, color=color, label=f'{wavelength[0]} Å', ls=styles[0], lw=.8)
ax2.plot(hum_range, delta_hum2, color=color, label=f'{wavelength[1]} Å', ls=styles[1], lw=.8)
ax2.plot(hum_range, delta_hum3, color=color, label=f'{wavelength[2]} Å', ls=styles[2], lw=.8)
ax2.tick_params(axis='y', labelcolor=color)
ax2.set_ylim(top=1.4)


# Get handles and labels for all data sets
handles1, _ = ax1.get_legend_handles_labels()
handles2, labels2 = ax2.get_legend_handles_labels()

# get the solid color lines from the 1st and 4th handles of ax1 and the 1st handle of ax2
color_handles = [handles1[0], handles1[3], handles2[0]]

# make new lines for linestyle legend
ls_handle1 = plt.Line2D([0], [0], ls='-', color='black', markersize=10)
ls_handle2 = plt.Line2D([0], [0], ls='--', color='black', markersize=10)
ls_handle3 = plt.Line2D([0], [0], ls='-.', color='black', markersize=10)

# setup legends 
color_legend = ax1.legend(handles=color_handles, labels=color_labels, loc='upper center')

linestyle_legend = ax2.legend([ls_handle1, ls_handle2, ls_handle3], labels2, loc='upper right')

ax1.add_artist(color_legend)
ax2.add_artist(linestyle_legend)

# save plot as png
fig.tight_layout() 
plt.savefig('hum_summit_air_vs_old.png')
plt.close()


""" Varying Humidity Graph at Summit Air Conditions """

# get refractive index and change in wavelength 
n_hum1 = []
n_hum2 = []
n_hum3 = []

delta_hum1 = []
delta_hum2 = []
delta_hum3 = []
for humidity in hum_range:
    # temperature, pressure, humidity from KR.20240505.22262.90_icubes.fits
    n_prop = n_prop_func(275.25, 61785, humidity, 426.9)
    n_hum1.append(n_prop[0])
    n_hum2.append(n_prop[1])
    n_hum3.append(n_prop[2])

    delta_wavelength = np.abs(wavelength * n_prop - wavelength)
    delta_hum1.append(delta_wavelength[0])
    delta_hum2.append(delta_wavelength[1])
    delta_hum3.append(delta_wavelength[2])

# print slopes
slopes = []
slopes.append(np.polyfit(hum_range, delta_hum1, 1)[0])
slopes.append(np.polyfit(hum_range, delta_hum2, 1)[0])
slopes.append(np.polyfit(hum_range, delta_hum3, 1)[0])
print(f'Slopes of delta_hum/hum_range: {slopes}')

# setup plot
fig, ax1 = plt.subplots()
plt.ticklabel_format(axis='y', style='plain', useOffset=False)
plt.title("Varying Fractional Humidity at Summit Air Conditions")

# plot average condition range on summit
plt.axvline(x=.1, color='darkgrey')
plt.axvline(x=.3, color='darkgrey')

# graph the original n vals
n_axs_hum = np.zeros(len(hum_range)) + n_axs[0]
ax1.plot(hum_range, n_axs_hum, color='g',label=wavelength[0], ls=styles[0], lw=.8)

n_axs_hum = np.zeros(len(hum_range)) + n_axs[1]
ax1.plot(hum_range, n_axs_hum, color='g', label=wavelength[1], ls=styles[1], lw=.8)

n_axs_hum = np.zeros(len(hum_range)) + n_axs[2]
ax1.plot(hum_range, n_axs_hum, color='g', label=wavelength[2], ls=styles[2], lw=.8)

# graph the Refracive Index vs Humidity for each wavelength
color = 'tab:orange'
ax1.set_xlabel('Fractional Humidity')
ax1.set_ylabel('Refractive Index', color=color)
ax1.plot(hum_range, n_hum1, color=color, label=wavelength[0], ls=styles[0], lw=.8)
ax1.plot(hum_range, n_hum2, color=color, label=wavelength[1], ls=styles[1], lw=.8)
ax1.plot(hum_range, n_hum3, color=color, label=wavelength[2], ls=styles[2], lw=.8)
ax1.tick_params(axis='y', labelcolor=color)
ax1.set_ylim(top=1.00034)

ax2 = ax1.twinx()

# graph the change in wavelength vs Humidity for each wavelength
color = 'tab:blue'
ax2.set_ylabel('Change in Wavelength [Å]', color=color)
ax2.plot(hum_range, delta_hum1, color=color, label=f'{wavelength[0]} Å', ls=styles[0], lw=.8)
ax2.plot(hum_range, delta_hum2, color=color, label=f'{wavelength[1]} Å', ls=styles[1], lw=.8)
ax2.plot(hum_range, delta_hum3, color=color, label=f'{wavelength[2]} Å', ls=styles[2], lw=.8)
ax2.tick_params(axis='y', labelcolor=color)
ax2.set_ylim(top=1.4)

# Get handles and labels for all data sets
handles1, _ = ax1.get_legend_handles_labels()
handles2, labels2 = ax2.get_legend_handles_labels()

# get the solid color lines from the 1st and 4th handles of ax1 and the 1st handle of ax2
color_handles = [handles1[0], handles1[3], handles2[0]]

# make new lines for linestyle legend
ls_handle1 = plt.Line2D([0], [0], ls='-', color='black', markersize=10)
ls_handle2 = plt.Line2D([0], [0], ls='--', color='black', markersize=10)
ls_handle3 = plt.Line2D([0], [0], ls='-.', color='black', markersize=10)

# setup legends 
color_legend = ax1.legend(handles=color_handles, labels=color_labels, loc='upper center')

linestyle_legend = ax2.legend([ls_handle1, ls_handle2, ls_handle3], labels2, loc='upper right')

ax1.add_artist(color_legend)
ax2.add_artist(linestyle_legend)

# save plot as png
fig.tight_layout() 
plt.savefig('hum_summit_air.png')
plt.close()


""" Varying Humidity Graph at Summit Air Conditions @ 3000 Å """
# setup plot

fig, ax1 = plt.subplots()
plt.ticklabel_format(axis='y', style='plain', useOffset=False)
plt.title("Varying Fractional Humidity at Summit Air Conditions at 3000 Å")

# plot average condition range on summit
plt.axvline(x=.1, color='darkgrey')
plt.axvline(x=.3, color='darkgrey')

# graph the original n vals
n_axs_hum = np.zeros(len(hum_range)) + n_axs[0]
# ax1.plot(hum_range, n_axs_hum, color='g',label=wavelength[0], ls=styles[0], lw=.8)

# graph the Refracive Index vs Humidity for each wavelength
color = 'tab:orange'
ax1.set_xlabel('Fractional Humidity')
ax1.set_ylabel('Refractive Index', color=color)
ax1.plot(hum_range, n_hum1, color=color, label=wavelength[0], ls=styles[0], lw=.8)

ax1.tick_params(axis='y', labelcolor=color)


ax2 = ax1.twinx()

# graph the change in wavelength vs Humidity for each wavelength
color = 'tab:blue'
ax2.set_ylabel('Change in Wavelength [Å]', color=color)
ax2.plot(hum_range, delta_hum1, color=color, label=f'{wavelength[0]} Å', ls=(0, (4,4)), lw=.8)
ax2.tick_params(axis='y', labelcolor=color)


# Get handles and labels for all data sets
handles1, _ = ax1.get_legend_handles_labels()
handles2, labels2 = ax2.get_legend_handles_labels()

# get the solid color lines from handles of ax1 and ax2
color_handles = [handles1[0], handles2[0]]

# make new lines for linestyle legend
ls_handle1 = plt.Line2D([0], [0], ls='-', color='black', markersize=10)

# setup legends 
color_legend = ax1.legend(handles=color_handles, labels=[color_labels[1], color_labels[2]], loc='upper center')

linestyle_legend = ax2.legend([ls_handle1], labels2, loc='upper right')

ax1.add_artist(color_legend)
ax2.add_artist(linestyle_legend)

# save plot as png
fig.tight_layout() 
plt.savefig('hum_summit_air_3000A.png')
plt.close()


""" Varying Humidity Graph at Summit Air Conditions @ 6000A """
# setup plot

fig, ax1 = plt.subplots()
plt.ticklabel_format(axis='y', style='plain', useOffset=False)
plt.title("Varying Fractional Humidity at Summit Air Conditions at 6000 Å")

# plot average condition range on summit
plt.axvline(x=.1, color='darkgrey')
plt.axvline(x=.3, color='darkgrey')

# graph the original n vals
n_axs_hum = np.zeros(len(hum_range)) + n_axs[1]
# ax1.plot(hum_range, n_axs_hum, color='g',label=wavelength[0], ls=styles[0], lw=.8)

# graph the Refracive Index vs Humidity for each wavelength
color = 'tab:orange'
ax1.set_xlabel('Fractional Humidity')
ax1.set_ylabel('Refractive Index', color=color)
ax1.plot(hum_range, n_hum2, color=color, label=wavelength[1], ls=styles[0], lw=.8)

ax1.tick_params(axis='y', labelcolor=color)


ax2 = ax1.twinx()

# graph the change in wavelength vs Humidity for each wavelength
color = 'tab:blue'
ax2.set_ylabel('Change in Wavelength [Å]', color=color)
ax2.plot(hum_range, delta_hum2, color=color, label=f'{wavelength[1]} Å', ls=(0, (4,4)), lw=.8)
ax2.tick_params(axis='y', labelcolor=color)


# Get handles and labels for all data sets
handles1, _ = ax1.get_legend_handles_labels()
handles2, labels2 = ax2.get_legend_handles_labels()

# get the solid color lines from handles of ax1 and ax2
color_handles = [handles1[0], handles2[0]]

# make new lines for linestyle legend
ls_handle1 = plt.Line2D([0], [0], ls='-', color='black', markersize=10)

# setup legends 
color_legend = ax1.legend(handles=color_handles, labels=[color_labels[1], color_labels[2]], loc='upper center')

linestyle_legend = ax2.legend([ls_handle1], labels2, loc='upper right')

ax1.add_artist(color_legend)
ax2.add_artist(linestyle_legend)

# save plot as png
fig.tight_layout() 
plt.savefig('hum_summit_air_6000A.png')
plt.close()


""" Varying Humidity Graph at Summit Air Conditions @ 9000A """
# setup plot

fig, ax1 = plt.subplots()
plt.ticklabel_format(axis='y', style='plain', useOffset=False)
plt.title("Varying Fractional Humidity at Summit Air Conditions at 9000 Å")

# plot average condition range on summit
plt.axvline(x=.1, color='darkgrey')
plt.axvline(x=.3, color='darkgrey')

# graph the original n vals
n_axs_hum = np.zeros(len(hum_range)) + n_axs[2]
# ax1.plot(hum_range, n_axs_hum, color='g',label=wavelength[0], ls=styles[0], lw=.8)

# graph the Refracive Index vs Humidity for each wavelength
color = 'tab:orange'
ax1.set_xlabel('Fractional Humidity')
ax1.set_ylabel('Refractive Index', color=color)
ax1.plot(hum_range, n_hum3, color=color, label=wavelength[2], ls=styles[0], lw=.8)

ax1.tick_params(axis='y', labelcolor=color)


ax2 = ax1.twinx()

# graph the change in wavelength vs Humidity for each wavelength
color = 'tab:blue'
ax2.set_ylabel('Change in Wavelength [Å]', color=color)
ax2.plot(hum_range, delta_hum3, color=color, label=f'{wavelength[2]} Å', ls=(0, (4,4)), lw=.8)
ax2.tick_params(axis='y', labelcolor=color)


# Get handles and labels for all data sets
handles1, _ = ax1.get_legend_handles_labels()
handles2, labels2 = ax2.get_legend_handles_labels()

# get the solid color lines from handles of ax1 and ax2
color_handles = [handles1[0], handles2[0]]

# make new lines for linestyle legend
ls_handle1 = plt.Line2D([0], [0], ls='-', color='black', markersize=10)

# setup legends 
color_legend = ax1.legend(handles=color_handles, labels=[color_labels[1], color_labels[2]], loc='upper center')

linestyle_legend = ax2.legend([ls_handle1], labels2, loc='upper right')

ax1.add_artist(color_legend)
ax2.add_artist(linestyle_legend)

# save plot as png
fig.tight_layout() 
plt.savefig('hum_summit_air_9000A.png')
plt.close()


""" Varying CO2 Graph at Standard Dry Air Conditions vs Old Refractive Index """

co2_range = np.arange(0, 1000, 10)

# get refractive index and change in wavelength 
n_co2_1 = []
n_co2_2 = []
n_co2_3 = []

delta_co2_1 = []
delta_co2_2 = []
delta_co2_3 = []
for co2 in co2_range:
    # standard dry air has T=288.15[K], 101325 [Pa], 0% humidity, 450 ppm CO2
    n_prop = n_prop_func(288.15, 101325, 0, co2)
    n_co2_1.append(n_prop[0])
    n_co2_2.append(n_prop[1])
    n_co2_3.append(n_prop[2])

    delta_wavelength = np.abs(wavelength * n_prop - wavelength * n_axs)
    delta_co2_1.append(delta_wavelength[0])
    delta_co2_2.append(delta_wavelength[1])
    delta_co2_3.append(delta_wavelength[2])

# setup plot
fig, ax1 = plt.subplots()
plt.ticklabel_format(axis='y', style='plain', useOffset=False)
plt.title("New vs Old Refractive Index w/ Varying ppm CO2 at Std. Air")


# graph the original n vals
n_axs_co2 = np.zeros(len(co2_range)) + n_axs[0]
ax1.plot(co2_range, n_axs_co2, color='g',label=wavelength[0], ls=styles[0], lw=.8)

n_axs_co2 = np.zeros(len(co2_range)) + n_axs[1]
ax1.plot(co2_range, n_axs_co2, color='g', label=wavelength[1], ls=styles[1], lw=.8)

n_axs_co2 = np.zeros(len(co2_range)) + n_axs[2]
ax1.plot(co2_range, n_axs_co2, color='g', label=wavelength[2], ls=styles[2], lw=.8)

# graph the Refracive Index vs ppm of CO2 for each wavelength
color = 'tab:orange'
ax1.set_xlabel('ppm of CO2')
ax1.set_ylabel('Refractive Index', color=color)
ax1.plot(co2_range, n_co2_1, color=color, label=wavelength[0], ls=styles[0], lw=.8)
ax1.plot(co2_range, n_co2_2, color=color, label=wavelength[1], ls=styles[1], lw=.8)
ax1.plot(co2_range, n_co2_3, color=color, label=wavelength[2], ls=styles[2], lw=.8)
ax1.tick_params(axis='y', labelcolor=color)
ax1.set_ylim(top=1.0003)

ax2 = ax1.twinx()

# graph the change in wavelength vs ppm of CO2 for each wavelength
color = 'tab:blue'
ax2.set_ylabel('Old and New Model Wavelength Adjustments Differnece [Å]', fontsize=9, color=color)
ax2.plot(co2_range, delta_co2_1, color=color, label=f'{wavelength[0]} Å', ls=styles[0], lw=.8)
ax2.plot(co2_range, delta_co2_2, color=color, label=f'{wavelength[1]} Å', ls=styles[1], lw=.8)
ax2.plot(co2_range, delta_co2_3, color=color, label=f'{wavelength[2]} Å', ls=styles[2], lw=.8)
ax2.tick_params(axis='y', labelcolor=color)
ax2.set_ylim(top=.0007)


# Get handles and labels for all data sets
handles1, _ = ax1.get_legend_handles_labels()
handles2, labels2 = ax2.get_legend_handles_labels()

# get the solid color lines from the 1st and 4th handles of ax1 and the 1st handle of ax2
color_handles = [handles1[0], handles1[3], handles2[0]]

# make new lines for linestyle legend
ls_handle1 = plt.Line2D([0], [0], ls='-', color='black', markersize=10)
ls_handle2 = plt.Line2D([0], [0], ls='--', color='black', markersize=10)
ls_handle3 = plt.Line2D([0], [0], ls='-.', color='black', markersize=10)

# setup legends 
color_legend = ax1.legend(handles=color_handles, labels=color_labels, loc='upper center')

linestyle_legend = ax2.legend([ls_handle1, ls_handle2, ls_handle3], labels2, loc='upper right')

ax1.add_artist(color_legend)
ax2.add_artist(linestyle_legend)

# save plot as png
fig.tight_layout() 
plt.savefig('co2_standard_air_vs_old.png')
plt.close()


""" Varying CO2 Graph at Standard Dry Air Conditions """

co2_range = np.arange(0, 1000, 10)

# get refractive index and change in wavelength 
n_co2_1 = []
n_co2_2 = []
n_co2_3 = []

delta_co2_1 = []
delta_co2_2 = []
delta_co2_3 = []
for co2 in co2_range:
    # standard dry air has T=288.15[K], 101325 [Pa], 0% humidity, 450 ppm CO2
    n_prop = n_prop_func(288.15, 101325, 0, co2)
    n_co2_1.append(n_prop[0])
    n_co2_2.append(n_prop[1])
    n_co2_3.append(n_prop[2])

    delta_wavelength = np.abs(wavelength * n_prop - wavelength)
    delta_co2_1.append(delta_wavelength[0])
    delta_co2_2.append(delta_wavelength[1])
    delta_co2_3.append(delta_wavelength[2])

# setup plot
fig, ax1 = plt.subplots()
plt.ticklabel_format(axis='y', style='plain', useOffset=False)
plt.title("Varying ppm of CO2 at Standard Dry Air Conditions")


# graph the original n vals
n_axs_co2 = np.zeros(len(co2_range)) + n_axs[0]
ax1.plot(co2_range, n_axs_co2, color='g',label=wavelength[0], ls=styles[0], lw=.8)

n_axs_co2 = np.zeros(len(co2_range)) + n_axs[1]
ax1.plot(co2_range, n_axs_co2, color='g', label=wavelength[1], ls=styles[1], lw=.8)

n_axs_co2 = np.zeros(len(co2_range)) + n_axs[2]
ax1.plot(co2_range, n_axs_co2, color='g', label=wavelength[2], ls=styles[2], lw=.8)

# graph the Refracive Index vs ppm of CO2 for each wavelength
color = 'tab:orange'
ax1.set_xlabel('ppm of CO2')
ax1.set_ylabel('Refractive Index', color=color)
ax1.plot(co2_range, n_co2_1, color=color, label=wavelength[0], ls=styles[0], lw=.8)
ax1.plot(co2_range, n_co2_2, color=color, label=wavelength[1], ls=styles[1], lw=.8)
ax1.plot(co2_range, n_co2_3, color=color, label=wavelength[2], ls=styles[2], lw=.8)
ax1.tick_params(axis='y', labelcolor=color)
ax1.set_ylim(top=1.0003)

ax2 = ax1.twinx()

# graph the change in wavelength vs ppm of CO2 for each wavelength
color = 'tab:blue'
ax2.set_ylabel('Change in Wavelength [Å]', color=color)
ax2.plot(co2_range, delta_co2_1, color=color, label=f'{wavelength[0]} Å', ls=styles[0], lw=.8)
ax2.plot(co2_range, delta_co2_2, color=color, label=f'{wavelength[1]} Å', ls=styles[1], lw=.8)
ax2.plot(co2_range, delta_co2_3, color=color, label=f'{wavelength[2]} Å', ls=styles[2], lw=.8)
ax2.tick_params(axis='y', labelcolor=color)
ax2.set_ylim(top=3)

# Get handles and labels for all data sets
handles1, _ = ax1.get_legend_handles_labels()
handles2, labels2 = ax2.get_legend_handles_labels()

# get the solid color lines from the 1st and 4th handles of ax1 and the 1st handle of ax2
color_handles = [handles1[0], handles1[3], handles2[0]]

# make new lines for linestyle legend
ls_handle1 = plt.Line2D([0], [0], ls='-', color='black', markersize=10)
ls_handle2 = plt.Line2D([0], [0], ls='--', color='black', markersize=10)
ls_handle3 = plt.Line2D([0], [0], ls='-.', color='black', markersize=10)

# setup legends 
color_legend = ax1.legend(handles=color_handles, labels=color_labels, loc='upper center')

linestyle_legend = ax2.legend([ls_handle1, ls_handle2, ls_handle3], labels2, loc='upper right')

ax1.add_artist(color_legend)
ax2.add_artist(linestyle_legend)

# save plot as png
fig.tight_layout() 
plt.savefig('co2_standard_air.png')
plt.close()


""" Varying ppm of CO2 Graph at Summit Air Conditions vs Old Refractive Index """


# get refractive index and change in wavelength 
n_co2_1 = []
n_co2_2 = []
n_co2_3 = []

delta_co2_1 = []
delta_co2_2 = []
delta_co2_3 = []
for co2 in co2_range:
    # temperature, pressure, humidity from KR.20240505.22262.90_icubes.fits
    n_prop = n_prop_func(275.25, 61785, .154, co2)
    n_co2_1.append(n_prop[0])
    n_co2_2.append(n_prop[1])
    n_co2_3.append(n_prop[2])

    delta_wavelength = np.abs(wavelength * n_prop - wavelength * n_axs)
    delta_co2_1.append(delta_wavelength[0])
    delta_co2_2.append(delta_wavelength[1])
    delta_co2_3.append(delta_wavelength[2])

# setup plot
fig, ax1 = plt.subplots()
plt.ticklabel_format(axis='y', style='plain', useOffset=False)
plt.title("New vs Old Refractive Index w/ Varying ppm CO2 at Summit Air")

# plot average condition range on summit
plt.axvline(x=375, color='darkgrey')
plt.axvline(x=500, color='darkgrey')

# graph the original n vals
n_axs_co2 = np.zeros(len(co2_range)) + n_axs[0]
ax1.plot(co2_range, n_axs_co2, color='g',label=wavelength[0], ls=styles[0], lw=.8)

n_axs_co2 = np.zeros(len(co2_range)) + n_axs[1]
ax1.plot(co2_range, n_axs_co2, color='g', label=wavelength[1], ls=styles[1], lw=.8)

n_axs_co2 = np.zeros(len(co2_range)) + n_axs[2]
ax1.plot(co2_range, n_axs_co2, color='g', label=wavelength[2], ls=styles[2], lw=.8)

# graph the Refracive Index vs CO2 for each wavelength
color = 'tab:orange'
ax1.set_xlabel('ppm of CO2')
ax1.set_ylabel('Refractive Index', color=color)
ax1.plot(co2_range, n_co2_1, color=color, label=wavelength[0], ls=styles[0], lw=.8)
ax1.plot(co2_range, n_co2_2, color=color, label=wavelength[1], ls=styles[1], lw=.8)
ax1.plot(co2_range, n_co2_3, color=color, label=wavelength[2], ls=styles[2], lw=.8)
ax1.tick_params(axis='y', labelcolor=color)
ax1.set_ylim(top=1.00034)

ax2 = ax1.twinx()

# graph the change in wavelength vs ppm of CO2 for each wavelength
color = 'tab:blue'
ax2.set_ylabel('Old and New Model Wavelength Adjustments Differnece [Å]', fontsize=9, color=color)
ax2.plot(co2_range, delta_co2_1, color=color, label=f'{wavelength[0]} Å', ls=styles[0], lw=.8)
ax2.plot(co2_range, delta_co2_2, color=color, label=f'{wavelength[1]} Å', ls=styles[1], lw=.8)
ax2.plot(co2_range, delta_co2_3, color=color, label=f'{wavelength[2]} Å', ls=styles[2], lw=.8)
ax2.tick_params(axis='y', labelcolor=color)
ax2.set_ylim(top=1.3)


# Get handles and labels for all data sets
handles1, _ = ax1.get_legend_handles_labels()
handles2, labels2 = ax2.get_legend_handles_labels()

# get the solid color lines from the 1st and 4th handles of ax1 and the 1st handle of ax2
color_handles = [handles1[0], handles1[3], handles2[0]]

# make new lines for linestyle legend
ls_handle1 = plt.Line2D([0], [0], ls='-', color='black', markersize=10)
ls_handle2 = plt.Line2D([0], [0], ls='--', color='black', markersize=10)
ls_handle3 = plt.Line2D([0], [0], ls='-.', color='black', markersize=10)

# setup legends 
color_legend = ax1.legend(handles=color_handles, labels=color_labels, loc='upper center')

linestyle_legend = ax2.legend([ls_handle1, ls_handle2, ls_handle3], labels2, loc='upper left')

ax1.add_artist(color_legend)
ax2.add_artist(linestyle_legend)

# save plot as png
fig.tight_layout() 
plt.savefig('co2_summit_air_vs_old.png')
plt.close()


""" Varying ppm of CO2 Graph at Summit Air Conditions """

# get refractive index and change in wavelength 
n_co2_1 = []
n_co2_2 = []
n_co2_3 = []

delta_co2_1 = []
delta_co2_2 = []
delta_co2_3 = []
for co2 in co2_range:
    # temperature, pressure, humidity from KR.20240505.22262.90_icubes.fits
    n_prop = n_prop_func(275.25, 61785, .154, co2)
    n_co2_1.append(n_prop[0])
    n_co2_2.append(n_prop[1])
    n_co2_3.append(n_prop[2])

    delta_wavelength = np.abs(wavelength * n_prop - wavelength)
    delta_co2_1.append(delta_wavelength[0])
    delta_co2_2.append(delta_wavelength[1])
    delta_co2_3.append(delta_wavelength[2])

# print slopes
slopes = []
slopes.append(np.polyfit(co2_range, delta_co2_1, 1)[0])
slopes.append(np.polyfit(co2_range, delta_co2_2, 1)[0])
slopes.append(np.polyfit(co2_range, delta_co2_3, 1)[0])
print(f'Slopes of delta_co2/co2_range: {slopes}')

# setup plot
fig, ax1 = plt.subplots()
plt.ticklabel_format(axis='y', style='plain', useOffset=False)
plt.title("Varying ppm of CO2 at Summit Air Conditions")

# plot average condition range on summit
plt.axvline(x=375, color='darkgrey')
plt.axvline(x=500, color='darkgrey')

# graph the original n vals
n_axs_co2 = np.zeros(len(co2_range)) + n_axs[0]
ax1.plot(co2_range, n_axs_co2, color='g',label=wavelength[0], ls=styles[0], lw=.8)

n_axs_co2 = np.zeros(len(co2_range)) + n_axs[1]
ax1.plot(co2_range, n_axs_co2, color='g', label=wavelength[1], ls=styles[1], lw=.8)

n_axs_co2 = np.zeros(len(co2_range)) + n_axs[2]
ax1.plot(co2_range, n_axs_co2, color='g', label=wavelength[2], ls=styles[2], lw=.8)

# graph the Refracive Index vs CO2 for each wavelength
color = 'tab:orange'
ax1.set_xlabel('ppm of CO2')
ax1.set_ylabel('Refractive Index', color=color)
ax1.plot(co2_range, n_co2_1, color=color, label=wavelength[0], ls=styles[0], lw=.8)
ax1.plot(co2_range, n_co2_2, color=color, label=wavelength[1], ls=styles[1], lw=.8)
ax1.plot(co2_range, n_co2_3, color=color, label=wavelength[2], ls=styles[2], lw=.8)
ax1.tick_params(axis='y', labelcolor=color)
ax1.set_ylim(top=1.00034)

ax2 = ax1.twinx()

# graph the change in wavelength vs ppm of CO2 for each wavelength
color = 'tab:blue'
ax2.set_ylabel('Change in Wavelength [Å]', color=color)
ax2.plot(co2_range, delta_co2_1, color=color, label=f'{wavelength[0]} Å', ls=styles[0], lw=.8)
ax2.plot(co2_range, delta_co2_2, color=color, label=f'{wavelength[1]} Å', ls=styles[1], lw=.8)
ax2.plot(co2_range, delta_co2_3, color=color, label=f'{wavelength[2]} Å', ls=styles[2], lw=.8)
ax2.tick_params(axis='y', labelcolor=color)
ax2.set_ylim(top=2)

# Get handles and labels for all data sets
handles1, _ = ax1.get_legend_handles_labels()
handles2, labels2 = ax2.get_legend_handles_labels()

# get the solid color lines from the 1st and 4th handles of ax1 and the 1st handle of ax2
color_handles = [handles1[0], handles1[3], handles2[0]]

# make new lines for linestyle legend
ls_handle1 = plt.Line2D([0], [0], ls='-', color='black', markersize=10)
ls_handle2 = plt.Line2D([0], [0], ls='--', color='black', markersize=10)
ls_handle3 = plt.Line2D([0], [0], ls='-.', color='black', markersize=10)

# setup legends 
color_legend = ax1.legend(handles=color_handles, labels=color_labels, loc='upper center')

linestyle_legend = ax2.legend([ls_handle1, ls_handle2, ls_handle3], labels2, loc='upper left')

ax1.add_artist(color_legend)
ax2.add_artist(linestyle_legend)

# save plot as png
fig.tight_layout() 
plt.savefig('co2_summit_air.png')
plt.close()


""" Varying ppm of CO2 Graph at Summit Air Conditions @ 3000 A """


# setup plot
fig, ax1 = plt.subplots()
plt.ticklabel_format(axis='y', style='plain', useOffset=False)
plt.title("Varying ppm of CO2 at Summit Air Conditions at 3000 Å")

# plot average condition range on summit
plt.axvline(x=375, color='darkgrey')
plt.axvline(x=500, color='darkgrey')

# graph the original n vals
n_axs_co2 = np.zeros(len(co2_range)) + n_axs[0]
#ax1.plot(co2_range, n_axs_co2, color='g',label=wavelength[0], ls=styles[0], lw=.8)

# graph the Refracive Index vs CO2 for each wavelength
color = 'tab:orange'
ax1.set_xlabel('ppm of CO2')
ax1.set_ylabel('Refractive Index', color=color)
ax1.plot(co2_range, n_co2_1, color=color, label=wavelength[0], ls=styles[0], lw=.8)
ax1.tick_params(axis='y', labelcolor=color)
# ax1.set_ylim(top=1.00044)

ax2 = ax1.twinx()

# graph the change in wavelength vs ppm of CO2 for each wavelength
color = 'tab:blue'
ax2.set_ylabel('Change in Wavelength [Å]', color=color)
ax2.plot(co2_range, delta_co2_1, color=color, label=f'{wavelength[0]} Å', ls=(0,(4,4)), lw=.8)
ax2.tick_params(axis='y', labelcolor=color)
# ax2.set_ylim(top=1.2018)
ax2.ticklabel_format(useOffset=False)

# Get handles and labels for all data sets
handles1, _ = ax1.get_legend_handles_labels()
handles2, labels2 = ax2.get_legend_handles_labels()

# get the solid color lines from the handles of ax1 and ax2
color_handles = [handles1[0], handles2[0]]

# make new lines for linestyle legend
ls_handle1 = plt.Line2D([0], [0], ls='-', color='black', markersize=10)

# setup legends 
color_legend = ax1.legend(handles=color_handles, labels=[color_labels[1], color_labels[2]], loc='upper center')

linestyle_legend = ax2.legend([ls_handle1, ls_handle2, ls_handle3], labels2, loc='upper left')

ax1.add_artist(color_legend)
ax2.add_artist(linestyle_legend)

# save plot as png
fig.tight_layout() 
plt.savefig('co2_summit_air_3000A.png')
plt.close()


""" Varying ppm of CO2 Graph at Summit Air Conditions @ 6000 A """

# setup plot
fig, ax1 = plt.subplots()
plt.ticklabel_format(axis='y', style='plain', useOffset=False)
plt.title("Varying ppm of CO2 at Summit Air Conditions at 6000 Å")

# plot average condition range on summit
plt.axvline(x=375, color='darkgrey')
plt.axvline(x=500, color='darkgrey')

# graph the original n vals
n_axs_co2 = np.zeros(len(co2_range)) + n_axs[1]
# ax1.plot(co2_range, n_axs_co2, color='g',label=wavelength[1], ls=styles[0], lw=.8)

# graph the Refracive Index vs CO2 for each wavelength
color = 'tab:orange'
ax1.set_xlabel('ppm of CO2')
ax1.set_ylabel('Refractive Index', color=color)
ax1.plot(co2_range, n_co2_2, color=color, label=wavelength[1], ls=styles[0], lw=.8)
ax1.tick_params(axis='y', labelcolor=color)

ax2 = ax1.twinx()

# graph the change in wavelength vs ppm of CO2 for each wavelength
color = 'tab:blue'
ax2.set_ylabel('Change in Wavelength [Å]', color=color)
ax2.plot(co2_range, delta_co2_2, color=color, label=f'{wavelength[1]} Å', ls=(0,(4,4)), lw=.8)
ax2.tick_params(axis='y', labelcolor=color)
ax2.ticklabel_format(useOffset=False)

# Get handles and labels for all data sets
handles1, _ = ax1.get_legend_handles_labels()
handles2, labels2 = ax2.get_legend_handles_labels()

# get the solid color lines from the handles of ax1 and ax2
color_handles = [handles1[0], handles2[0]]

# make new lines for linestyle legend
ls_handle1 = plt.Line2D([0], [0], ls='-', color='black', markersize=10)

# setup legends 
color_legend = ax1.legend(handles=color_handles, labels=[color_labels[1], color_labels[2]], loc='upper center')

linestyle_legend = ax2.legend([ls_handle1], labels2, loc='upper left')

ax1.add_artist(color_legend)
ax2.add_artist(linestyle_legend)

# save plot as png
fig.tight_layout() 
plt.savefig('co2_summit_air_6000A.png')
plt.close()


""" Varying ppm of CO2 Graph at Summit Air Conditions @ 9000 A """

# setup plot
fig, ax1 = plt.subplots()
plt.ticklabel_format(axis='y', style='plain', useOffset=False)
plt.title("Varying ppm of CO2 at Summit Air Conditions at 9000 Å")

# plot average condition range on summit
plt.axvline(x=375, color='darkgrey')
plt.axvline(x=500, color='darkgrey')

# graph the original n vals
n_axs_co2 = np.zeros(len(co2_range)) + n_axs[2]
# ax1.plot(co2_range, n_axs_co2, color='g',label=wavelength[2], ls=styles[0], lw=.8)

# graph the Refracive Index vs CO2 for each wavelength
color = 'tab:orange'
ax1.set_xlabel('ppm of CO2')
ax1.set_ylabel('Refractive Index', color=color)
ax1.plot(co2_range, n_co2_3, color=color, label=wavelength[2], ls=styles[0], lw=.8)
ax1.tick_params(axis='y', labelcolor=color)

ax2 = ax1.twinx()

# graph the change in wavelength vs ppm of CO2 for each wavelength
color = 'tab:blue'
ax2.set_ylabel('Change in Wavelength [Å]', color=color)
ax2.plot(co2_range, delta_co2_3, color=color, label=f'{wavelength[2]} Å', ls=(0,(4,4)), lw=.8)
ax2.tick_params(axis='y', labelcolor=color)
ax2.ticklabel_format(useOffset=False)

# Get handles and labels for all data sets
handles1, _ = ax1.get_legend_handles_labels()
handles2, labels2 = ax2.get_legend_handles_labels()

# get the solid color lines from the handles of ax1 and ax2
color_handles = [handles1[0], handles2[0]]

# make new lines for linestyle legend
ls_handle1 = plt.Line2D([0], [0], ls='-', color='black', markersize=10)

# setup legends 
color_legend = ax1.legend(handles=color_handles, labels=[color_labels[1], color_labels[2]], loc='upper center')

linestyle_legend = ax2.legend([ls_handle1], labels2, loc='upper left')

ax1.add_artist(color_legend)
ax2.add_artist(linestyle_legend)

# save plot as png
fig.tight_layout() 
plt.savefig('co2_summit_air_9000A.png')
plt.close()

print(n_axs)

""" 3D Plots for 2 Variables """


''' Temperature and Pressure vs Change in Wavelength'''

n_axs = 1.0002769832562914
n_ws = 1.00000309

temp_range = np.linspace(265, 285, 20, endpoint=True)
pres_range = np.linspace(50000, 70000, 20, endpoint=True)

X, Y = np.meshgrid(temp_range, pres_range)
Z = n_prop_func(X, Y, .159, 450)
Z = Z * 6000 - 6000

fig = plt.figure(figsize=(12, 6))
ax = plt.axes(projection='3d')

color = ax.plot_surface(X, Y, Z, cmap='viridis', alpha=0.8)

ax.set_title('Effects of Temperature and Pressure on Change in Wavelength at 6000 Å', fontsize=16)
ax.set_xlabel('Temperature [K]', fontsize=12)
ax.set_ylabel('Pressure [Pa]', fontsize=12)
ax.set_zlabel('Change in Wavelength [Å]', fontsize=12)
ax.view_init(azim=-40)

fig.colorbar(color)

plt.tight_layout()
plt.savefig('temp_pres_heatmap.png')


''' Temperature and Humidiy vs Change in Wavelength'''


temp_range = np.linspace(265, 285, 20, endpoint=True)
hum_range = np.linspace(0, 1, 20, endpoint=True)

X, Y = np.meshgrid(temp_range, hum_range)
Z = n_prop_func(X, 61585, Y, 450)
Z = Z * 6000 - 6000

fig = plt.figure(figsize=(12, 6))
ax = plt.axes(projection='3d')

color = ax.plot_surface(X, Y, Z, cmap='plasma', alpha=0.8)

ax.set_title('Effects of Temperature and Humidity on Change in Wavelength at 6000 Å', fontsize= 16)
ax.set_xlabel('Temperature [K]', fontsize=12)
ax.set_ylabel('Fractional Humidity', fontsize=12)
ax.set_zlabel('Change in Wavelength [Å]', fontsize=12)

ax.view_init(azim=-115)
fig.colorbar(color)

plt.tight_layout()
plt.savefig('temp_hum_heatmap.png')


''' Pressure and Humidiy vs Change in Wavelength'''


pres_range = np.linspace(50000, 70000, 20, endpoint=True)
hum_range = np.linspace(0, 1, 20, endpoint=True)


X, Y = np.meshgrid(pres_range, hum_range)
Z = n_prop_func(275.5, X, Y, 450)
Z = Z * 6000 - 6000

fig = plt.figure(figsize=(12, 6))
ax = plt.axes(projection='3d')

color = ax.plot_surface(X, Y, Z, cmap='cool', alpha=0.8)

ax.set_title('Effects of Pressure and Humidity on Change in Wavelength at 6000 Å', fontsize=16)
ax.set_xlabel('Pressure [Pa]', fontsize=12)
ax.set_ylabel('Fractional Humidity', fontsize=12)
ax.set_zlabel('Change in Wavelength [Å]', fontsize=12)

ax.view_init(azim=-115)
fig.colorbar(color)

plt.tight_layout()
plt.savefig('pres_hum_heatmap.png')


''' Pressure and Humidiy vs Change in Wavelength'''


hum_range = np.linspace(0, 1, 20, endpoint=True)
co2_range = np.linspace(0, 800, 20, endpoint=True)

X, Y = np.meshgrid(hum_range, co2_range)
Z = n_prop_func(275.5, 61585, X, Y)
Z = Z * 6000 - 6000

fig = plt.figure(figsize=(12, 6))
ax = plt.axes(projection='3d')

color = ax.plot_surface(X, Y, Z, cmap='magma', alpha=0.8)

ax.set_title('Effects of Humidity and ppm of CO2 on Change in Wavelength at 6000 Å', fontsize=16)
ax.set_xlabel('Fractional Humidity', fontsize=12)
ax.set_ylabel('ppm of CO2', fontsize=12)
ax.set_zlabel('Change in Wavelength [Å]', fontsize=12)

ax.view_init(azim=55)
fig.colorbar(color)

plt.tight_layout()
plt.savefig('hum_co2_heatmap.png')