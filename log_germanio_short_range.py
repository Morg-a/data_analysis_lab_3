# diode I-V semilog plot with fit on short range
import numpy as np
import matplotlib.pyplot as plt
from scipy.odr import ODR, Model, RealData

PM = '\u00B1' # plus-minus symbol

def diode_eq(parameters, x):
    """Definisco la funzione di fit (y = a*(np.exp(x/b) - 1)."""
    y = parameters[0]*(np.exp(x/parameters[1]) - 1)
    return y


def plot_data(x_vector, y_vector, x_err_vector, y_err_vector):
    """Definisco la funzione per plottare i dati in scala semilogaritmica. """
    plt.errorbar(x_vector, y_vector, yerr=y_err_vector, xerr=x_err_vector,
                 fmt='o', color='black')
    plt.ylabel('I (mA)')
    plt.xlabel('V (mV)')
    plt.yscale('log')
    plt.show()


def print_R_squared(y, y_fit):
    """Definisco la funzione per calcolare R^2 e stamparlo."""
    delta = y - y_fit
    R_squared = 1 - (np.var(delta)/np.var(y))
    print("R^2=", R_squared)


#definisco la funzione per fittare i dati

def fit(x_short_vector, y_short_vector, x_err_vector, y_err_vector,
        initial_parameters, x_vector, y_vector):
    """Definisco la funzione per fittare i dati."""
    set_data = RealData(x_short_vector, y_short_vector, sx=x_err_vector,
                        sy=y_err_vector)
    linear_model = Model(diode_eq)
    odr=ODR(set_data, linear_model, beta0=initial_parameters)
    #odr.set_job(fit_type=2)
    fit = odr.run()
    current_I0,eta_Vt = fit.beta
    current_I0_err,eta_Vt_err=fit.sd_beta
    print("I0 = " + str(current_I0) + "      I0 err = " + str(current_I0_err)
          + "\netaVT = " + str(eta_Vt) + "    etaVT err = " + str(eta_Vt_err))
    y_fit=diode_eq(fit.beta, x_vector)

    # calcolo R^2 e lo stampo
    print_R_squared(y_vector, y_fit)

    return y_fit

def plot_fit(x_vector, y_vector, x_err_vector, y_err_vector, y_fit_vector,
             x_max_fit, x_min_fit):
    """Definisco la funzione per plottare la linea di fit sui dati."""
    plt.errorbar(x_vector, y_vector, y_err_vector, x_err_vector, fmt='o',
                 color='black', label='Data')
    plt.plot(x_vector, y_fit_vector, color='red', label='Fit')
    plt.ylabel('I (mA)')
    plt.xlabel('V (mV)')
    plt.yscale('log')
    plt.legend()
    plt.axvline(x=x_max_fit, color='blue', linestyle='--')
    plt.axvline(x=x_min_fit, color='blue', linestyle='--')
    plt.title('Caratteristica I-V del diodo al Germanio')
    plt.savefig('germanio.pdf')
    plt.show()


###################################################################

# prendo da file i dati e li plotto usando la funzione plot_data

voltage_osc, v_osc_err, current, c_err = np.loadtxt('presa_dati_ge.txt',
                                                    delimiter=',',
                                                    unpack=True, comments='#')

voltage_osc_short, v_osc_short_err, current_short, c_short_err = np.loadtxt(
                'dati_fit_ge.txt', delimiter=',', unpack=True, comments='#')

slope = 0.982
slope_err = 0.00658
intercept = -0.03691
intercept_err = 1.5057
ab_cov = -0.00799
voltage = voltage_osc/0.982
voltage_short = voltage_osc_short/0.982
v_err = (1/slope)*np.sqrt((voltage_osc*slope_err/(slope**2))**2 +
                          + v_osc_err**2 + (intercept_err*slope)**2 +
                          + 2*voltage_osc*v_osc_err*ab_cov)
v_short_err = (1/slope)*np.sqrt((voltage_osc_short*slope_err/(slope**2))**2 +
                          + v_osc_short_err**2 + (intercept_err*slope)**2 +
                          + 2*voltage_osc_short*v_osc_short_err*ab_cov)

plot_data(voltage, current, v_err, c_err)

initial_parameters_ge = [0.001, 26]
y_fit = fit(voltage_short, current_short, v_short_err, c_short_err,
            initial_parameters_ge, voltage, current)

# imposta range di fit intorno a 600-800 mV
voltage_max_fit = voltage_short[-1]
voltage_min_fit = voltage_short[0]
plot_fit(voltage, current, v_err, c_err, y_fit, voltage_max_fit,
         voltage_min_fit)

print()
