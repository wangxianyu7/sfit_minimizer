import sfit_minimizer

import MulensModel as mm
import os.path

import time

mm.utils.MAG_ZEROPOINT = 18.

data_path = os.path.join(sfit_minimizer.DATA_PATH, 'MMTest')
datafiles = ['FSPL_Obs_1_I.pho', 'FSPL_Obs_2_V.pho']
datasets = []
for i, filename in enumerate(datafiles):
    if i == 0:
        band = 'I'
    else:
        band = 'V'

    data = mm.MulensData(
        file_name=os.path.join(data_path, filename), phot_fmt='mag',
        bandpass=band)

    datasets.append(data)

gammas = {'I': 0.44,
          'V': 0.72,
          'H': 0.26}
model = mm.Model({'t_0': 2451697.2, 'u_0': 0.00600, 't_E': 25.00000, 'rho': 0.0065})
n_t_star = 10.
t_star = model.parameters.rho * model.parameters.t_E
model.set_magnification_methods([
    model.parameters.t_0 - n_t_star * t_star,
    'finite_source_LD_Yoo04_direct',
    model.parameters.t_0 + n_t_star * t_star])
for band, value in gammas.items():
    model.set_limb_coeff_gamma(band, value)

event = mm.Event(datasets=datasets, model=model)
print(event)
print(model)
parameters_to_fit = ['t_0', 'u_0', 't_E', 'rho']
initial_guess = [2451697.19995,  0.00600, 25.00000,  0.00650, 1.30000,  0.00000,
                 1.00000,  0.00000]

my_func = sfit_minimizer.mm_funcs.PointLensSFitFunction(event, parameters_to_fit)

print('Time to update:')
start = time.time()
my_func.update_all(initial_guess)
stop = time.time()
print(stop - start)

print('Time to fit_fluxes:')
start = time.time()
my_func.event.fit_fluxes()
stop = time.time()
print(stop - start)

# print('Time to calculate residuals:')
# start = time.time()
# my_func.calc_residuals()
# stop = time.time()
# print(stop - start)

print('Time to calculate df:')
start = time.time()
my_func.calc_df()
stop = time.time()
print(stop - start)

for i in range(len(datasets)):
    print('Time to get_magnification for {0}:'.format(i))
    start = time.time()
    event.fits[i].get_data_magnification(bad=False)
    # print(
    #     'Time to initialize FSPLDerivs for {0}:'.format(
    #         i))
    # start = time.time()
    # derivs = mm.FitData.FSPLDerivs(event.fits[i])
    stop = time.time()
    print(stop - start)
    # print('Time to calculate get_gradient for {0}:'.format(i))
    # start = time.time()
    # derivs.get_gradient(parameters_to_fit)
    # stop = time.time()
    # print(stop - start)


# print('Time to calculate bmat:')
# start = time.time()
# my_func.calc_bmat()
# stop = time.time()
# print(stop - start)

# print('Time to calculate cmat:')
# start = time.time()
# my_func.calc_cmat()
# stop = time.time()
# print(stop - start)