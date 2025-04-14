import numpy as np
from tqdm import tqdm
import sys

sys.path.append('../src')

from dataReading import *
from inference import *

print("")
print("---------------------------- Ordered lattice ---------------------------")
print("")

L = 8            # System linear size
Tmin = 1.670     # Minimum temperature
Tmax = 30.000    # Maximum temperature
M = 40           # Number of temperatures
isDis = False    # Is disorder present?
nSamples = 15000 # Number of samples per temperature

betaArray = np.linspace(1/Tmax, 1/Tmin, M)
TArray = 1/betaArray

latticeOrd_MF_gammaArray = np.zeros(M)
latticeOrd_PL_gammaArray = np.zeros(M)
latticeOrd_PLWithLasso_gammaArray = np.zeros(M)

for i in range(M):
    print("")
    print(f"Processing temperature {i+1}/{M}...")
    print("")

    beta = betaArray[i]
    configs = load_lattice_configs(L, Tmin, Tmax, M, isDis, i)
    selected_configs = configs[np.random.choice(configs.shape[0], nSamples, replace=False)]
    
    trueCouplings = load_lattice_couplings(L, Tmin, Tmax, M, isDis)
    true_gamma = gamma(trueCouplings, trueCouplings)

    MF_inferredCouplings = evaluate_MF_inference(selected_configs, beta)
    MF_gamma = gamma(trueCouplings, MF_inferredCouplings)
    
    PL_inferredCouplings = evaluate_PL_inference(selected_configs, beta, learning_rate=10., epochs=5000, decay=0.999, l1_amplitude=0.0)
    PL_gamma = gamma(trueCouplings, PL_inferredCouplings)

    PLWithLasso_inferredCouplings = evaluate_PL_inference(selected_configs, beta, learning_rate=10., epochs=5000, decay=0.999, l1_amplitude=0.1)
    PLWithLasso_gamma = gamma(trueCouplings, PLWithLasso_inferredCouplings)

    latticeOrd_MF_gammaArray[i] = MF_gamma
    latticeOrd_PL_gammaArray[i] = PL_gamma
    latticeOrd_PLWithLasso_gammaArray[i] = PLWithLasso_gamma

# Save the ordered lattice results to a CSV file
np.savetxt(f'../data/inference/latticeOrd_gammaArray_nSamples{nSamples}.csv', 
           np.column_stack((betaArray, latticeOrd_MF_gammaArray, latticeOrd_PL_gammaArray, latticeOrd_PLWithLasso_gammaArray)), 
           delimiter=',', 
           fmt='%.3f', 
           header='beta,MF_gamma,PL_gamma,PLWithLasso_gamma', 
           comments='')

print("")
print("---------------------------- Ordered rrg ---------------------------")
print("")

L = 8            # System linear size
Tmin = 1.670     # Minimum temperature
Tmax = 30.000    # Maximum temperature
M = 40           # Number of temperatures
isDis = False    # Is disorder present?
nSamples = 15000 # Number of samples per temperature

betaArray = np.linspace(1/Tmax, 1/Tmin, M)
TArray = 1/betaArray

rrgOrd_MF_gammaArray = np.zeros(M)
rrgOrd_PL_gammaArray = np.zeros(M)
rrgOrd_PLWithLasso_gammaArray = np.zeros(M)

for i in range(M):

    print("")
    print(f"Processing temperature {i+1}/{M}...")
    print("")

    beta = betaArray[i]
    configs = load_rrg_configs(L, Tmin, Tmax, M, isDis, i)
    selected_configs = configs[np.random.choice(configs.shape[0], nSamples, replace=False)]
    
    trueCouplings = load_rrg_couplings(L, Tmin, Tmax, M, isDis)
    true_gamma = gamma(trueCouplings, trueCouplings)

    MF_inferredCouplings = evaluate_MF_inference(selected_configs, beta)
    MF_gamma = gamma(trueCouplings, MF_inferredCouplings)
    
    PL_inferredCouplings = evaluate_PL_inference(selected_configs, beta, learning_rate=10., epochs=5000, decay=0.999, l1_amplitude=0.0)
    PL_gamma = gamma(trueCouplings, PL_inferredCouplings)

    PLWithLasso_inferredCouplings = evaluate_PL_inference(selected_configs, beta, learning_rate=10., epochs=5000, decay=0.999, l1_amplitude=0.1)
    PLWithLasso_gamma = gamma(trueCouplings, PLWithLasso_inferredCouplings)

    rrgOrd_MF_gammaArray[i] = MF_gamma
    rrgOrd_PL_gammaArray[i] = PL_gamma
    rrgOrd_PLWithLasso_gammaArray[i] = PLWithLasso_gamma

# Save the ordered rrg results to a CSV file
np.savetxt(f'../data/inference/rrgOrd_gammaArray_nSamples{nSamples}.csv', 
           np.column_stack((betaArray, rrgOrd_MF_gammaArray, rrgOrd_PL_gammaArray, rrgOrd_PLWithLasso_gammaArray)), 
           delimiter=',', 
           fmt='%.3f', 
           header='beta,MF_gamma,PL_gamma,PLWithLasso_gamma', 
           comments='')

print("")
print("---------------------------- Disordered lattice ---------------------------")
print("")

L = 8            # System linear size
Tmin = 1.670     # Minimum temperature
Tmax = 30.000    # Maximum temperature
M = 40           # Number of temperatures
isDis = True     # Is disorder present? (Disordered version)
nSamples = 15000 # Number of samples per temperature

betaArray = np.linspace(1/Tmax, 1/Tmin, M)
TArray = 1/betaArray

latticeDis_MF_gammaArray = np.zeros(M)
latticeDis_PL_gammaArray = np.zeros(M)
latticeDis_PLWithLasso_gammaArray = np.zeros(M)

for i in range(M):

    print("")
    print(f"Processing temperature {i+1}/{M}...")
    print("")

    beta = betaArray[i]
    configs = load_lattice_configs(L, Tmin, Tmax, M, isDis, i)
    selected_configs = configs[np.random.choice(configs.shape[0], nSamples, replace=False)]
    
    trueCouplings = load_lattice_couplings(L, Tmin, Tmax, M, isDis)
    true_gamma = gamma(trueCouplings, trueCouplings)

    MF_inferredCouplings = evaluate_MF_inference(selected_configs, beta)
    MF_gamma = gamma(trueCouplings, MF_inferredCouplings)
    
    PL_inferredCouplings = evaluate_PL_inference(selected_configs, beta, learning_rate=10., epochs=5000, decay=0.999, l1_amplitude=0.0)
    PL_gamma = gamma(trueCouplings, PL_inferredCouplings)

    PLWithLasso_inferredCouplings = evaluate_PL_inference(selected_configs, beta, learning_rate=10., epochs=5000, decay=0.999, l1_amplitude=0.1)
    PLWithLasso_gamma = gamma(trueCouplings, PLWithLasso_inferredCouplings)

    latticeDis_MF_gammaArray[i] = MF_gamma
    latticeDis_PL_gammaArray[i] = PL_gamma
    latticeDis_PLWithLasso_gammaArray[i] = PLWithLasso_gamma

# Save the disordered lattice results to a CSV file
np.savetxt(f'../data/inference/latticeDis_gammaArray_nSamples{nSamples}.csv', 
           np.column_stack((betaArray, latticeDis_MF_gammaArray, latticeDis_PL_gammaArray, latticeDis_PLWithLasso_gammaArray)), 
           delimiter=',', 
           fmt='%.3f', 
           header='beta,MF_gamma,PL_gamma,PLWithLasso_gamma', 
           comments='')

print("")
print("---------------------------- Disordered rrg ---------------------------")
print("")

L = 8            # System linear size
Tmin = 1.670     # Minimum temperature
Tmax = 30.000    # Maximum temperature
M = 40           # Number of temperatures
isDis = True     # Is disorder present? (Disordered version)
nSamples = 15000 # Number of samples per temperature

betaArray = np.linspace(1/Tmax, 1/Tmin, M)
TArray = 1/betaArray

rrgDis_MF_gammaArray = np.zeros(M)
rrgDis_PL_gammaArray = np.zeros(M)
rrgDis_PLWithLasso_gammaArray = np.zeros(M)

for i in range(M):

    print("")
    print(f"Processing temperature {i+1}/{M}...")
    print("")

    beta = betaArray[i]
    configs = load_rrg_configs(L, Tmin, Tmax, M, isDis, i)
    selected_configs = configs[np.random.choice(configs.shape[0], nSamples, replace=False)]
    
    trueCouplings = load_rrg_couplings(L, Tmin, Tmax, M, isDis)
    true_gamma = gamma(trueCouplings, trueCouplings)

    MF_inferredCouplings = evaluate_MF_inference(selected_configs, beta)
    MF_gamma = gamma(trueCouplings, MF_inferredCouplings)
    
    PL_inferredCouplings = evaluate_PL_inference(selected_configs, beta, learning_rate=10., epochs=5000, decay=0.999, l1_amplitude=0.0)
    PL_gamma = gamma(trueCouplings, PL_inferredCouplings)

    PLWithLasso_inferredCouplings = evaluate_PL_inference(selected_configs, beta, learning_rate=10., epochs=5000, decay=0.999, l1_amplitude=0.1)
    PLWithLasso_gamma = gamma(trueCouplings, PLWithLasso_inferredCouplings)

    rrgDis_MF_gammaArray[i] = MF_gamma
    rrgDis_PL_gammaArray[i] = PL_gamma
    rrgDis_PLWithLasso_gammaArray[i] = PLWithLasso_gamma

# Save the disordered rrg results to a CSV file
np.savetxt(f'../data/inference/rrgDis_gammaArray_nSamples{nSamples}.csv', 
           np.column_stack((betaArray, rrgDis_MF_gammaArray, rrgDis_PL_gammaArray, rrgDis_PLWithLasso_gammaArray)), 
           delimiter=',', 
           fmt='%.3f', 
           header='beta,MF_gamma,PL_gamma,PLWithLasso_gamma', 
           comments='')