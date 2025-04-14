import numpy as np

def load_lattice_configs(L, Tmin, Tmax, M, isDis, TRank):
    """
    Load lattice configurations from a binary file.

    Args:
        L (int): Lattice size (LxL).
        Tmin (float): Minimum temperature.
        Tmax (float): Maximum temperature.
        M (int): Number of samples.
        isDis (int): Indicator for disorder (1 for disordered, 0 otherwise).
        TRank (int): Rank of the temperature.

    Returns:
        np.ndarray: Array of lattice configurations with shape (-1, L*L).
    """
    # Construct the file path for the lattice configurations
    file_name = "../data/configs/ising_lat_L{:d}_Tmin{:.3f}_Tmax{:.3f}_M{:d}_{:d}/config_rank{:d}.bin".format(L, Tmin, Tmax, M, isDis, TRank)
    dataType = np.int32  # Data type for the binary file
    # Read and reshape the binary data into the desired format
    data = np.fromfile(file_name, dtype=dataType).reshape(-1, L*L)
    return data

def load_lattice_couplings(L, Tmin, Tmax, M, isDis):
    """
    Load lattice couplings (interaction matrix) from a binary file.

    Args:
        L (int): Lattice size (LxL).
        Tmin (float): Minimum temperature.
        Tmax (float): Maximum temperature.
        M (int): Number of samples.
        isDis (int): Indicator for disorder (1 for disordered, 0 otherwise).

    Returns:
        np.ndarray: Array of lattice couplings with shape (-1, L*L).
    """
    # Construct the file path for the lattice couplings
    file_name = "../data/configs/ising_lat_L{:d}_Tmin{:.3f}_Tmax{:.3f}_M{:d}_{:d}/Jmat.bin".format(L, Tmin, Tmax, M, isDis)
    # Use double precision for disordered systems, otherwise use int32
    dataType = np.double if isDis else np.int32
    # Read and reshape the binary data into the desired format
    data = np.fromfile(file_name, dtype=dataType).reshape(-1, L*L)
    return data

def load_rrg_configs(L, Tmin, Tmax, M, isDis, TRank):
    """
    Load random regular graph (RRG) configurations from a binary file.

    Args:
        L (int): Graph size.
        Tmin (float): Minimum temperature.
        Tmax (float): Maximum temperature.
        M (int): Number of samples.
        isDis (int): Indicator for disorder (1 for disordered, 0 otherwise).
        TRank (int): Rank of the temperature.

    Returns:
        np.ndarray: Array of RRG configurations with shape (-1, L*L).
    """
    # Construct the file path for the RRG configurations
    file_name = "../data/configs/ising_rrg_L{:d}_Tmin{:.3f}_Tmax{:.3f}_M{:d}_{:d}/config_rank{:d}.bin".format(L, Tmin, Tmax, M, isDis, TRank)
    dataType = np.int32  # Data type for the binary file
    # Read and reshape the binary data into the desired format
    data = np.fromfile(file_name, dtype=dataType).reshape(-1, L*L)
    return data

def load_rrg_couplings(L, Tmin, Tmax, M, isDis):
    """
    Load random regular graph (RRG) couplings (interaction matrix) from a binary file.

    Args:
        L (int): Graph size.
        Tmin (float): Minimum temperature.
        Tmax (float): Maximum temperature.
        M (int): Number of samples.
        isDis (int): Indicator for disorder (1 for disordered, 0 otherwise).

    Returns:
        np.ndarray: Array of RRG couplings with shape (-1, L*L).
    """
    # Construct the file path for the RRG couplings
    file_name = "../data/configs/ising_rrg_L{:d}_Tmin{:.3f}_Tmax{:.3f}_M{:d}_{:d}/Jmat.bin".format(L, Tmin, Tmax, M, isDis)
    # Use double precision for disordered systems, otherwise use int32
    dataType = np.double if isDis else np.int32
    # Read and reshape the binary data into the desired format
    data = np.fromfile(file_name, dtype=dataType).reshape(-1, L*L)
    return data
