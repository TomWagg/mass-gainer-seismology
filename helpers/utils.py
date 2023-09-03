import matplotlib.pyplot as plt
import numpy as np
import astropy.units as u
import astropy.constants as const
import os
import pandas as pd

# matplotlib settings
fs = 24
params = {'figure.figsize': (12, 8),
          'legend.fontsize': 0.6*fs,
          'axes.labelsize': 0.8*fs,
          'axes.titlesize': 0.5*fs,
          'xtick.labelsize': 0.6 * fs,
          'ytick.labelsize': 0.6 * fs,
          'axes.linewidth': 1.1,
          'xtick.major.size': 7,
          'xtick.minor.size': 4,
          'ytick.major.size': 7,
          'ytick.minor.size': 4}

mass_gainer_col = plt.get_cmap("magma")(0.3)
single_col = plt.get_cmap("magma")(0.8)

__all__ = ["fs", "params", "mass_gainer_col", "single_col", "set_styles", "find_closest_model_number",
           "get_pi_0", "asymptotic_period_spacing", "get_radius", "get_eigenfunctions", "get_core_boundary",
           "get_delta_p", "append_surface_He_abundance", "create_GYRE_bash"]

def set_styles():
    plt.rc('font', family='serif')
    plt.rcParams['text.usetex'] = False
    plt.rcParams.update(params)

@np.vectorize
def find_closest_model_number(track, age=None, X_c=None):
    if age is None and X_c is None:
        raise ValueError("At least one of `age` or `X_c` must not be None")

    if X_c is None:
        return np.abs(track.history["star_age"] / 1e6 - age).argmin()
    else:
        return np.abs(track.history["center_h1"] - X_c).argmin()


def get_pi_0(profile):
    p_sort = profile.sort_values(by="logR")

    R = 10**(p_sort["logR"].values)
    r = R / np.max(R)
    N = p_sort["brunt_N"].values
    return 1 / np.trapz(y=N / r, x=r)


def asymptotic_period_spacing(profile, unit=u.day):
    l = 1
    L = np.sqrt(l * (l + 1))
    delta_P = 2 * np.pi**2 / L * get_pi_0(profile=profile) * u.s
    return delta_P.to(unit)   


def get_radius(history=None, L=None, Teff=None):
    if L is None or Teff is None:
        L = 10**(history["log_L"]).values * u.Lsun
        Teff = 10**(history["log_Teff"]).values * u.K
    return np.sqrt(L / (4 * np.pi * const.sigma_sb * Teff**4)).to(u.Rsun)

def get_eigenfunctions(track, profile_number):
    eig_dir = os.path.join(track.dir, 'profile' + str(profile_number) + '-freqs')
    eigs = {}
    for fname in os.listdir(eig_dir):
        if fname[:2] != '00':
            continue
        eig_i = pd.read_table(os.path.join(eig_dir, fname), sep='\s+', skiprows=5)
        sph_deg, rad_ord = fname.split('_')
        eigs[(int(sph_deg), 'M_r')] = eig_i['M_r']
        eigs[(int(sph_deg), 'x')] = eig_i['x']
        eigs[(int(sph_deg), int(rad_ord))] = eig_i[['Re(xi_r)', 'Re(xi_h)', 'dE_dx']]
    return eigs

def get_core_boundary(track, mod):
    p = track.profiles[mod - 1]
    logR_core = p[p["x_mass_fraction_H"] > p["x_mass_fraction_H"].min() * 1.01]["logR"].iloc[-1]
    logR_edge = p["logR"].iloc[0]
    
    return 10**(logR_core) / 10**(logR_edge)


def get_delta_p(track, mod=None, X_c=None, age=None, drop_duplicate_ng=True):
    if age is None and X_c is None and mod is None:
        raise ValueError("At least one of `age` or `X_c` or `mod` must not be None")
    if mod is None:
        mod = find_closest_model_number(track=track, age=age, X_c=X_c)
    df = track.freqs[mod - 1]
    df = df[(df["l"] == 1) & (df["m"] == 0)]# & (df["n_p"] == 0)]
    mixed_modes = df[df["n_p"] > 0]
    duplicates = df[df["n_g"].isin(mixed_modes["n_g"])]
    drop_these = duplicates[duplicates["n_p"] == 0].index

    df = df.drop(index=drop_these)

    if drop_duplicate_ng:
        df = df.drop_duplicates(subset="n_g", keep="last")
            
    periods = 1 / df["Re(freq)"].values * u.day
    ng = df["n_g"].values
    delta_p  = periods[:-1] - periods[1:]
    return periods[:-1], ng[:-1], delta_p


def append_surface_He_abundance(track):
    surface_he = [track.profiles[i].iloc[0]["y_mass_fraction_He"] for i in range(len(track.profiles))]
    track.history.loc[:, "surface_he"] = surface_he


def create_GYRE_bash(mods=None, track=None, X_c=None, procs=6,
                     script="/afs/mpa/temp/tomwagg/kavli/GYRE_submitter.sh", change_folders=False):
    if mods is None:
        mods = find_closest_model_number(track=track, X_c=X_c)
    mods_strings = [f"profile{mod}.data.GYRE" for mod in mods]
    command = "echo -n '" + ','.join(mods_strings) + "' | xargs -d ',' -P "\
        + str(min(procs, len(mods_strings)))\
        + " -I {} " + script + " -i {} -t 1 -e"
    if change_folders:
        command = f"cd {track.dir}; " + command + "; cd -;"
    print(command)