import matplotlib.pyplot as plt
import numpy as np
import astropy.units as u
import astropy.constants as const
import os
import pandas as pd

plt.rc('font', family='serif')
plt.rcParams['text.usetex'] = False
fs = 24

# update various fontsizes to match
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
plt.rcParams.update(params)

mass_gainer_col = plt.get_cmap("magma")(0.3)
single_col = plt.get_cmap("magma")(0.8)


"""First some help functions that I use in plotting"""

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
        eigs[(int(sph_deg), 'x')] = eig_i['x']
        eigs[(int(sph_deg), int(rad_ord))] = eig_i[['Re(xi_r)', 'Re(xi_h)']]
    return eigs

def get_core_boundary(track, mod):
    p = track.profiles[mod - 1]
    logR_core = p[p["x_mass_fraction_H"] > p["x_mass_fraction_H"].min() * 1.01]["logR"].iloc[-1]
    logR_edge = p["logR"].iloc[0]
    
    return 10**(logR_core) / 10**(logR_edge)

"""Then the actual plotting functions"""


def simple_hr(track=None, df=None, ylabel=r'Luminosity $\log_{10}(\mathbf{L/L_{\odot}})$',
              cbar_var="center_he4", cbar_label=r"$X_{\rm He, center}$", trim_pre_ms=True,
              fig=None, ax=None, show=True, add_axes_info=False, plot_line=True, 
              cbar_loc=[0.38, 0.025, 0.6, 0.025], annotate_start=None, annotate_end=None, R_levels=None,
              **kwargs):
    new_fig = (fig is None or ax is None)
    if new_fig:
        fig, ax = plt.subplots(figsize=(8,6))
    
    if df is None:
        df = track.history

    if trim_pre_ms:
        df = df.loc[df.center_h1 <= df.center_h1.max() - 0.005]

    if 's' not in kwargs:
        kwargs['s'] = 10
    if 'cmap' not in kwargs and cbar_var is not None:
        kwargs['cmap'] = "copper"
    c = df[cbar_var] if cbar_var is not None else None
    
    if plot_line:
        ax.plot(df['log_Teff'], df['log_L'], color="lightgrey", zorder=-1)
    ax.scatter(df['log_Teff'], df['log_L'], c=c, **kwargs)

    if annotate_start is not None:
        ax.annotate(annotate_start, xy=(df['log_Teff'].iloc[0], df['log_L'].iloc[0]), color="lightgrey",
                    ha="right", va="top", fontsize=0.5*fs)
    if annotate_end is not None:
        ax.annotate(annotate_end, xy=(df['log_Teff'].iloc[-1], df['log_L'].iloc[-1]), color="lightgrey",
                    ha="left", va="bottom", fontsize=0.5*fs)
    
    if new_fig or add_axes_info:
        ax.invert_xaxis()
        ax.set_xlabel(r'Effective temperature $\log_{10}(\mathbf{T_{eff}/K})$')
        ax.set_ylabel(ylabel)

        if cbar_var is not None:
            inset_ax = ax.inset_axes(cbar_loc)
            fig.colorbar(ax.collections[0], label=cbar_label, cax=inset_ax, orientation="horizontal", location="top")
    
    if R_levels is not None:
        T_range = np.logspace(*ax.get_xlim(), 1000)
        L_range = np.logspace(*ax.get_ylim(), 1000)

        T, L = np.meshgrid(T_range, L_range)

        R = get_radius(Teff=T * u.K, L=L * u.Lsun)

        def fmt(x):
            return f"{x:1.2f} " + r"$\rm R_{\odot}$"

        CS = ax.contour(np.log10(T), np.log10(L), R.to(u.Rsun).value, levels=R_levels,
                        colors="lightgrey", zorder=-1, linewidths=1, linestyles="dotted", alpha=0.4)
        ax.clabel(CS, CS.levels, inline=True, fontsize=10, fmt=fmt)

    if show:
        plt.show()
        
    return fig, ax


def plot_X_H_profile(age=None, X_c=None, tracks=None, labels=["Mass-gainer", "Single"],
                    colours=[mass_gainer_col, single_col], fig=None, ax=None, show=True,
                     label_with="title", fill=False):
    if age is None and X_c is None:
        raise ValueError("At least one of `age` or `X_c` must not be None")
    if tracks is None:
        raise ValueError("`tracks` cannot be None")
    if fig is None or ax is None:
        fig, ax = plt.subplots(figsize=(7, 3))
    plt.cla()

    lw = 2

    for track, label, col in zip(tracks, labels, colours):
        mod = find_closest_model_number(track=track, age=age, X_c=X_c)
        ax.plot(track.profiles[mod]["mass"],
                track.profiles[mod]["x_mass_fraction_H"],
                label=label, zorder=3, color=col, lw=lw)

    ax.legend(loc="lower right", fontsize=0.4 * fs)

    ax.set_xlabel(r"Mass [$\rm M_{\odot}$]")
    ax.set_ylabel("H mass fraction", fontsize=0.7*fs)

    m_fin = tracks[0].history["star_mass"].iloc[-1]
    ax.annotate(f"Extends to {m_fin:1.1f} " + r"$\rm M_{\odot}$", xy=(1, 0.5), xytext=(0.95, 0.5), xycoords="axes fraction",
                 ha="right", va="center", color="lightgrey", arrowprops=dict(arrowstyle="-|>", color="lightgrey"))

    ax.set_ylim(0, 0.7)
    ax.set_xlim(0, 1.1)
    
    if label_with == "title":
        ax.set_title(r"$X_H$ profile for a final mass " + f"~{m_fin:1.1f} " + r"$\rm M_{\odot}$ star"\
                    + (f' at {age:1.1f} Myr' if X_c is None else r" with $X_c =$" + f' {X_c:1.2f}'))
    else:
        ax.annotate(f'Age = {age:1.1f} Myr' if X_c is None else r"$X_c =$" + f' {X_c:1.2f}',
                    xy=(0.02, 0.95), xycoords="axes fraction", va="top", fontsize=0.5*fs)

    plt.tight_layout()

    if show:
        plt.show()
    
    return fig, ax


def plot_BV_profile(age=None, X_c=None, tracks=None, labels=["Mass-gainer", "Single"],
                    colours=[mass_gainer_col, single_col],
                    lw=1, x_scale="linear", fractional_mass=False, fill=True,
                    fig=None, ax=None, show=True, label_with="title", legend_loc="upper right"):
    if age is None and X_c is None:
        raise ValueError("At least one of `age` or `X_c` must not be None")
    if tracks is None:
        raise ValueError("`tracks` cannot be None")
    if fig is None or ax is None:
        fig, ax = plt.subplots(figsize=(7, 3))
    plt.cla()

    ax.set_yscale("log")
    ax.set_xscale(x_scale)
    ax.set_ylim(5e1, 1e4)
    
    for track, tag, col in zip(tracks, labels, colours):
        mod = find_closest_model_number(track=track, age=age, X_c=X_c)
        m = track.profiles[mod - 1]["mass"]
        if fractional_mass:
            m = m / m.max()
        
        brunt = track.profiles[mod - 1]["brunt_N"] * ((2 * np.pi * u.Hz).to(u.day**(-1))).value
        ax.plot(m, brunt, lw=lw, color=col, label=tag, zorder=4)
        if fill:
            ax.fill_between(m, ax.get_ylim()[0], brunt, color=col, alpha=0.2, zorder=4)

    ax.set_ylabel("Brunt–Väisälä\n[Cycles per day]", fontsize=0.5 * fs)
    ax.set_xlabel(r"Mass [$\rm M_{\odot}$]")
    ax.legend(loc=legend_loc, ncol=2, fontsize=0.4 * fs)
    
    m_fin = tracks[0].history["star_mass"].iloc[-1]
    if label_with == "title":
        ax.set_title(r"Profile for a final mass " + f"~{m_fin:1.1f} " + r"$\rm M_{\odot}$ star"\
                        + (f' at {age:1.1f} Myr' if X_c is None else r" with $X_c =$" + f' {X_c:1.2f}'))
    else:
        ax.annotate(f'Age = {age:1.1f} Myr' if X_c is None else r"$X_c =$" + f' {X_c:1.2f}',
                    xy=(0.02, 0.95), xycoords="axes fraction", va="top", fontsize=0.5*fs)
    
    plt.tight_layout()
    
    if show:
        plt.show()
        
    return fig, ax


def plot_period_spacing(age=None, X_c=None, tracks=None, labels=["Mass-gainer", "Single"],
                        colours=[mass_gainer_col, single_col], legend_loc="upper left", label_with="an",
                        x_var="period", label_modes=False, xlims=None, ylims=None, divide_delta_n=False,
                        fig=None, ax=None, show=True):
    if age is None and X_c is None:
        raise ValueError("At least one of `age` or `X_c` must not be None")
    
    if fig is None or ax is None:
        fig, ax = plt.subplots(figsize=(7, 3))
    plt.cla()
    
    for track, tag, col in zip(tracks, labels, colours):
        
        mod = find_closest_model_number(track=track, age=age, X_c=X_c)

        
        df = track.freqs[mod - 1]
        df = df[(df["l"] == 1) & (df["m"] == 0)]# & (df["n_p"] == 0)]
        mixed_modes = df[df["n_p"] > 0]
        duplicates = df[df["n_g"].isin(mixed_modes["n_g"])]
        drop_these = duplicates[duplicates["n_p"] == 0].index

        df = df.drop(index=drop_these)
                
        periods = 1 / df["Re(freq)"].values * u.day
        ng = df["n_g"].values
        delta_p  = periods[:-1] - periods[1:]

        if divide_delta_n:
            delta_n = ng[:-1] - ng[1:]
            delta_n[delta_n == 0] = 1
            delta_p /= delta_n

        x_vals = periods if x_var == "period" else -ng
    
        if label_modes:
            ax.plot(x_vals[:-1], delta_p, label=tag, color=col)
        else:
            ax.plot(x_vals[:-1], delta_p, marker="o", label=tag, color=col)

        if label_modes:
            for i in range(len(delta_p)):
                ax.annotate(ng[i], xy=(x_vals[i], delta_p[i]), ha="center", va="center",
                            fontsize=0.25*fs, color="grey",
                            bbox=dict(boxstyle="circle", facecolor="white", ec=col))

        aps = asymptotic_period_spacing(track.profiles[mod - 1]).to(u.day).value
        ax.axhline(aps, color=col, alpha=0.5, linestyle="dotted", zorder=-1)#, label=f"Expected asymptotic spacing ({tag})")

    ax.set_xlabel(r"Period, $P \, [\rm days]$" if x_var == "period" else r"$k$ (g-modes)")
    ax.set_ylabel(r"$\Delta P \, [\rm days]$")

    # ax.set_title(f"\nAge = {age:1.1f} Myr, X_c = {X_c:1.2f}", fontsize=0.6*fs)
    
    if label_with == "title":
        ax.set_title(r"Period spacing for final mass ~3.5 $\rm M_{\odot}$ star"\
                     + (f' at {age:1.1f} Myr' if X_c is None else r" with $X_c =$" + f' {X_c:1.2f}'))
    else:
        ax.annotate(f'Age = {age:1.1f} Myr' if X_c is None else r"$X_c =$" + f' {X_c:1.2f}',
                    xy=(0.02, 0.95), xycoords="axes fraction", va="top", fontsize=0.5*fs)

    ax.annotate(r"$(l = 1, m = 0)$ g modes", xy=(0.02, 0.02), xycoords="axes fraction", va="bottom")

    if ylims is not None:
        if ylims == "auto":
            ax.set_ylim(aps - 0.04, aps + 0.04)
        else:
            ax.set_ylim(ylims)

    if xlims is not None:
        ax.set_xlim(xlims)

    ax.legend(loc=legend_loc, fontsize=0.5 * fs)
    
    if show:
        plt.show()

    return fig, ax, ng, periods


def plot_eigs(track, sph_deg, rad_ord, mod=None, age=None, X_c=None, show=True):
    if age is None and X_c is None and mod is None:
        raise ValueError("At least one of `age` or `X_c` or `mod` must not be None")
    
    if mod is None:
        mod = find_closest_model_number(track=track, age=age, X_c=X_c)
    
    fig, ax = plt.subplots(figsize=(9, 4))
    
    eigs = get_eigenfunctions(track, mod)
    
    x   = eigs[(sph_deg, 'x')]
    eig = eigs[(sph_deg, rad_ord)]
    
    ax.plot(x, eig['Re(xi_r)'], label='radial',     lw=2, zorder=11)
    ax.plot(x, eig['Re(xi_h)'], label='horizontal', lw=2, zorder=10)
    
    ax.set_xlim([0, 1])
    ax.set_ylim([-20, 20])
    
    core = get_core_boundary(track, mod)
    ax.axvspan(0, core, ls='dotted', color="lightgrey", alpha=0.5)
    ax.axvline(core, ls='dotted', color="lightgrey")
    ax.annotate("Convective\nCore", xy=(core / 2, ax.get_ylim()[-1] * 0.9), ha="center", va="top", color="grey")
    
    ax.set_xlabel(r'Fractional radius $r/R$')
    ax.set_ylabel('Normalized\neigenfunction ' + r'$\xi$')
    
    ax.axhline(0, ls='dotted', color="lightgrey")
    
    X_c = track.history.loc[mod - 1]["center_h1"]
    ax.set_title(r'Eigenfunctions at $X_c=$' + f'{X_c:1.3f} ' + r'($\ell = $' + f'{sph_deg}, ' + r'$n = $' + f'{rad_ord})')
    ax.legend(loc='lower left')

    if show:
        plt.show()
    
    return fig, ax