import matplotlib.pyplot as plt
import numpy as np
import astropy.units as u

plt.rc('font', family='serif')
plt.rcParams['text.usetex'] = False
fs = 24

# update various fontsizes to match
params = {'figure.figsize': (12, 8),
          'legend.fontsize': 0.6*fs,
          'axes.labelsize': 0.8*fs,
          'xtick.labelsize': 0.6 * fs,
          'ytick.labelsize': 0.6 * fs,
          'axes.linewidth': 1.1,
          'xtick.major.size': 7,
          'xtick.minor.size': 4,
          'ytick.major.size': 7,
          'ytick.minor.size': 4}
plt.rcParams.update(params)


def simple_hr(track, ylabel=r'Luminosity $\mathbf{\log_{10}(L/L_{\odot})}$',
              cbar_var="center_he4", cbar_label=r"$X_{\rm He, center}$", trim_pre_ms=True,
              fig=None, ax=None, show=True, add_axes_info=False, plot_line=True, **kwargs):
    new_fig = (fig is None or ax is None)
    if new_fig:
        fig, ax = plt.subplots(figsize=(8,6))
    
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
    
    if new_fig or add_axes_info:
        ax.invert_xaxis()
        ax.set_xlabel(r'Effective temperature $\mathbf{\log_{10}(T_{eff}/K)}$')
        ax.set_ylabel(ylabel)

        if cbar_var is not None:
            inset_ax = ax.inset_axes([0.38, 0.025, 0.6, 0.025])
            fig.colorbar(ax.collections[0], label=cbar_label, cax=inset_ax, orientation="horizontal", location="top")
    
    if show:
        plt.show()
        
    return fig, ax


def plot_X_H_profile(age=None, X_c=None, tracks=None, mt_index=1, ref_index=2, fig=None, ax=None, show=True):
    if age is None and X_c is None:
        raise ValueError("At least one of `age` or `X_c` must not be None")
    if tracks is None:
        raise ValueError("`tracks` cannot be None")
    if fig is None or ax is None:
        fig, ax = plt.subplots(figsize=(7, 4))
    plt.cla()

    if X_c is None:
        acc_mod = np.abs(tracks[mt_index].history["star_age"] / 1e6 - age).argmin()
        non_acc_mod = np.abs(tracks[ref_index].history["star_age"] / 1e6 - age).argmin()
    else:
        acc_mod = np.abs(tracks[mt_index].history["center_h1"] - X_c).argmin()
        non_acc_mod = np.abs(tracks[ref_index].history["center_h1"] - X_c).argmin()

    ax.plot(tracks[mt_index].profiles[acc_mod]["mass"],
            tracks[mt_index].profiles[acc_mod]["x_mass_fraction_H"],
             label="Had accretion", zorder=3)
    ax.plot(tracks[ref_index].profiles[non_acc_mod]["mass"],
            tracks[ref_index].profiles[non_acc_mod]["x_mass_fraction_H"],
             label="No accretion")

    ax.legend(loc="lower right")

    ax.set_xlabel(r"Mass [$\rm M_{\odot}$]")
    ax.set_ylabel("H mass fraction")

    m_fin = tracks[mt_index].history["star_mass"].iloc[-1]
    ax.annotate(f"Extends to {m_fin:1.1f} " + r"$\rm M_{\odot}$", xy=(1, 0.5), xytext=(0.95, 0.5), xycoords="axes fraction",
                 ha="right", va="center", color="lightgrey", arrowprops=dict(arrowstyle="-|>", color="lightgrey"))

    ax.set_ylim(0, 0.7)
    ax.set_xlim(0, 1.1)
    
    ax.set_title(r"$X_H$ profile for a final mass " + f"~{m_fin:1.1f} " + r"$\rm M_{\odot}$ star"\
                 + (f' at {age:1.1f} Myr' if X_c is None else r" with $X_c =$" + f' {X_c:1.2f}'))

    plt.tight_layout()

    if show:
        plt.show()
    
    return fig, ax


def plot_BV_profile(age=None, X_c=None, tracks=None, mt_index=3, ref_index=2,
                    lw=1, x_scale="linear", fractional_mass=False, fill=True,
                    fig=None, ax=None, show=True):
    if age is None and X_c is None:
        raise ValueError("At least one of `age` or `X_c` must not be None")
    if tracks is None:
        raise ValueError("`tracks` cannot be None")
    if fig is None or ax is None:
        fig, ax = plt.subplots(figsize=(7, 3))
    plt.cla()
    

    if X_c is None:
        print(f"Plotting BV profile for age = {age:1.1f} Myr")
        acc_mod = np.abs(tracks[mt_index].history["star_age"] / 1e6 - age).argmin()
        non_acc_mod = np.abs(tracks[ref_index].history["star_age"] / 1e6 - age).argmin()
    else:
        print(f"Plotting BV profile for X_c = {X_c:1.1f}")
        acc_mod = np.abs(tracks[mt_index].history["center_h1"] - X_c).argmin()
        non_acc_mod = np.abs(tracks[ref_index].history["center_h1"] - X_c).argmin()

    ax.set_yscale("log")
    ax.set_xscale(x_scale)
    ax.set_ylim(1e0, 1e6)
    
    for mod, track, tag, col in zip([acc_mod, non_acc_mod],
                                    [tracks[mt_index], tracks[ref_index]],
                                    ["Mass-gainer", "Single"],
                                    [plt.get_cmap("magma")(0.3), plt.get_cmap("magma")(0.8)]):
        m = track.profiles[mod - 1]["mass"]
        if fractional_mass:
            m = m / m.max()
        
        brunt = track.profiles[mod - 1]["brunt_N"] * ((2 * np.pi * u.Hz).to(u.day**(-1))).value
        ax.plot(m, brunt, lw=lw, color=col, label=tag, zorder=4)
        if fill:
            ax.fill_between(m, ax.get_ylim()[0], brunt, color=col, alpha=0.2, zorder=4)

    ax.set_ylabel("Brunt–Väisälä\n[Cycles per day]", fontsize=0.5 * fs)
    ax.set_xlabel(r"Mass $\mathbf{M/M_{\odot}}$")
    ax.legend(loc="best", ncol=2, fontsize=0.4 * fs)
    
    m_fin = tracks[mt_index].history["star_mass"].iloc[-1]
    ax.set_title(r"Profile for a final mass " + f"~{m_fin:1.1f} " + r"$\rm M_{\odot}$ star"\
                     + (f' at {age:1.1f} Myr' if X_c is None else r" with $X_c =$" + f' {X_c:1.2f}'))    
    
    plt.tight_layout()
    
    if show:
        plt.show()
        
    return fig, ax