import matplotlib.pyplot as plt
import numpy as np
import astropy.units as u

from utils import find_closest_model_number, mass_gainer_col, single_col, fs,\
    get_eigenfunctions, get_core_boundary

__all__ = ["plot_X_H_profile", "plot_BV_profile", "plot_eigs"]

def plot_X_H_profile(age=None, X_c=None, tracks=None, labels=["Mass-gainer", "Single"],
                     colours=[mass_gainer_col, single_col], fig=None, ax=None, show=True,
                     label_with="title", annotate_col="lightgrey"):
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
    ax.annotate(f"Extends to {m_fin:1.1f} " + r"$\rm M_{\odot}$", xy=(1, 0.5), xytext=(0.96, 0.5),
                xycoords="axes fraction", ha="right", va="center", color=annotate_col,
                arrowprops=dict(arrowstyle="-|>", color=annotate_col), fontsize=0.4*fs,
                bbox=dict(boxstyle="round", fc="white", ec="white", pad=0))

    ax.set_ylim(0, 0.7)
    ax.set_xlim(0, 1.1)
    
    if label_with == "title":
        ax.set_title(r"$X_H$ profile for a final mass " + f"~{m_fin:1.1f} " + r"$\rm M_{\odot}$ star"\
                    + (f' at {age:1.1f} Myr' if X_c is None else r" with $X_c =$" + f' {X_c:1.2f}'))
    elif label_with == "an":
        ax.annotate(f'Age = {age:1.1f} Myr' if X_c is None else r"$X_c =$" + f' {X_c:1.2f}',
                    xy=(0.02, 0.95), xycoords="axes fraction", va="top", fontsize=0.5*fs,
                    bbox=dict(boxstyle="round", pad=0, fc="white", ec="white"))

    plt.tight_layout()

    if show:
        plt.show()
    
    return fig, ax


def plot_BV_profile(age=None, X_c=None, tracks=None, labels=["Mass-gainer", "Single"],
                    colours=[mass_gainer_col, single_col],
                    lw=1, x_scale="linear", fractional_mass=False, radius=False, fractional_radius=False,
                    fill=True, fig=None, ax=None, show=True, label_with="title", legend_loc="upper right"):
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

        if radius:
            m = 10**(track.profiles[mod - 1]["logR"])

        if fractional_radius:
            m = 10**(track.profiles[mod - 1]["logR"]) / 10**(track.profiles[mod - 1]["logR"].max())
        
        brunt = track.profiles[mod - 1]["brunt_N"] * ((2 * np.pi * u.Hz).to(u.day**(-1))).value
        ax.plot(m, brunt, lw=lw, color=col, label=tag, zorder=4)
        if fill:
            ax.fill_between(m, ax.get_ylim()[0], brunt, color=col, alpha=0.2, zorder=4)

    ax.set_ylabel("Brunt–Väisälä\n[Cycles per day]", fontsize=0.5 * fs)

    if radius:
        ax.set_xlabel(r"Radius [$\rm R_{\odot}$]")
    elif fractional_radius:
        ax.set_xlabel(r"Fractional Radius, $r / R$")
    else:
        ax.set_xlabel(r"Mass [$\rm M_{\odot}$]")
    ax.legend(loc=legend_loc, ncol=2, fontsize=0.4 * fs)
    
    m_fin = tracks[0].history["star_mass"].iloc[-1]
    if label_with == "title":
        ax.set_title(r"Profile for a final mass " + f"~{m_fin:1.1f} " + r"$\rm M_{\odot}$ star"\
                        + (f' at {age:1.1f} Myr' if X_c is None else r" with $X_c =$" + f' {X_c:1.2f}'))
    else:
        ax.annotate(f'Age = {age:1.1f} Myr' if X_c is None else r"$X_c =$" + f' {X_c:1.2f}',
                    xy=(0.02, 0.95), xycoords="axes fraction", va="top", fontsize=0.5*fs,
                    bbox=dict(boxstyle="round", pad=0, fc="white", ec="white"), zorder=10)
    
    plt.tight_layout()
    
    if show:
        plt.show()
        
    return fig, ax


def plot_eigs(tracks, sph_deg, rad_ord, which="horizontal", colours=[mass_gainer_col, single_col],
              labels=["Mass-gainer", "Single"], difference=False,
              x_range=None, mod=None, age=None, X_c=None, lw=1,
              plot_core=True, fig=None, ax=None, show=True):
    if age is None and X_c is None and mod is None:
        raise ValueError("At least one of `age` or `X_c` or `mod` must not be None")
    
    if fig is None or ax is None:
        fig, ax = plt.subplots(figsize=(9, 4))

    if isinstance(rad_ord, int):
        rad_ord = [rad_ord]
    if not isinstance(tracks, (list, np.ndarray)):
        tracks = [tracks]
    
    for track, col, label in zip(tracks, colours, labels):
        mod = find_closest_model_number(track=track, age=age, X_c=X_c)
        eigs = get_eigenfunctions(track, mod)
        
        x = eigs[(sph_deg, 'x')]
        x_mask = np.repeat(True, len(x)) if x_range is None else (x > x_range[0]) & (x < x_range[-1])
        
        if difference:
            eig_low = eigs[(sph_deg, rad_ord[0])]
            eig_high = eigs[(sph_deg, rad_ord[1])]
            if which == "enorm":
                ax.plot(x[x_mask], eig_high['dE_dx'][x_mask] - eig_low['dE_dx'][x_mask],
                        lw=lw, zorder=11, color=col, label=label)
            else:
                if which == "radial" or which == "both":
                    ax.plot(x[x_mask], eig_high['Re(xi_r)'][x_mask] - eig_low['Re(xi_r)'][x_mask],
                            color=col, label=label, lw=lw, zorder=11)
                if which == "horizontal" or which == "both": 
                    ax.plot(x[x_mask], eig_high['Re(xi_h)'][x_mask] - eig_low['Re(xi_h)'][x_mask],
                            color=col, label=label, lw=lw, zorder=10)
        else:
            for n in rad_ord:
                eig = eigs[(sph_deg, n)]
                label = f'n = {n}'
                if which == "enorm":
                    ax.plot(x[x_mask], eig['dE_dx'][x_mask], lw=lw, zorder=11, color=col, label=label)
                    # ax.fill_between(x[x_mask], 0, eig['dE_dx'][x_mask], lw=lw, zorder=11, color=col, label=label, alpha=0.4)
                else:
                    if which == "radial" or which == "both":
                        ax.plot(x[x_mask], eig['Re(xi_r)'][x_mask],
                                label=label if which != "both" else label + " (radial)", lw=lw, zorder=11)
                    if which == "horizontal" or which == "both": 
                        ax.plot(x[x_mask], eig['Re(xi_h)'][x_mask],
                                label=label if which != "both" else label + " (horizontal)", lw=lw, zorder=10)
    
    if x_range is None:
        ax.set_xlim([0, 1])
    else:
        ax.set_xlim(x_range)
    # ax.set_ylim([-20, 20])
    
    if plot_core:
        core = get_core_boundary(track, mod)
        ax.axvspan(0, core, ls='dotted', color="lightgrey", alpha=0.5)
        ax.axvline(core, ls='dotted', color="lightgrey")
        ax.annotate("Convective\nCore", xy=(core / 2, ax.get_ylim()[-1] * 0.9), ha="center", va="top", color="grey", rotation=90)
    
    ax.set_xlabel(r'Fractional radius $r/R$')

    if which == "enorm":
        ax.set_ylabel('Differential Inertia', fontsize=0.5*fs)
    else:
        ax.set_ylabel('Normalized\neigenfunction ' + r'$\xi$', fontsize=0.5*fs)
        ax.axhline(0, ls='dotted', color="lightgrey")
    
    X_c = track.history.loc[mod - 1]["center_h1"]
    ax.annotate((f'Age = {age:1.1f} Myr' if X_c is None else rf"$X_c = {{{X_c:1.2f}}}$") + '' if len(rad_ord) > 1 else ('\n' + rf"$(\ell = {{{sph_deg}}}, n = {{{rad_ord[0] if len(rad_ord) == 1 else rad_ord}}})$"),
                xy=(0.5, 0.98), xycoords="axes fraction", ha="center", va="top", fontsize=0.5*fs)
    
    if which == "both" or len(rad_ord) > 1:
        ax.legend(loc='best')

    if show:
        plt.show()
    
    return fig, ax