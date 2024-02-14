import numpy as np
import matplotlib.pyplot as plt
import astropy.units as u

from utils import fs, get_radius, find_closest_model_number

__all__ = ["simple_hr", "add_singles_tracks", "get_hr_position"]


def get_hr_position(track, pos):
    df = track.history
    if pos == "zams":
        df = df.loc[df.center_h1 <= df.center_h1.max() - 0.005]
        return df.iloc[0]["log_Teff"], df.iloc[0]["log_L"]
    elif pos == "mt_start":
        df = df.loc[df["log_abs_mdot"] > -8]
        return df.iloc[0]["log_Teff"], df.iloc[0]["log_L"]
    elif pos == "closest_approach":
        df = df.loc[df["log_abs_mdot"] > -8]
        df = df.loc[df["log_L"] == df["log_L"].min()]
        return df.iloc[0]["log_Teff"], df.iloc[0]["log_L"]
    elif pos == "mt_end":
        df = df.loc[df["log_abs_mdot"] > -8]
        return df.iloc[-1]["log_Teff"], df.iloc[-1]["log_L"]
    elif pos == "cheb":
        mt_end_ind = df.loc[df["log_abs_mdot"] > -8].index[-1] + 1
        df = df.iloc[mt_end_ind:]
        df = df.loc[df.log_L == df.log_L.min()]
        return df.iloc[0]["log_Teff"], df.iloc[0]["log_L"]
    else:
        raise ValueError("Invalid pos")

def add_singles_tracks(fig, ax, tracks, Ms=None, colour="lightgrey", linestyle="-", an_every=0.5):
    if Ms is None:
        Ms = np.arange(3.0, 6.1, 0.1)
    for M in Ms:
        folder = f"M_{M:1.1f}"
        ans = fr"$M = {{{M:1.1f}}} {{\rm M_\odot}}$" if M.round(1) % an_every == 0.0 else None
        fig, ax = simple_hr(tracks[folder], fig=fig, ax=ax, show=False, cbar_var=None,
                            annotate_start=ans, s=0, line_colour=colour, rasterized=True,
                            linestyle=linestyle, line_z=-10)
    return fig, ax

def simple_hr(track=None, df=None, ylabel=r'Luminosity $\log_{10}(L/{\rm L_{\odot}})$',
              cbar_var="center_he4", cbar_label=r"$X_{\rm He, center}$", trim_pre_ms=True,
              fig=None, ax=None, show=True, add_axes_info=False,
              plot_line=True, line_colour="grey", line_z=-1, linestyle="-",
              cbar_loc=[0.38, 0.025, 0.6, 0.025], inset_cbar=True,
              annotate_start=None, annotate_end=None, R_levels=None, mod_range=None, time_step=None,
              **kwargs):
    new_fig = (fig is None or ax is None)
    if new_fig:
        fig, ax = plt.subplots(figsize=(8,6))

    if df is None:
        df = track.history

    if trim_pre_ms:
        df = df.loc[df.center_h1 <= df.center_h1.max() - 0.005]

    if mod_range is not None:
        df = df.loc[(df.index >= mod_range[0]) & (df.index < mod_range[1])]

    if time_step is not None:
        mods = find_closest_model_number(track, age=np.arange(df["star_age"].min() / 1e6, df["star_age"].max() / 1e6, time_step))
        df = df.loc[mods]


    if 's' not in kwargs:
        kwargs['s'] = 10
    if 'cmap' not in kwargs and cbar_var is not None:
        kwargs['cmap'] = "copper"
    c = df[cbar_var] if cbar_var is not None else None
    
    if plot_line:
        ax.plot(df['log_Teff'], df['log_L'], color=line_colour, zorder=line_z, linestyle=linestyle)
    ax.scatter(df['log_Teff'], df['log_L'], c=c, **kwargs)

    if annotate_start is not None:
        ax.annotate(annotate_start, xy=(df['log_Teff'].iloc[0], df['log_L'].iloc[0]), color=line_colour,
                    ha="right", va="top", fontsize=0.5*fs)
    if annotate_end is not None:
        ax.annotate(annotate_end, xy=(df['log_Teff'].iloc[-1], df['log_L'].iloc[-1]), color=line_colour,
                    ha="left", va="bottom", fontsize=0.5*fs)
    
    if new_fig or add_axes_info:
        ax.invert_xaxis()
        ax.set_xlabel(r'Effective temperature $\log_{10}(T_{\rm eff}/{\rm K})$')
        ax.set_ylabel(ylabel)

        if cbar_var is not None:
            if inset_cbar:
                inset_ax = ax.inset_axes(cbar_loc)
                cb = fig.colorbar(ax.collections[0], cax=inset_ax,
                                  orientation="horizontal", location="top", extend="min")
                cb.set_label(label=cbar_label, fontsize=0.6*fs)
            else:
                fig.colorbar(ax.collections[0], label=cbar_label)
    
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
