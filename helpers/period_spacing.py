import matplotlib.pyplot as plt
import numpy as np

from utils import get_delta_p, mass_gainer_col, single_col,\
    find_closest_model_number, asymptotic_period_spacing, fs

__all__ = ["plot_period_spacing"]


def plot_period_spacing(age=None, X_c=None, tracks=None, labels=["Mass-gainer", "Single"],
                        colours=[mass_gainer_col, single_col], legend_loc="upper left", label_with="an",
                        x_var="period", label_modes=False, xlims=None, ylims=None, divide_delta_n=False,
                        drop_duplicate_ng=True, fig=None, ax=None, show=True, ylim_auto_fac=2):
    if age is None and X_c is None:
        raise ValueError("At least one of `age` or `X_c` must not be None")
    
    if fig is None or ax is None:
        fig, ax = plt.subplots(figsize=(7, 3))
    plt.cla()

    delta_p_max = -1
    
    for track, tag, col in zip(tracks, labels, colours):
        mod = find_closest_model_number(track=track, age=age, X_c=X_c)
        periods, ng, delta_p = get_delta_p(track=track, mod=mod, drop_duplicate_ng=drop_duplicate_ng)

        if divide_delta_n:
            delta_n = np.diff(ng)
            delta_n[delta_n == 0] = 1
            delta_n = np.concatenate((delta_n, [-1]))
            delta_p /= -delta_n

        x_vals = periods if x_var == "period" else -ng
    
        if label_modes:
            ax.plot(x_vals, delta_p, label=tag, color=col)
        else:
            ax.plot(x_vals, delta_p, marker="o", label=tag, color=col)

        if label_modes:
            for i in range(len(delta_p)):
                ax.annotate(ng[i], xy=(x_vals[i], delta_p[i]), ha="center", va="center",
                            fontsize=0.25*fs, color="grey",
                            bbox=dict(boxstyle="circle", facecolor="white", ec=col))

        aps = asymptotic_period_spacing(track.profiles[mod - 1]).to(u.day).value
        delta_p_max = max(delta_p_max, np.max(np.abs(delta_p.value - aps)))
        ax.axhline(aps, color=col, alpha=0.5, linestyle="dotted", zorder=-1)#, label=f"Expected asymptotic spacing ({tag})")

    ax.set_xlabel(r"Period, $P \, [\rm days]$" if x_var == "period" else r"$k$ (g-modes)")
    ax.set_ylabel(r"$\Delta P \, [\rm days]$")
    
    if label_with == "title":
        ax.set_title(r"Period spacing for final mass ~3.5 $\rm M_{\odot}$ star"\
                     + (f' at {age:1.1f} Myr' if X_c is None else r" with $X_c =$" + f' {X_c:1.2f}'))
    else:
        ax.annotate(f'Age = {age:1.1f} Myr' if X_c is None else r"$X_c =$" + f' {X_c:1.2f}',
                    xy=(0.02, 0.95), xycoords="axes fraction", va="top", fontsize=0.5*fs)

    ax.annotate(r"$(l = 1, m = 0)$ g modes", xy=(0.02, 0.02), xycoords="axes fraction", va="bottom")

    if ylims is not None:
        if ylims == "auto":
            ax.set_ylim(aps - delta_p_max * ylim_auto_fac, aps + delta_p_max * ylim_auto_fac)
        else:
            ax.set_ylim(ylims)

    if xlims is not None:
        ax.set_xlim(xlims)

    ax.legend(loc=legend_loc, fontsize=0.5 * fs, ncols=2)
    
    if show:
        plt.show()

    return fig, ax, ng, periods