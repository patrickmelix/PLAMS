__all__ = ['plot_band_structure']

def plot_band_structure(x, y, labels, fermi_energy=None, zero=None, show=False):
    """
    Plots an electronic band structure from DFTB or BAND with matplotlib.

    To control the appearance of the plot you need to call ``plt.ylim(bottom, top)``, ``plt.title(title)``, etc. 
    manually outside this function.

    x: list of float
        Returned by AMSResults.get_band_structure()

    y: 2D numpy array of float
        Returned by AMSResults.get_band_structure()

    labels: list of str
        Returned by AMSResults.get_band_structure()

    fermi_energy: float
        Returned by AMSResults.get_band_structure(). Should have the same unit as ``y``.

    zero: None or float or one of 'fermi', 'vbmax', 'cbmin'
        Shift the curves so that y=0 is at the specified value. If None, no shift is performed. 'fermi', 'vbmax', and 'cbmin' require that the ``fermi_energy`` is not None. Note: 'vbmax' and 'cbmin' calculate the zero as the highest (lowest) eigenvalue smaller (greater) than or equal to ``fermi_energy``. This is NOT necessarily equal to the valence band maximum or conduction band minimum as calculated by the compute engine.

    show: bool
        If True, call plt.show() at the end
    """
    import matplotlib.pyplot as plt
    import numpy as np
    if zero is None:
        zero = 0
    elif zero == 'fermi':
        assert fermi_energy is not None
        zero = fermi_energy 
    elif zero == 'vbmax':
        assert fermi_energy is not None
        zero = y[y <= fermi_energy].max()
    elif zero == 'cbmin':
        assert fermi_energy is not None
        zero = y[y >= fermi_energy].min()

    plt.plot(x, y-zero, '-')
    tick_x = []
    tick_labels = []
    for xx, ll in zip(x, labels):
        if ll:
            if len(tick_x) == 0:
                tick_x.append(xx)
                tick_labels.append(ll)
                continue
            if np.isclose(xx, tick_x[-1]):
                if ll != tick_labels[-1]:
                    tick_labels[-1] += f',{ll}'
            else:
                tick_x.append(xx)
                tick_labels.append(ll)

    for xx in tick_x:
        plt.axvline(xx)

    if fermi_energy is not None:
        plt.axhline(fermi_energy-zero, linestyle='--')

    plt.xticks(ticks=tick_x, labels=tick_labels)

    if show:
        plt.show()
