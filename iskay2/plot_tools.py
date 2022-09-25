import matplotlib.pyplot as plt

def plot_map(themap, show=True, fnameout=None):
    '''Plots themap, axis show degrees.'''
    ax = plt.subplot(projection=themap.wcs)
    ax.coords[0].set_major_formatter('d.d')
    ax.coords[1].set_major_formatter('d.d')
    plt.imshow(themap, vmin=-150, vmax=150)
    
    if show:
        plt.show()
    elif fnamout is not None:
        plt.savefig("%s" % fnameout)
