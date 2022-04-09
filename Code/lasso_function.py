from matplotlib.widgets import LassoSelector
from matplotlib.path import Path

# - Below function modification of stack overflow user Adrian Martin.

class SelectFromCollection(object):
    # Select indices from a matplotlib collection using `LassoSelector`.

    def __init__(self, ax, collection, alpha_other=0.3, facecolors=None):
        self.canvas = ax.figure.canvas
        self.collection = collection
        self.alpha_other = alpha_other

        self.xys = collection.get_offsets()
        self.Npts = len(self.xys)

        # Ensure that we have separate colors for each object
        self.fc = collection.get_facecolors()
        if len(self.fc) == 0:
            raise ValueError('Collection must have a facecolor')
        elif len(self.fc) == 1:
            self.fc = np.tile(self.fc, (self.Npts, 1))

        if facecolors is not None: self.fc = facecolors

        self.lasso = LassoSelector(ax, onselect=self.onselect)
        self.ind = []

    def onselect(self, verts):
        path = Path(verts)
        self.ind = np.nonzero(path.contains_points(self.xys))[0]
        self.fc[:, -1] = self.alpha_other
        self.fc[self.ind, -1] = 1
        self.collection.set_facecolors(self.fc)
        self.canvas.draw_idle()

    def disconnect(self):
        self.lasso.disconnect_events()
        self.fc[:, -1] = 1
        self.collection.set_facecolors(self.fc)
        self.canvas.draw_idle()


if __name__ == '__main__':
    import matplotlib.pyplot as plt

    x = scanned_dataframe['Mobs']; y = scanned_dataframe['Mreal']

    fig, ax = plt.subplots()
    pts = ax.scatter(x, y, s = 10)
    ax.set_xlabel('Instrumental Magnitude, $M_{inst.}$')
    ax.set_ylabel('Real Magnitude, $M_{real}$')
    # facecolors = plt.cm.jet(data.flatten())
    selector = SelectFromCollection(ax, pts)

    def accept(event):
        if event.key == "enter":
            global selected_point_idx
            selected_point_idx = selector.ind
            selector.disconnect()
            ax.set_title("")
            fig.canvas.draw()


    fig.canvas.mpl_connect("key_press_event", accept)
    ax.set_title("Press enter to accept selected points.")
    plt.show()