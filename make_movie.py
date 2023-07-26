import os

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation



# amp = 1.0
amp = 1e7
dt = 0.0001

def load_matrix(filename):
    with open(filename, 'r') as f:
        data = np.array([d.split(",")[:-1] for d in f.readlines()], dtype=float)

    return data


def create_animation(files, show=True, savefn=None):
    data = load_matrix(files[0])
    fig, ax = plt.subplots(1, 1)
    im = ax.imshow(data, cmap='seismic', vmin=-amp, vmax=amp)
    ax.set_xlabel("X [m]")
    ax.set_ylabel("Z [m]")
    fig.colorbar(im)

    def update_image(i):
        data = load_matrix(files[i+1])
        print(np.max(data))
        im.set_data(data)
        ax.set_title(f"2D Elastic Wave Propagation\nTime: {i*dt:10.5f} [s]")

    ani = FuncAnimation(fig, update_image, len(files)-1, interval=50)

    if show:
        plt.show()

    if savefn is not None:
        ani.save(savefn + ".mp4")

    return ani


def main():
    # folder = "pressure"
    folder = "vz"
    files = [os.path.join(folder, f) for f in sorted(os.listdir(folder))]

    create_animation(files)


if __name__ == "__main__":
    main()
