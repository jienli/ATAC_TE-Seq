#!/usr/bin/env python3
import sys
import numpy as np
from scipy.ndimage import gaussian_filter1d
import matplotlib.pyplot as plt

def load_bedgraph(path):
    """Load a 2-column bedGraph: col1=position, col2=coverage."""
    data = np.loadtxt(path)
    x = data[:,0]
    y = data[:,1]
    return x, y

def main(bedgraph_path):
    x, y = load_bedgraph(bedgraph_path)

    # 1) Smooth with a Gaussian kernel
    sigma = 100  # tune this
    y_smooth = gaussian_filter1d(y, sigma=sigma)

    # 2) Lift the floor
    baseline = np.percentile(y_smooth, 10)
    y_smooth += baseline * 0.2

    # 3) Squash the peaks
    scale = (y.max() * 0.8) / y_smooth.max()
    y_smooth *= scale

    # 4) Plot
    fig, (ax1, ax2) = plt.subplots(2,1, sharex=True, figsize=(8,3))
    ax1.fill_between(x, y, color="navy");    ax1.set_title("Original")
    ax2.fill_between(x, y_smooth, color="teal"); ax2.set_title("Smoothed")
    ax2.set_xlabel("Genomic position")
    plt.tight_layout()
    plt.show()

if __name__ == "__main__":
    if len(sys.argv)!=2:
        print("Usage: python3 smooth_track.py <file.bedGraph>")
        sys.exit(1)
    main(sys.argv[1])