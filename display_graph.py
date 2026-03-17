#!/usr/bin/env python3
"""
Simple script to plot two-column data with linear or log-log scale.
Usage: python plot.py <input_file> [options]
"""

import numpy as np
import matplotlib.pyplot as plt
import argparse
import sys

def main():
    parser = argparse.ArgumentParser(description='Plot two-column data')
    parser.add_argument('input_file', help='Text file with two columns of data')
    parser.add_argument('-o', '--output', default='plot.jpg', help='Output filename (default: plot.jpg)')
    parser.add_argument('-t', '--title', default='Data Plot', help='Plot title')
    parser.add_argument('-x', '--xlabel', default='X', help='X-axis label')
    parser.add_argument('-y', '--ylabel', default='Y', help='Y-axis label')
    parser.add_argument('-c', '--color', default='blue', help='Line color (default: blue)')
    parser.add_argument('--loglog', action='store_true', help='Use log-log scale')
    parser.add_argument('--grid', action='store_true', help='Add grid')
    
    args = parser.parse_args()
    
    # Read data
    try:
        data = np.loadtxt(args.input_file)
        if data.ndim != 2 or data.shape[1] < 2:
            print("Error: File must have at least 2 columns")
            sys.exit(1)
    except Exception as e:
        print(f"Error reading file: {e}")
        sys.exit(1)
    
    # Extract columns
    x = data[:, 0]
    y = data[:, 1]
    
    # Create plot
    plt.figure(figsize=(10, 6))
    plt.plot(x, y, color=args.color, linewidth=2)
    
    # Set scale
    if args.loglog:
        plt.xscale('log')
        plt.yscale('log')
    
    # Labels and title
    plt.title(args.title)
    plt.xlabel(args.xlabel)
    plt.ylabel(args.ylabel)
    
    # Grid
    if args.grid:
        plt.grid(True, alpha=0.3)
    
    # Save
    plt.savefig(args.output, dpi=300, bbox_inches='tight')
    print(f"Plot saved as: {args.output}")
    plt.close()

if __name__ == "__main__":
    main()