#!/usr/bin/env python3
"""
Generate comparison plots for scIB benchmark metrics.
Creates barplots comparing methods across metrics and heatmaps showing differences.
"""

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import argparse
from pathlib import Path
import numpy as np
import sys


def parse_args():
    parser = argparse.ArgumentParser(
        description='Generate comparison plots for benchmark metrics'
    )
    parser.add_argument(
        '--input',
        type=str,
        required=True,
        help='Path to combined_metrics.csv file'
    )
    parser.add_argument(
        '--reference',
        type=str,
        default=None,
        help='Path to reference combined_metrics.csv for comparison (optional)'
    )
    parser.add_argument(
        '--output-dir',
        type=str,
        default='.',
        help='Directory for output plots (default: current directory)'
    )
    parser.add_argument(
        '--format',
        type=str,
        default='png',
        choices=['png', 'pdf', 'svg'],
        help='Output plot format (default: png)'
    )
    parser.add_argument(
        '--dpi',
        type=int,
        default=300,
        help='DPI for raster formats (default: 300)'
    )
    return parser.parse_args()


def plot_metric_barplots(df, output_dir, fmt='png', dpi=300):
    """
    Create individual barplots for each metric across methods.
    """
    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)
    
    # Set style
    sns.set_style("whitegrid")
    
    for metric in df.index:
        fig, ax = plt.subplots(figsize=(10, 6))
        
        metric_values = df.loc[metric].dropna()
        
        # Create barplot
        bars = ax.bar(range(len(metric_values)), metric_values.values, 
                      color=sns.color_palette("husl", len(metric_values)))
        
        ax.set_xticks(range(len(metric_values)))
        ax.set_xticklabels(metric_values.index, rotation=45, ha='right')
        ax.set_ylabel('Value', fontsize=12)
        ax.set_title(f'{metric}', fontsize=14, fontweight='bold')
        ax.grid(axis='y', alpha=0.3)
        
        # Add value labels on bars
        for bar in bars:
            height = bar.get_height()
            ax.text(bar.get_x() + bar.get_width()/2., height,
                   f'{height:.3f}',
                   ha='center', va='bottom', fontsize=9)
        
        plt.tight_layout()
        
        # Save plot
        filename = f"barplot_{metric.replace(' ', '_')}.{fmt}"
        plt.savefig(output_path / filename, dpi=dpi, bbox_inches='tight')
        plt.close()
        
        print(f"  Created: {filename}", file=sys.stderr)


def plot_combined_barplot(df, output_dir, fmt='png', dpi=300):
    """
    Create a combined barplot with all metrics, grouped by method.
    """
    output_path = Path(output_dir)
    
    # Prepare data for grouped barplot
    df_melted = df.T.reset_index().melt(id_vars='index', var_name='Metric', value_name='Value')
    df_melted.rename(columns={'index': 'Method'}, inplace=True)
    
    # Create plot
    fig, ax = plt.subplots(figsize=(14, 8))
    
    sns.barplot(data=df_melted, x='Method', y='Value', hue='Metric', ax=ax)
    
    ax.set_xlabel('Integration Method', fontsize=12, fontweight='bold')
    ax.set_ylabel('Metric Value', fontsize=12, fontweight='bold')
    ax.set_title('Comparison of Integration Methods Across Metrics', 
                fontsize=14, fontweight='bold')
    ax.legend(title='Metric', bbox_to_anchor=(1.05, 1), loc='upper left')
    plt.xticks(rotation=45, ha='right')
    plt.grid(axis='y', alpha=0.3)
    
    plt.tight_layout()
    
    filename = f"barplot_all_metrics.{fmt}"
    plt.savefig(output_path / filename, dpi=dpi, bbox_inches='tight')
    plt.close()
    
    print(f"  Created: {filename}", file=sys.stderr)


def plot_difference_heatmap(df1, df2, output_dir, fmt='png', dpi=300, diff_type='absolute'):
    """
    Create heatmap showing differences between two metric tables.
    
    Parameters:
    -----------
    df1 : DataFrame
        New/current metrics
    df2 : DataFrame
        Reference/original metrics
    diff_type : str
        Type of difference: 'absolute', 'fold_change', or 'log2fc'
    """
    output_path = Path(output_dir)
    
    # Align dataframes (common metrics and methods)
    common_metrics = df1.index.intersection(df2.index)
    common_methods = df1.columns.intersection(df2.columns)
    
    if len(common_metrics) == 0 or len(common_methods) == 0:
        print(f"Warning: No common metrics/methods found for comparison", file=sys.stderr)
        return
    
    df1_aligned = df1.loc[common_metrics, common_methods]
    df2_aligned = df2.loc[common_metrics, common_methods]
    
    # Calculate differences
    if diff_type == 'absolute':
        diff = df1_aligned - df2_aligned
        title = 'Absolute Difference (New - Reference)'
        cmap = 'RdBu_r'
        center = 0
    elif diff_type == 'fold_change':
        diff = df1_aligned / df2_aligned
        title = 'Fold Change (New / Reference)'
        cmap = 'RdBu_r'
        center = 1
    elif diff_type == 'log2fc':
        # Add small epsilon to avoid log(0)
        epsilon = 1e-10
        diff = np.log2((df1_aligned + epsilon) / (df2_aligned + epsilon))
        title = 'Log2 Fold Change'
        cmap = 'RdBu_r'
        center = 0
    else:
        raise ValueError(f"Unknown diff_type: {diff_type}")
    
    # Create heatmap
    fig, ax = plt.subplots(figsize=(12, 8))
    
    sns.heatmap(diff, annot=True, fmt='.3f', cmap=cmap, center=center,
                cbar_kws={'label': title}, linewidths=0.5, ax=ax)
    
    ax.set_title(f'{title}\n(Metrics × Methods)', fontsize=14, fontweight='bold')
    ax.set_xlabel('Integration Method', fontsize=12, fontweight='bold')
    ax.set_ylabel('Metric', fontsize=12, fontweight='bold')
    
    plt.tight_layout()
    
    filename = f"heatmap_diff_{diff_type}.{fmt}"
    plt.savefig(output_path / filename, dpi=dpi, bbox_inches='tight')
    plt.close()
    
    print(f"  Created: {filename}", file=sys.stderr)


def main():
    args = parse_args()
    
    try:
        # Read current metrics
        print(f"Reading metrics from: {args.input}", file=sys.stderr)
        df_current = pd.read_csv(args.input, index_col=0)
        print(f"  Loaded: {df_current.shape[0]} metrics × {df_current.shape[1]} methods\n", 
              file=sys.stderr)
        
        output_path = Path(args.output_dir)
        output_path.mkdir(parents=True, exist_ok=True)
        
        # Generate individual metric barplots
        print("Generating individual metric barplots...", file=sys.stderr)
        plot_metric_barplots(df_current, args.output_dir, args.format, args.dpi)
        
        # Generate combined barplot
        print("\nGenerating combined barplot...", file=sys.stderr)
        plot_combined_barplot(df_current, args.output_dir, args.format, args.dpi)
        
        # If reference provided, generate comparison plots
        if args.reference:
            print(f"\nReading reference metrics from: {args.reference}", file=sys.stderr)
            df_reference = pd.read_csv(args.reference, index_col=0)
            print(f"  Loaded: {df_reference.shape[0]} metrics × {df_reference.shape[1]} methods\n", 
                  file=sys.stderr)
            
            print("Generating comparison heatmaps...", file=sys.stderr)
            plot_difference_heatmap(df_current, df_reference, args.output_dir, 
                                  args.format, args.dpi, diff_type='absolute')
            plot_difference_heatmap(df_current, df_reference, args.output_dir, 
                                  args.format, args.dpi, diff_type='fold_change')
            plot_difference_heatmap(df_current, df_reference, args.output_dir, 
                                  args.format, args.dpi, diff_type='log2fc')
        
        print(f"\n✓ All plots saved to: {args.output_dir}", file=sys.stderr)
        return 0
        
    except Exception as e:
        print(f"Error: {e}", file=sys.stderr)
        import traceback
        traceback.print_exc()
        return 1


if __name__ == '__main__':
    sys.exit(main())
