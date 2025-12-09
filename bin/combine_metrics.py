#!/usr/bin/env python3
"""
Combine individual method metric CSV files into a single comparison table.
Each input CSV has metrics as rows and a single method column.
Output CSV will have metrics as rows and methods as columns.
"""

import pandas as pd
import argparse
from pathlib import Path
import sys


def parse_args():
    parser = argparse.ArgumentParser(
        description='Combine individual metric CSV files into a single comparison table'
    )
    parser.add_argument(
        '--input-dir',
        type=str,
        required=True,
        help='Directory containing individual *_metrics.csv files'
    )
    parser.add_argument(
        '--output',
        type=str,
        default='combined_metrics.csv',
        help='Output filename for combined metrics table (default: combined_metrics.csv)'
    )
    parser.add_argument(
        '--pattern',
        type=str,
        default='*_metrics.csv',
        help='Glob pattern for input files (default: *_metrics.csv)'
    )
    return parser.parse_args()


def read_method_metrics(filepath):
    """
    Read a single method's metrics CSV file.
    Expected format: first column is metric name, second column is value.
    Returns the method name and a Series of metric values.
    """
    try:
        df = pd.read_csv(filepath, index_col=0)
        
        # Get method name from filename (remove _metrics.csv suffix)
        method_name = Path(filepath).stem.replace('_metrics', '')
        
        # Assuming second column contains the values
        if df.shape[1] >= 1:
            metric_values = df.iloc[:, 0]
            metric_values.name = method_name
            return method_name, metric_values
        else:
            print(f"Warning: {filepath} has unexpected format, skipping", file=sys.stderr)
            return None, None
            
    except Exception as e:
        print(f"Error reading {filepath}: {e}", file=sys.stderr)
        return None, None


def combine_metrics(input_dir, pattern='*_metrics.csv'):
    """
    Combine all metric CSV files in the input directory.
    Returns a DataFrame with metrics as rows and methods as columns.
    """
    input_path = Path(input_dir)
    
    if not input_path.exists():
        raise FileNotFoundError(f"Input directory does not exist: {input_dir}")
    
    # Find all matching CSV files
    metric_files = sorted(input_path.glob(pattern))
    
    if not metric_files:
        raise FileNotFoundError(f"No files matching pattern '{pattern}' found in {input_dir}")
    
    print(f"Found {len(metric_files)} metric files to combine", file=sys.stderr)
    
    # Read and combine all method metrics
    all_metrics = {}
    for filepath in metric_files:
        method_name, metrics = read_method_metrics(filepath)
        if method_name is not None:
            all_metrics[method_name] = metrics
            print(f"  - Loaded metrics for: {method_name}", file=sys.stderr)
    
    if not all_metrics:
        raise ValueError("No valid metric files could be read")
    
    # Combine into a single DataFrame
    combined_df = pd.DataFrame(all_metrics)
    
    # Sort columns (methods) alphabetically for consistency
    combined_df = combined_df.sort_index(axis=1)
    
    print(f"\nCombined table shape: {combined_df.shape[0]} metrics × {combined_df.shape[1]} methods", 
          file=sys.stderr)
    
    return combined_df


def main():
    args = parse_args()
    
    try:
        # Combine metrics
        combined_df = combine_metrics(args.input_dir, args.pattern)
        
        # Write output
        output_path = Path(args.input_dir) / args.output
        combined_df.to_csv(output_path)
        
        print(f"\n✓ Combined metrics saved to: {output_path}", file=sys.stderr)
        print(f"  Metrics: {', '.join(combined_df.index.tolist())}", file=sys.stderr)
        print(f"  Methods: {', '.join(combined_df.columns.tolist())}", file=sys.stderr)
        
        return 0
        
    except Exception as e:
        print(f"Error: {e}", file=sys.stderr)
        return 1


if __name__ == '__main__':
    sys.exit(main())
