import csv
import matplotlib.pyplot as plt
from pymol import cmd
import os

def load_residue_rates_from_csv(csv_path, rate_column='c_rate'):
    """
    Load residue identifiers and their rate values from a CSV file.
    Expects a header row with 'chain', 'residue', and rate columns.

    Parameters:
    - csv_path: Path to the CSV file
    - rate_column: Which column to use for values ('rate', 'category', or 'c_rate')

    Returns:
    - A dictionary mapping (chain, residue) tuples to numerical values.
    """
    residue_value_dict = {}
    skipped_count = 0
    
    with open(csv_path, newline='', encoding='utf-8-sig') as csvfile:
        reader = csv.DictReader(csvfile)
        for row in reader:
            try:
                chain = row['chain'].strip()
                residue = row['residue'].strip()
                
                # Try to get the rate value
                value_str = row[rate_column].strip()
                if not value_str:
                    skipped_count += 1
                    continue
                    
                value = float(value_str)
                residue_value_dict[(chain, residue)] = value
            except (KeyError, ValueError) as e:
                skipped_count += 1
                continue
    
    if skipped_count > 0:
        print(f"Skipped {skipped_count} residues with missing or invalid {rate_column} values")
    
    return residue_value_dict

def color_residues_by_value(residue_value_dict, prefix, cmap_name="RdYlBu_r", 
                            min_val=None, max_val=None, selection="all",
                            binary_mode=False, threshold=None, 
                            low_color="blue", high_color="red"):
    """
    Colors residues based on a dictionary of (chain, residue)-to-value mappings 
    using either a color gradient or binary coloring.
    
    Parameters:
    - residue_value_dict: Dictionary mapping (chain, residue) tuples to numerical values.
    - prefix: PyMOL object name (e.g., 'structure').
    - cmap_name: Name of a Matplotlib colormap (e.g., "RdYlBu_r", "viridis", "coolwarm").
    - min_val: Minimum value for color normalization (optional).
    - max_val: Maximum value for color normalization (optional).
    - selection: Additional selection criteria (e.g., "name CA" to color only CA atoms)
    - binary_mode: If True, use two colors based on threshold
    - threshold: Dividing value for binary coloring (values < threshold get low_color)
    - low_color: Color for values below threshold (color name or RGB tuple)
    - high_color: Color for values above/equal threshold (color name or RGB tuple)
    """
    if not residue_value_dict:
        print("No residues provided for coloring.")
        return

    colored_count = 0
    missing_count = 0
    low_count = 0
    high_count = 0
    
    if binary_mode:
        # Binary coloring mode
        if threshold is None:
            print("✗ Binary mode requires a threshold value")
            return
        
        # Parse color specifications
        def parse_color(color_spec):
            if isinstance(color_spec, str):
                # Try to parse as RGB tuple string like "(1,0,0)" or "1,0,0"
                if ',' in color_spec:
                    color_spec = color_spec.strip('()').strip()
                    try:
                        rgb = tuple(float(x.strip()) for x in color_spec.split(','))
                        if len(rgb) == 3:
                            return rgb
                    except ValueError:
                        pass
                # Otherwise treat as PyMOL color name
                return color_spec
            return color_spec
        
        low_color_parsed = parse_color(low_color)
        high_color_parsed = parse_color(high_color)
        
        # Set up the two colors
        if isinstance(low_color_parsed, tuple):
            cmd.set_color("binary_low", low_color_parsed)
            low_color_name = "binary_low"
        else:
            low_color_name = low_color_parsed
            
        if isinstance(high_color_parsed, tuple):
            cmd.set_color("binary_high", high_color_parsed)
            high_color_name = "binary_high"
        else:
            high_color_name = high_color_parsed
        
        for (chain, residue), value in residue_value_dict.items():
            # Choose color based on threshold
            if value < threshold:
                color_name = low_color_name
                low_count += 1
            else:
                color_name = high_color_name
                high_count += 1
            
            # Create selection for this specific residue
            res_selection = f"{prefix} and chain {chain} and resi {residue}"
            if selection != "all":
                res_selection += f" and {selection}"
            
            atom_count = cmd.count_atoms(res_selection)
            
            if atom_count > 0:
                cmd.color(color_name, res_selection)
                colored_count += 1
            else:
                missing_count += 1
                if missing_count <= 5:
                    print(f"⚠ Residue not found: chain {chain}, residue {residue}")
        
        if missing_count > 5:
            print(f"⚠ ... and {missing_count - 5} more residues not found")
        
        print(f"\n✓ Colored {colored_count} residues using binary coloring.")
        print(f"  Threshold: {threshold}")
        print(f"  {low_count} residues < {threshold} ({low_color_name})")
        print(f"  {high_count} residues ≥ {threshold} ({high_color_name})")
        if missing_count > 0:
            print(f"  {missing_count} residues from CSV not found in structure")
    
    else:
        # Gradient coloring mode (original behavior)
        # Determine normalization range
        if min_val is None:
            min_val = min(residue_value_dict.values())
        if max_val is None:
            max_val = max(residue_value_dict.values())

        value_range = max_val - min_val or 1  # Avoid division by zero
        cmap = plt.get_cmap(cmap_name)
        
        for i, ((chain, residue), value) in enumerate(residue_value_dict.items()):
            # Normalize value to [0, 1]
            normalized_value = (value - min_val) / value_range
            rgb = cmap(normalized_value)[:3]  # ignore alpha
            color_name = f"rate_color_{i}"

            cmd.set_color(color_name, rgb)
            
            # Create selection for this specific residue
            res_selection = f"{prefix} and chain {chain} and resi {residue}"
            if selection != "all":
                res_selection += f" and {selection}"
            
            atom_count = cmd.count_atoms(res_selection)
            
            if atom_count > 0:
                cmd.color(color_name, res_selection)
                colored_count += 1
            else:
                missing_count += 1
                if missing_count <= 5:  # Only print first few missing residues
                    print(f"⚠ Residue not found: chain {chain}, residue {residue}")

        if missing_count > 5:
            print(f"⚠ ... and {missing_count - 5} more residues not found")
        
        print(f"\n✓ Colored {colored_count} residues using the '{cmap_name}' colormap.")
        print(f"  Value range: {min_val:.4f} to {max_val:.4f}")
        if missing_count > 0:
            print(f"  {missing_count} residues from CSV not found in structure")
            
def create_rate_legend(cmap_name="RdYlBu_r", min_val=0, max_val=1, 
                       title="Substitution Rate", n_colors=10):
    """
    Create a color legend/scale bar showing the rate categories.
    
    Parameters:
    - cmap_name: Colormap name
    - min_val: Minimum rate value
    - max_val: Maximum rate value
    - title: Title for the legend
    - n_colors: Number of discrete colors to show
    """
    import numpy as np
    
    cmap = plt.get_cmap(cmap_name)
    value_range = max_val - min_val or 1
    
    # Create objects representing the color scale
    for i in range(n_colors):
        # Value from min to max
        value = min_val + (i / (n_colors - 1)) * value_range
        normalized = (value - min_val) / value_range
        rgb = cmap(normalized)[:3]
        
        color_name = f"legend_{i}"
        cmd.set_color(color_name, rgb)
        
        # Create a pseudoatom for visualization
        cmd.pseudoatom(f"legend_scale", pos=[i * 2, 0, 0], color=color_name, 
                      vdw=1.5, label=f"{value:.2f}")
    
    print(f"\n✓ Created legend scale object 'legend_scale'")
    print(f"  Use 'hide everything, legend_scale' to hide it")
    print(f"  Use 'show spheres, legend_scale' to show it")

def pymol_rate_colors(csv_path, prefix, rate_column='c_rate', 
                      cmap_name="RdYlBu_r", min_val=None, max_val=None,
                      selection="all", create_legend=False):
    """
    Main function for use in PyMOL to color residues by site rates.
    
    Parameters:
    - csv_path: Path to CSV file with 'chain', 'residue', and rate columns.
    - prefix: Molecule object name (e.g., 'structure').
    - rate_column: Which rate column to use ('rate', 'category', or 'c_rate'). Default: 'c_rate'
    - cmap_name: Name of a matplotlib colormap. Default: 'RdYlBu_r' (red=fast, blue=slow)
    - min_val: Optional minimum value for normalization.
    - max_val: Optional maximum value for normalization.
    - selection: Additional PyMOL selection (e.g., "name CA"). Default: "all"
    - create_legend: Whether to create a visual legend. Default: False
    
    Example usage in PyMOL:
        pymol_rate_colors buS10_mapping.csv, structure, c_rate
        pymol_rate_colors buS10_mapping.csv, structure, category, viridis
        pymol_rate_colors buS10_mapping.csv, structure, rate, coolwarm, 0, 2
    """
    if not os.path.isfile(csv_path):
        print(f"✗ File not found: {csv_path}")
        return

    if not prefix:
        print("✗ You must specify a molecule object name (e.g., 'structure').")
        return
    
    # Check if object exists
    if prefix not in cmd.get_names('objects'):
        print(f"✗ Object '{prefix}' not found in PyMOL")
        print(f"  Available objects: {cmd.get_names('objects')}")
        return

    # Convert min and max values from strings if provided
    try:
        min_val = float(min_val) if min_val else None
        max_val = float(max_val) if max_val else None
    except (ValueError, TypeError):
        print("✗ Invalid min or max value. Please provide numeric values.")
        return

    # Validate rate column
    valid_columns = ['rate', 'category', 'c_rate']
    if rate_column not in valid_columns:
        print(f"✗ Invalid rate column '{rate_column}'. Must be one of: {valid_columns}")
        return

    print(f"\n{'='*70}")
    print(f"Coloring residues by {rate_column}")
    print(f"{'='*70}")
    
    # Load and color residues
    residue_rates = load_residue_rates_from_csv(csv_path, rate_column)
    
    if not residue_rates:
        print("✗ No valid rate data found in CSV file")
        return
    
    print(f"Loaded {len(residue_rates)} residues with {rate_column} values")
    
    color_residues_by_value(residue_rates, prefix, cmap_name, min_val, max_val, selection)
    
    # Create legend if requested
    if create_legend:
        actual_min = min_val if min_val is not None else min(residue_rates.values())
        actual_max = max_val if max_val is not None else max(residue_rates.values())
        create_rate_legend(cmap_name, actual_min, actual_max, 
                          title=f"{rate_column.upper()} Scale")
    
    print(f"{'='*70}\n")

# Convenience functions for common use cases
def color_by_rate_category(csv_path, prefix):
    """
    Color residues by their rate category (discrete values 0-10).
    Uses a diverging colormap where blue=slow, red=fast.
    """
    pymol_rate_colors(csv_path, prefix, rate_column='category', 
                     cmap_name='RdYlBu_r', min_val=0, max_val=10)

def color_by_c_rate(csv_path, prefix):
    """
    Color residues by their categorical rate (c_rate).
    Uses a diverging colormap where blue=slow, red=fast.
    """
    pymol_rate_colors(csv_path, prefix, rate_column='c_rate', 
                     cmap_name='RdYlBu_r')

def color_by_mean_rate(csv_path, prefix):
    """
    Color residues by their posterior mean rate.
    Uses a diverging colormap where blue=slow, red=fast.
    """
    pymol_rate_colors(csv_path, prefix, rate_column='rate', 
                     cmap_name='RdYlBu_r')

# Register functions for PyMOL
cmd.extend("pymol_rate_colors", pymol_rate_colors)
cmd.extend("color_by_rate_category", color_by_rate_category)
cmd.extend("color_by_c_rate", color_by_c_rate)
cmd.extend("color_by_mean_rate", color_by_mean_rate)

# Print usage instructions when script is loaded
print("\n" + "="*70)
print("PyMOL Rate Coloring Script Loaded")
print("="*70)
print("\nAvailable commands:")
print("  pymol_rate_colors CSV_PATH, OBJECT_NAME, RATE_COLUMN, [COLORMAP], [MIN], [MAX]")
print("  color_by_rate_category CSV_PATH, OBJECT_NAME")
print("  color_by_c_rate CSV_PATH, OBJECT_NAME")
print("  color_by_mean_rate CSV_PATH, OBJECT_NAME")
print("\nExamples:")
print("  pymol_rate_colors buS10_mapping.csv, testing, c_rate")
print("  pymol_rate_colors buS10_mapping.csv, testing, category, viridis")
print("  color_by_c_rate buS10_mapping.csv, testing")
print("\nColormap options: RdYlBu_r (default), viridis, coolwarm, plasma, etc.")
print("="*70 + "\n")