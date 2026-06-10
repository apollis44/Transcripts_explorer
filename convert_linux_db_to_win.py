import shelve
import dbm.dumb
import os

def convert_linux_to_universal(source_db_name):
    print(f"Opening native Linux database: {source_db_name}")
    
    # 1. Open the true Linux database file (e.g., genes or genes.db)
    # On Linux, shelve automatically handles the native binary driver
    with shelve.open(source_db_name, flag='r') as src_shelf:
        data = {key: src_shelf[key] for key in src_shelf.keys()}
        
    print(f"Extracted {len(data)} items. Saving to universal format...")

    # 2. Write it out using dbm.dumb, which splits it into portable .dat and .dir files
    # We use a temporary target name
    target_name = f"{source_db_name}"
    
    # Clean up old target files if they exist to prevent corruption
    for ext in ['.dat', '.dir', '.bak']:
        if os.path.exists(target_name + ext):
            os.remove(target_name + ext)
            
    dst_db = dbm.dumb.open(target_name, flag='c')
    with shelve.Shelf(dst_db) as dst_shelf:
        for key, value in data.items():
            dst_shelf[key] = value

    print(f"Success! Created portable files: {target_name}.dat and {target_name}.dir\n")

if __name__ == "__main__":
    # If your files are in the current directory
    current_dir = os.path.dirname(os.path.abspath(__file__))
    plot_dir = os.path.join(current_dir, "/files_for_plots")
    
    shelves_to_convert = [
        "deeploc2_output", 
        "TCGA_GTEx_plotting_data", 
        "transcripts_to_isoforms_mapping", 
        "membrane_topology_objects"
    ]
    
    for name in shelves_to_convert:
        # Strip '.db' extension if you typed it in the list, 
        # shelve.open expects the base filename
        base_name = os.path.join(plot_dir, name.replace(".db", ""))
        print(base_name)        
        # Check if the file (or file.db) exists
        if os.path.exists(base_name) or os.path.exists(f"{base_name}.db"):
            try:
                convert_linux_to_universal(base_name)
            except Exception as e:
                print(f"Error converting {name}: {e}\n")
        else:
            print(f"Skipping {name}: Base file not found.\n")