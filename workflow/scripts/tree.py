import tskit
import numpy as np
import click
import os
import tqdm
from concurrent.futures import ProcessPoolExecutor, as_completed

def process_single_simulation(sim_id, tree_name, input_folder, n_samples):
    """
    Worker function to process a single tree sequence file.
    """
     # Load the tree sequence
    ts_path = os.path.join(input_folder, f"{tree_name}_{sim_id}.trees")
    ts = tskit.load(ts_path)
        
    rng = np.random.default_rng()
    
    sample_nodes = [ind.nodes[0] for ind in ts.individuals()]

    tmrca_list = []

    tree = ts.first()

    # Randomly select pairs and calculate TMRCA
    for _ in tqdm.tqdm(range(n_samples), total=n_samples, desc=f"Sim {sim_id}"):
        chosen_nodes = rng.choice(sample_nodes, size=2, replace=False)

        try:
            tmrca = tree.tmrca(chosen_nodes[0], chosen_nodes[1])
            tmrca_list.append(tmrca)
        except ValueError:
            continue # Skip if no common ancestor (usually in multiple-tree TS)

    return sim_id, tmrca_list

@click.command()
@click.option("--tree-name", "-t", required=True, help="The name of the tree sequence file (without the .trees extension).")
@click.option("--input-folder", "-i", required=True, help="The folder containing the tree sequence file.")
@click.option("--output", "-o", default="tmp/results", help="The folder where the output will be saved.")
@click.argument("n_iter", type=int)
@click.argument("n_samples", type=int, default=100)
@click.argument("max_workers", type=int, default=max(1, (os.cpu_count() or 1) - 1), required=False)

def generate_tmrca_list(n_iter, tree_name, input_folder, n_samples, output, max_workers):
    """
    Parallelized generation of TMRCA values using ProcessPoolExecutor.
    """
    os.makedirs(output, exist_ok=True)
    output_path = os.path.join(output, f"{tree_name}.txt")
    
    # Initialize the file with a header
    with open(output_path, "w") as f:
        f.write("\t".join(["sim_id", "tmrca_list"]) + "\n")

    print(f"Starting parallel processing with {max_workers or 'all available'} cores...")

    # Use ProcessPoolExecutor for CPU-bound tasks
    with ProcessPoolExecutor(max_workers=max_workers) as executor:
        # Submit all tasks
        futures = {
            executor.submit(process_single_simulation, i, tree_name, input_folder, n_samples): i 
            for i in range(n_iter)
        }

        # Process results as they complete
        for future in as_completed(futures):
            sim_id, result = future.result()
            
            if isinstance(result, str) and result.startswith("Error"):
                print(f"Simulation {sim_id} failed: {result}")
                continue

            # Write to file sequentially to avoid race conditions
            tmrca_str = ", ".join(map(str, result))
            with open(output_path, "a") as f:
                f.write(f"{sim_id}\t{tmrca_str}\n")
            
            print(f"Finished simulation {sim_id}")


if __name__ == "__main__":
    generate_tmrca_list()