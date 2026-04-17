import tskit
import numpy as np
import click
import os

@click.command()
@click.option("--tree-name", "-t", required=True, help="The name of the tree sequence file (without the .trees extension).")
@click.option("--input-folder", "-i", required=True, help="The folder containing the tree sequence file.")
@click.option("--output", "-o", default="tmp/results/", help="The folder where the output will be saved.")
@click.argument("n_iter", type=int)
@click.argument("n_samples", type=int, default=100)
def generate_tmrca_list(sim_id, tree_name, input_folder, output="tmp/results/", n_samples=100):
    """
    Generate a list of TMRCA values for pairs of samples from a tree sequence.

    Parameters
    ----------
    tree_name : The name of the tree sequence file (without the .trees extension).
    input_folder : The folder containing the tree sequence file.
    output : The folder where the output will be saved.
    n_samples : The number of TMRCA values to generate.
    """
    for sim_id in range(n_iter):
        ts = tskit.load(f"{input_folder}/{tree_name}_{sim_id}.trees")

        real_nodes = [ind.nodes[0] for ind in ts.individuals()]

        reduced_ts = ts.simplify(real_nodes)

        tree = reduced_ts.first()

        rng = np.random.default_rng()

        tmrca_list = []
        sample_nodes = []

        for ind in reduced_ts.individuals():
            sample_nodes.append(ind.nodes[0])

        for _ in range(n_samples):
            chosen_nodes = rng.choice(sample_nodes, size=2, replace=False)
            tmrca = tree.tmrca(chosen_nodes[0], chosen_nodes[1])
            tmrca_list.append(tmrca)

        file_exists = os.path.isfile(f"{output}/{tree_name}.txt")
        with open(f"{output}/{tree_name}.txt", "a") as f:
            header = ["sim_id", "tmrca_list"]
            if not file_exists:
                f.write("\t".join(header) + "\n")
            tmrca_str = ", ".join(map(str, tmrca_list))
            line = [str(sim_id), tmrca_str]
            f.write("\t".join(line) + "\n")

if __name__ == "__main__":
    generate_tmrca_list()