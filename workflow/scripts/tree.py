import tskit
import numpy as np

ts = tskit.load("tmp/results/normal_N5000_U0.006_s0.002_sigma0.001.trees")

print(ts)

reduced_ts = ts.simplify()

svg = reduced_ts.draw_svg()

with open("tmp/results/tree.svg", "w") as f:
    f.write(svg)

np.random.seed(42)

pairwise_samples = np.random.choice(ts.samples(), size=2, replace=False)

node_load = {node: 0.0 for node in pairwise_samples}

# Iterate over all sites that have mutations
for variant in ts.variants(samples=pairwise_samples):
    # variant.genotypes is an array of 0s and 1s 
    # corresponding to the order of sample_nodes
    
    # Get the 's' for the mutation(s) at this site
    # (Simplified: assuming one mutation per site for this example)
    mut_id = variant.site.mutations[0].id
    s = ts.mutation(mut_id).metadata["mutation_list"][0]["selection_coeff"]
    
    # Add 's' to the load of every individual who carries the mutation (genotype == 1)
    for i, has_mutation in enumerate(variant.genotypes):
        if has_mutation > 0:
            node_id = pairwise_samples[i]
            node_load[node_id] += s

# Print results
for node, load in list(node_load.items())[:10]:
    print(f"Node {node} total load: {load}")