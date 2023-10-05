import warnings

from ForgiIncidenceProducer import ForgiIncidenceProducer
from FornaRnaStats import FornaRnaStats

if __name__ == "__main__":
    warnings.filterwarnings("ignore")  # non vengono scritti i warning di pandas
    rna_stats = FornaRnaStats(ForgiIncidenceProducer("examples/mus_musculus.json"))

    rna_stats.plot_hypergraph()
    print(f"modularity={rna_stats.modularity()}")
    # print("-" * 70)

    # for partition, conductance in rna_stats.get_partitions_conductance():
    #     print(f"{partition=}, {conductance=}")
    rna_stats.plot_partitions_conductance()

    # print("-" * 70)
    # for edge, centrality in rna_stats.get_n_between_centrality(n=1).items():
    #     print(f"{edge=}, {centrality=}")
    rna_stats.plot_n_between_centrality(n=1)
    rna_stats.plot_n_between_centrality(n=2)
    pass
