import time
import warnings

from incidence_producers.ForgiIncidenceProducer import ForgiIncidenceProducer
from incidence_producers.FornaIncidenceProducer import FornaIncidenceProducer
from rna_stats.FornaRnaStats import FornaRnaStats


def timefunc(function, *args):
    start_time = time.time()
    function(*args)
    print("--- %s seconds ---" % (time.time() - start_time))


if __name__ == "__main__":
    warnings.filterwarnings("ignore")  # non vengono scritti i warning di pandas
    forna_file = "examples/mus_musculus.json"

    for producer in [
        FornaIncidenceProducer(forna_file),
        ForgiIncidenceProducer(forna_file),
    ]:
        rna_stats = FornaRnaStats(producer)

        rna_stats.plot_hypergraph()
        # print(f"modularity={rna_stats.modularity()}")
        # print("-" * 70)
        #
        # rna_stats.plot_partitions_conductance()
        # # for partition, conductance in enumerate(rna_stats.partitions_conductance()):
        # #     print(f"{partition=}, {conductance=}")
        # # print("-" * 70)
        # # for edge, centrality in rna_stats.n_between_centrality(n=1).items():
        # #     print(f"{edge=}, {centrality=}")
        # rna_stats.plot_n_between_centrality(n=1)
        # rna_stats.plot_n_between_centrality(n=2)
    pass
