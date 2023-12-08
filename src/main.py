import time
import warnings
from pathlib import Path

from Bio import SeqIO

from hypergraph_folding.rna_folder import RNAFolder
from hypergraph_folding.temperature_hypergraph import TemperatureFoldingHypergraph, MemoryOptimizedFoldingHypergraph
from incidence_producers.vienna_incidence_producer import ViennaIncidenceProducer
from rna_stats.forna_rna_stats import TemperatureFoldingStats, RnaStats


def timefunc(function, *args):
    start_time = time.perf_counter()
    result = function(*args)
    end_time = time.perf_counter()
    total_time = end_time - start_time
    print(f"--- {total_time:.2f} seconds ---")
    return total_time, result


if __name__ == "__main__":
    warnings.filterwarnings("ignore")  # non vengono scritti i warning di pandas
    fasta_file = Path("../examples/sequence (3).fasta")
    sequence = str(next(SeqIO.parse(fasta_file, "fasta")).seq)
    # sequence = "CCACGCACCUAAUCGAUAUAACUGCUUAUGG"
    print(f"{len(sequence)=}")

    thg = TemperatureFoldingHypergraph(ViennaIncidenceProducer(RNAFolder(sequence)), MemoryOptimizedFoldingHypergraph())
    thg.insert_temperature_range(80, 100)
    stats = RnaStats(thg.get_hypergraph(80))
    t_stats = TemperatureFoldingStats(thg)

    stats.s_between_centrality()
    stats.connection_differences(thg.get_hypergraph(100))
    stats.get_nucleotides_change_structure(thg.get_hypergraph(100))
    stats.partitions_conductance()
    stats.structure_differences(thg.get_hypergraph(100))

    t_stats.get_nucleotide_sensibility_to_changes(80, 100)
    t_stats.get_structure_differences(80, 100)
pass
