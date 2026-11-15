#include <filesystem>
#include <sdsl/suffix_arrays.hpp>

#include "fms_index.h"
#include "fms_index2.h"
#include "doubled.h"
#include "query_pipeline.h"
#include "memory_usage.h"
#include "parser.h"
#include "version.h"
#include "compact.h"

#include <fstream>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <getopt.h>
#include <math.h>
#include <thread>

static int usage() {
  std::cerr << std::endl;
  std::cerr << "Program: FMSI 2 - a space-efficient index of k-mer sets based on masked superstrings." << std::endl;
  std::cerr << "Version: " << VERSION << std::endl;
  std::cerr << "Contact: Ondrej Sladky (ondra.sladky@gmail.com)" << std::endl << std::endl;
  std::cerr << "Usage:   fmsi <command> [options]" << std::endl << std::endl;
  std::cerr << "Commands:" << std::endl;
  std::cerr << "    index   - Index a masked superstring: one file <name>.fmsi (presets --extra-small, --small, --fast)." << std::endl;
  std::cerr << "    query   - Membership queries: for every k-mer of the queries, whether it is in the set." << std::endl;
  std::cerr << "    lookup  - Dictionary queries: the hash of every k-mer of the queries (-1 if it is not in the set)." << std::endl;
  std::cerr << "    access  - Dictionary queries: the k-mer of every hash, the inverse of lookup." << std::endl;
  std::cerr << "    export  - Print the masked superstring of an index." << std::endl;
  std::cerr << "    info    - Print what an index file holds." << std::endl;
  std::cerr << "    convert - Rewrite an index of an older format as a .fmsi file." << std::endl << std::endl;
  std::cerr << "Commands of FMSI 0.4 (its index of several files; `fmsi legacy` lists them):" << std::endl;
  std::cerr << "    legacy <command> - index, query, lookup, export, and the experimental set operations" << std::endl;
  std::cerr << "                       (union, inter, diff, symdiff, merge, compact)." << std::endl << std::endl;
  return 1;
}

static int usage_legacy() {
  std::cerr << std::endl;
  std::cerr << "Usage:   fmsi legacy <command> [options]" << std::endl << std::endl;
  std::cerr << "The commands of FMSI 0.4, on its index of several files (<prefix>.fmsi.*)." << std::endl << std::endl;
  std::cerr << "Command (stable):" << std::endl;
  std::cerr << "    index   - Creates a BWT based index of the given masked superstring." << std::endl;
  std::cerr << "    query   - Queries k-mers against an index." << std::endl;
  std::cerr << "    lookup  - Return unique hashes of present k-mers." << std::endl;
  std::cerr << "    export  - Print the underlying masked superstring to stdout." << std::endl << std::endl;
  std::cerr << "Command (experimental, using f-MS framework):" << std::endl;
  std::cerr << "    union   - Compute union of k-mers from several indices." << std::endl;
  std::cerr << "    inter   - Compute intersection of k-mers from several indices." << std::endl;
  std::cerr << "    diff    - Compute set difference of k-mers from several indices." << std::endl;
  std::cerr << "    symdiff - Compute symmetric difference of k-mers from several indices." << std::endl;
  std::cerr << "    merge   - Merges several indices." << std::endl;
  std::cerr << "    compact - Compacts the given index." << std::endl << std::endl;
  return 1;
}

static int usage_index() {
  std::cerr << std::endl;
  std::cerr << "Usage:   fmsi legacy index [options] <masked-superstring-input>" << std::endl << std::endl;
  std::cerr << "Options:" << std::endl;
  std::cerr << "    -k INT  - size of k-mers [recommended, default: number of mask trailing zeros - 1]"
            << std::endl;
  std::cerr << "    -x      - do not compute the kLCP array used for faster streaming queries."
            << std::endl;
  std::cerr << "    -Q      - construct the suffix array by QSufSort instead of libsais."
            << std::endl << std::endl;
  std::cerr << "Note: `fmsi legacy index` accepts only masked superstrings - these can be computed e.g. by KmerCamel from any FASTA file." << std::endl << std::endl;
  return 1;
}

static int usage_merge() {
  std::cerr << std::endl << 
      "Usage:   fmsi merge [options]"
      << std::endl << std::endl;
  std::cerr << std::endl << "Options:" << std::endl;
  std::cerr << "  `-p path_to_fasta`  - The path to the fasta file with "
                "input sets. Can be provided multiple times. It is "
                "expected that it appears at least twice."
            << std::endl;
    std::cerr << "  `-r path_of_result` - The path where the result should be "
                 "stored. Required."
              << std::endl;
    std::cerr << "  `-k value_of_k` - The size of queried k-mers (only used to check with the index one)."
              << std::endl;
    std::cerr << std::endl;
  return 1;
}

static int usage_op(std::string op) {
  std::cerr << std::endl << 
      "Usage:   fmsi " << op << " [options]"
      << std::endl << std::endl;
  std::cerr << std::endl << "Options:" << std::endl;
  std::cerr << "  `-p path_to_fasta`  - The path to the fasta file with "
                "input sets. Can be provided multiple times. It is "
                "expected that it appears at least twice."
            << std::endl;
    std::cerr << "  `-r path_of_result` - The path where the result should be "
                 "stored. Required."
              << std::endl;
    std::cerr << "  `-k value_of_k` - The size of queried k-mers (only used to check with the index one)."
              << std::endl;
    std::cerr << std::endl;
    return 1;
}

void usage_functions() {
  std::cerr << "  -f FUNCTION - Demasking function to determine k-mer presence; recognized functions:"
            << std::endl;
  std::cerr << "    or      - represented when at least 1 ON occurrence [default]"
            << std::endl;
  std::cerr << "    all     - all occurrence are either ON or "
               "OFF (equivalent to -O flag for queries)"
            << std::endl;
  std::cerr << "    and     - represented when no OFF occurrence"
            << std::endl;
  std::cerr << "    xor     - represented when an odd number "
               "of ON occurrences"
            << std::endl;
  std::cerr << "    INT-INT - represented when in the bounds"
            << std::endl;
}

static int usage_query() {
  std::cerr << std::endl;
  std::cerr << "Usage:   fmsi query [options] <index-prefix>" << std::endl << std::endl;
  std::cerr << "Options (stable):" << std::endl;
  std::cerr << "  -q FILE - Path to FASTA/FASTQ with queries [default: stdin]"
            << std::endl;
  std::cerr << "  -k INT  - Size of k-mers [default: infer automatically from index]"
            << std::endl;
  std::cerr << "  -S      - Use kLCP array for streamed queries (increses memory consumption)" << std::endl;
  std::cerr << "  -O      - FMSI uses properties of max-one masked superstrings to speed up queries" << std::endl;
  std::cerr << "            Use only if a masked superstring with maximum number of ones is indexed." << std::endl;
  std::cerr << "Parameters (experimental, using f-MS framework):" << std::endl;
  usage_functions();
  std::cerr << std::endl;
  return 1;
}

static int usage_lookup() {
  std::cerr << std::endl;
  std::cerr << "Usage:   fmsi lookup [options] <index-prefix>" << std::endl << std::endl;
  std::cerr << "Options (stable):" << std::endl;
  std::cerr << "  -q FILE - Path to FASTA/FASTQ with queries [default: stdin]"
            << std::endl;
  std::cerr << "  -k INT  - Size of k-mers [default: infer automatically from index]"
            << std::endl;
  std::cerr << "  -S      - Use kLCP array for streamed queries (increses memory consumption)" << std::endl;
  std::cerr << std::endl;
  return 1;
}

static int usage_query(bool lookup) {
  if (lookup) return usage_lookup();
  return usage_query();
}

static int usage_normalize() {
  std::cerr << std::endl << 
      "Usage:   fmsi compact [options] <index-prefix>"
      << std::endl << std::endl;
  std::cerr << std::endl << "Options (all experimental):" << std::endl;
  std::cerr << std::endl << "The recognized arguments are:" << std::endl;
  std::cerr << "    `-k value_of_k` - The size of queried k-mers."
            << std::endl;
  std::cerr << "    `-s` - Only print the compacted masked superstring and do not "
               "compact the index"
            << std::endl;
  usage_functions();
  std::cerr << std::endl;
  return 1;
}

static int usage_export() {
  std::cerr << std::endl << 
      "Usage:   fmsi export <index-prefix>"
      << std::endl << std::endl;
    return 1;
}

static int usage_clean() {
  std::cerr << std::endl << 
  "Usage:   fmsi clean <index-prefix>"
  << std::endl << std::endl;
  return 1;
}

static int version() {
  std::cout << VERSION << std::endl;
  return 0;
}

int ms_index(int argc, char *argv[]) {
  int c;
  bool usage = false;
  std::string fn;

  if (argc > 1 && std::string(argv[argc - 1]) != "-h") {
    fn = argv[argc - 1];
    argc--;
  }

  int k = 0;
  bool no_streaming = false;
  bool use_qsufsort = false;
  while ((c = getopt(argc, argv, "hk:xQ")) >= 0) {
    switch (c) {
    case 'h':
      usage = true;
      break;
    case 'k':
      k = atoi(optarg);
      break;
    case 'x':
      no_streaming = true;
      break;
    case 'Q':
      use_qsufsort = true;
      break;
    default:
      return usage_index();
    }
  }
  if (usage) {
    usage_index();
    return 0;
  } else if (fn.empty()) {
    std::cerr << "ERROR: Path to the masked superstring is a required argument." << std::endl;
    return usage_index();
  }

  std::cerr << "Starting " << fn << std::endl;
  auto ms = read_masked_superstring(fn);
  if (ms.size() == 0) {
    std::cerr << "ERROR: The file '" << fn
              << "' is in incorrect format. It is supposed to be a fasta file "
                 "with a single entry, the masked superstring"
              << std::endl;
    return usage_index();
  }
  std::cerr << "Read masked superstring of length " << ms.size() << std::endl;
  int inferred_k = infer_k(ms);
  if (k == 0) {
    k = inferred_k;
    std::cerr << "Inferred k from the masked case convention: " << k << std::endl;
  }
  if (k != inferred_k) {
      std::cerr << "WARNING: The provided k (" << k << ") does not match the k inferred from the mask convention (" << inferred_k << "). The provided k is used but we recommend double checking that it is correct." << std::endl;
  }
  if (k > 64 && !no_streaming) {
      std::cerr << "WARNING: Construction of kLCP array for streaming support is only available for k <= 64. The index will be constructed without streaming support, which results in slower positive streaming queries." << std::endl;
      no_streaming = true;
  }
  fms_index index;
  if (k <= 32) index = construct<uint64_t>(ms, k, !no_streaming, use_qsufsort);
  else index = construct<__uint128_t>(ms, k, !no_streaming, use_qsufsort);
  std::cerr << "Constructed index" << std::endl;
  dump_index(index, fn);
  std::cerr << "Written index" << std::endl;
  return 0;
}

int ms_query(int argc, char *argv[], bool output_orders) {
  bool usage = false;
  int c;
  int k = 0;
  std::string fn;

  if (argc > 1 && std::string(argv[argc - 1]) != "-h") {
    fn = argv[argc - 1];
    argc--;
  }

  std::string query_fn = "-";
  std::string f_name = "or";
  std::function<bool(size_t, size_t)> f = mask_function("or");
  bool has_klcp = false;
  while ((c = getopt(argc, argv, "f:hq:k:OS")) >= 0) {
    switch (c) {
    case 'f':
      try {
        f_name = optarg;
        f = mask_function(optarg);
      } catch (std::invalid_argument &) {
        std::cerr << "ERROR: Function '" << optarg << "' not recognized." << std::endl;
        return usage_query();
      }
      break;
    case 'h':
      usage = true;
      break;
    case 'q':
      query_fn = optarg;
      break;
    case 'k':
      k = atoi(optarg);
      break;
    case 'O':
        if (f_name != "or") {
            std::cerr << "WARNING: Parameter -O is ignored when parameter -f is specified." << std::endl;
        } else {
            f_name = "all";
            f = mask_function("all");
        }
        break;
    case 'S':
      has_klcp = true;
      break;
    default:
      return usage_query(output_orders);
    }
  }
  if (usage) {
    usage_query(output_orders);
    return 0;
  } else if (fn.empty()) {
    std::cerr << "ERROR: Path to the fasta file is a required argument." << std::endl;
    return usage_query(output_orders);
  }
  if (output_orders && f_name == "all") {
    std::cerr << "WARNING: The current version of FMSI has speed benefits only if output as (minimum) perfect hash function is not used. Additionally, if you desire minimum perfect hash function, please minimize the number of ones in the mask." << std::endl;
  } else if (f_name != "or" && output_orders) {
    std::cerr << "ERROR: FMSI as Minimum Perfect Hash Function is not allowed with f-masked superstrings in the current version." << std::endl;
    return usage_query(output_orders);
  }

  fms_index index = load_index(fn, has_klcp);

  if (index.sa_transformed_mask.size() == 0) {
    std::cerr << "ERROR: index not correctly loaded. Ensure that you correctly call `fmsi index` before." << std::endl;
    return usage_query(output_orders);
  }

  if (has_klcp != (index.klcp.size() > 0)) {
    std::cerr << "ERROR: kLCP array was not constructed for the given index. Either construct it again without the `-s` flag or use `query -s` which slows down streaming queries." << std::endl;
    return usage_query(output_orders);
  }
  int index_k = index.k;
  if (k != 0 && k != index_k) {
    std::cerr << "ERROR: Mismatch. Provided k (" << k << ") does not match the k of the index (" << index_k << ")." << std::endl;
    return usage_query(output_orders);
  }
  if (k == 0) {
    k = index_k;
  }

  gzFile fp = OpenFile(query_fn);
  kseq_t *seq = kseq_init(fp);
  int64_t sequence_length = 0;

  std::cin.tie(&std::cout);

  while ((sequence_length = kseq_read(seq)) >= 0) {
    // Small overhead for the chunking (while gaining superior time from prediction).
    int64_t max_sequence_chunk_length = 400;
    max_sequence_chunk_length = k + std::max((int64_t)10, std::min(max_sequence_chunk_length, 2*(int64_t)std::sqrt(sequence_length)));

    std::cout << seq->name.s << "\t";

    auto sequence = seq->seq.s;
    bool output_comma = false;
    while (sequence_length > 0) {
        int64_t current_length = next_invalid_character_or_end(sequence, sequence_length);
        
        while (current_length >= k) {
            if (output_orders && output_comma) std::cout << ",";
            output_comma = true;
            int64_t chunk_length = std::min(current_length, max_sequence_chunk_length);
            if (f_name == "or") {
                query_kmers<query_mode::orr>(index, sequence, chunk_length, k, has_klcp, std::cout, output_orders);
            } else if (f_name == "all") {
                query_kmers<query_mode::all>(index, sequence, chunk_length, k, has_klcp, std::cout, output_orders);
            } else {
                query_kmers<query_mode::general>(index, sequence, chunk_length, k, has_klcp, std::cout, output_orders, f);
            }
            sequence += chunk_length - k + 1;
            current_length -= chunk_length - k + 1;
            sequence_length -= chunk_length - k + 1;
        }
        // Skip also the next character.
        sequence_length -= current_length + 1;
        sequence += current_length + 1;
        // Print 0 on invalid k-mers.
        if (sequence_length >= 0) {
          for (int64_t i = 0; i < std::min((int64_t) k, current_length + 1); ++i) {
            if (output_orders) {
              if (output_comma) std::cout << ",";
              output_comma = true;
              std::cout << "-1";
            }
            else {
              std::cout << "0";
            }
          }
        }
    }
    std::cout << "\n";
  }
  return 0;
}

static int usage_index2() {
  std::cerr << std::endl;
  std::cerr << "Usage:   fmsi index [options] <masked-superstring>" << std::endl << std::endl;
  std::cerr << "Indexes a masked superstring (FASTA, lower case = masked out; e.g. from KmerCamel) for membership" << std::endl;
  std::cerr << "queries (fmsi query) and dictionary queries (fmsi lookup, fmsi access). The index is one file." << std::endl << std::endl;
  std::cerr << "Presets:" << std::endl;
  std::cerr << "    --extra-small - size first: the masked superstring alone, loaded in the most compact layout" << std::endl;
  std::cerr << "                    (about the memory of FMSI 0.4; queries about 1.5 times slower than with --small)" << std::endl;
  std::cerr << "    --small       - [default] the masked superstring alone, loaded in the fastest layout (the same file" << std::endl;
  std::cerr << "                    as --extra-small, about a sixth more memory)" << std::endl;
  std::cerr << "    --fast        - speed first: the masked superstring and its reverse complement (twice the space), so" << std::endl;
  std::cerr << "                    that a k-mer is searched as it is only, never as its reverse complement" << std::endl;
  std::cerr << "    The layout in memory is chosen when the index is loaded: the presets --extra-small and --small write" << std::endl;
  std::cerr << "    the same index and differ in the layout that the file recommends (fmsi query --extra-small, --small)." << std::endl << std::endl;
  std::cerr << "Options:" << std::endl;
  std::cerr << "    -o FILE - the index file [default: <masked-superstring>.fmsi]; .fmsi is added if missing" << std::endl;
  std::cerr << "    -k INT  - size of k-mers [default: number of trailing masked-out characters + 1]" << std::endl;
  std::cerr << "    -l      - also write the kLCP array (<index file>.klcp), which speeds up streamed queries (query -S)" << std::endl;
  std::cerr << "              of sequences that do not follow the superstring, e.g. of RNA-seq reads" << std::endl;
  std::cerr << "    -z, --compact-file - store the mask along the superstring, where it is a few runs: a file smaller by" << std::endl;
  std::cerr << "              0.2 to 0.7 bits per character, the same index once loaded; loading takes about 10 ns per" << std::endl;
  std::cerr << "              character on one thread instead of 2 (fmsi query -t loads on several)" << std::endl;
  std::cerr << "    -t INT  - number of threads [default: 1]; the suffix array is constructed on "
            << (SUFFIX_ARRAY_THREADED ? "all of them (compiled with OpenMP)" : "one thread (compiled without OpenMP)")
            << std::endl << std::endl;
  std::cerr << "Advanced:" << std::endl;
  std::cerr << "    -c INT  - layout in memory recommended by the file, which fmsi query -c overrides: bytes per block" << std::endl;
  std::cerr << "              of the BWT, 64 (the fastest), 128 or 256 (the smallest) [default: 256 for --extra-small, else 64]" << std::endl;
  std::cerr << "    -m STR  - likewise the encoding of the mask: plain (the fastest), packed (smaller) or rrr (slow)" << std::endl;
  std::cerr << "              [default: packed for --extra-small, else plain]" << std::endl;
  std::cerr << "    -M STR  - the type of the mask: min-one (every k-mer set once: lookup gives a minimal perfect hash)," << std::endl;
  std::cerr << "              max-one (every occurrence of a k-mer of the set is set), or auto: as the header of the masked" << std::endl;
  std::cerr << "              superstring states it (mask=min-one, as KmerCamel writes), else as computed [default: auto];" << std::endl;
  std::cerr << "              a stated type that the mask contradicts is replaced by the computed one, with a warning" << std::endl;
  std::cerr << "    --lean  - with --fast: merge the suffix arrays of the superstring and of its reverse complement instead" << std::endl;
  std::cerr << "              of constructing that of both (less memory; superstrings below 2^31 characters; temporary" << std::endl;
  std::cerr << "              files <index file>.tmp.sa_*); the index is the same" << std::endl;
  std::cerr << "    -Q      - construct the suffix array by QSufSort instead of libsais (an index of one strand only)" << std::endl << std::endl;
  return 1;
}

static int usage_query2(bool lookup) {
  std::cerr << std::endl;
  std::cerr << "Usage:   fmsi " << (lookup ? "lookup" : "query") << " [options] <index.fmsi>" << std::endl << std::endl;
  if (lookup) {
    std::cerr << "Dictionary queries: prints, for every sequence of the queries, the hashes of its k-mers (-1 for a k-mer" << std::endl;
    std::cerr << "that is not in the set). With a min-one mask the hashes are a minimal perfect hash, 0 to (number of" << std::endl;
    std::cerr << "k-mers - 1); fmsi access inverts them." << std::endl << std::endl;
  } else {
    std::cerr << "Membership queries: prints, for every sequence of the queries, a string of 1 (in the set) and 0 for its" << std::endl;
    std::cerr << "k-mers. Works with any mask." << std::endl << std::endl;
  }
  std::cerr << "Options:" << std::endl;
  std::cerr << "  -q FILE - FASTA/FASTQ with the queries, may be gzipped [default: stdin]" << std::endl;
  std::cerr << "  -S      - streamed queries: search the k-mers of a sequence by extending the match of the previous one" << std::endl;
  std::cerr << "            (for reads and longer sequences); uses the kLCP array if the index has one" << std::endl;
  std::cerr << "  -L      - do not use the kLCP array" << std::endl;
  std::cerr << "  -t INT  - number of threads searching; with more than one, the queries are also read and the output printed" << std::endl;
  std::cerr << "            by a thread each [default: 1]" << std::endl << std::endl;
  std::cerr << "Layout in memory (by default the one that the index file recommends, i.e. that of its preset):" << std::endl;
  std::cerr << "  --extra-small - the most compact layout (-c 256 -m packed)" << std::endl;
  std::cerr << "  --small       - the fastest layout (-c 64 -m plain)" << std::endl << std::endl;
  std::cerr << "Advanced:" << std::endl;
  std::cerr << "  -c INT  - bytes per block of the BWT in memory: 64 (the fastest), 128 or 256 (the smallest)" << std::endl;
  std::cerr << "            [default: as the index file recommends]" << std::endl;
  std::cerr << "  -m STR  - encoding of the mask in memory: plain (the fastest), packed (smaller) or rrr (slow)" << std::endl;
  std::cerr << "            [default: as the index file recommends]" << std::endl;
  std::cerr << "  -b INT  - number of k-mers (or sequences if -S) searched at the same time; 1 turns batching off [default: 32]" << std::endl;
  std::cerr << "  -p INT  - length of the patterns with tabulated SA intervals (8 * 4^INT bytes); 0 turns the table off" << std::endl;
  std::cerr << "            [default: by the size and the layout of the index, at most 12]" << std::endl;
  std::cerr << "  -B INT  - size of a block of sequences handed to a thread, in bytes of the query file [default: 131072]" << std::endl;
  std::cerr << "  -O      - accepted and ignored (the option of FMSI 0.4 for max-one masks; any mask is tested in constant time)" << std::endl;
  std::cerr << "  -k INT  - size of k-mers, checked against the index" << std::endl;
  std::cerr << std::endl;
  return 1;
}

/// The path given to -o or the masked superstring, with the extension of index files.
static std::string index_file_name(const std::string &name) {
  const std::string extension = fmsi2::INDEX_EXTENSION;
  if (name.size() >= extension.size() && name.compare(name.size() - extension.size(), extension.size(), extension) == 0) return name;
  return name + extension;
}

/// The options -c and -m of a layout; false if the value is not valid.
static bool parse_layout_option(int option, const char *value, fmsi2::layout &layout) {
  if (option == 'c') {
    layout.bwt_block_bytes = atoi(value);
    if (!fmsi2::quad_sequence::valid_block_bytes(layout.bwt_block_bytes)) {
      std::cerr << "ERROR: The blocks of the BWT can have 64, 128 or 256 bytes." << std::endl;
      return false;
    }
  } else {
    fmsi2::mask_kind kind;
    if (!fmsi2::parse_mask_name(value, kind)) {
      std::cerr << "ERROR: Mask encoding '" << value << "' not recognized." << std::endl;
      return false;
    }
    layout.mask = (int)kind;
  }
  return true;
}

int ms_index2(int argc, char *argv[]) {
  int c;
  bool usage = false;
  std::string fn;

  if (argc > 1 && std::string(argv[argc - 1]) != "-h") {
    fn = argv[argc - 1];
    argc--;
  }

  int k = 0;
  bool use_qsufsort = false;
  std::string prefix;
  // The layout in memory that the file recommends: that of the preset unless -c or -m is given.
  fmsi2::layout layout;
  std::string mask_type = "auto";
  bool use_klcp = false;
  int threads = 1;
  // Presets: --extra-small and --small (one strand, the compact or the fast layout) or --fast (both strands).
  bool fast = false, lean = false, extra_small = false;
  // The mask stored along the superstring: a smaller file, loaded more slowly.
  bool compact_file = false;
  enum { OPTION_SMALL = 256, OPTION_FAST, OPTION_LEAN, OPTION_EXTRA_SMALL };
  static const struct option long_options[] = {
      {"compact-file", no_argument, nullptr, 'z'},
      {"extra-small", no_argument, nullptr, OPTION_EXTRA_SMALL},
      {"small", no_argument, nullptr, OPTION_SMALL},
      {"fast", no_argument, nullptr, OPTION_FAST},
      {"lean", no_argument, nullptr, OPTION_LEAN},
      {nullptr, 0, nullptr, 0},
  };
  while ((c = getopt_long(argc, argv, "hk:m:c:lo:Qt:M:z", long_options, nullptr)) >= 0) {
    switch (c) {
    case OPTION_EXTRA_SMALL:
      fast = false;
      extra_small = true;
      break;
    case OPTION_SMALL:
      fast = extra_small = false;
      break;
    case OPTION_FAST:
      fast = true;
      extra_small = false;
      break;
    case OPTION_LEAN:
      lean = true;
      break;
    case 'h':
      usage = true;
      break;
    case 'k':
      k = atoi(optarg);
      break;
    case 'c':
    case 'm':
      if (!parse_layout_option(c, optarg, layout)) return usage_index2();
      break;
    case 'M':
      mask_type = optarg;
      if (mask_type != "auto" && mask_type != "min-one" && mask_type != "max-one") {
        std::cerr << "ERROR: The type of the mask is min-one, max-one or auto." << std::endl;
        return usage_index2();
      }
      break;
    case 'l':
      use_klcp = true;
      break;
    case 'z':
      compact_file = true;
      break;
    case 'o':
      prefix = optarg;
      break;
    case 'Q':
      use_qsufsort = true;
      break;
    case 't':
      threads = atoi(optarg);
      if (threads < 1) {
        std::cerr << "ERROR: The number of threads must be positive." << std::endl;
        return usage_index2();
      }
      break;
    default:
      return usage_index2();
    }
  }
  if (usage) {
    usage_index2();
    return 0;
  } else if (fn.empty()) {
    std::cerr << "ERROR: Path to the masked superstring is a required argument." << std::endl;
    return usage_index2();
  }
  const std::string index_path = index_file_name(prefix.empty() ? fn : prefix);
  const std::string preset = fast ? "fast" : extra_small ? "extra-small" : "small";
  const fmsi2::layout preset_layout = fmsi2::preset_layout(preset);
  if (layout.bwt_block_bytes == 0) layout.bwt_block_bytes = preset_layout.bwt_block_bytes;
  if (layout.mask < 0) layout.mask = preset_layout.mask;
  if (fast && use_qsufsort) {
    std::cerr << "ERROR: QSufSort (-Q) constructs the index of one strand only (--small, --extra-small)." << std::endl;
    return usage_index2();
  }
  if (lean && !fast) {
    std::cerr << "ERROR: The lean construction (--lean) is the one of the --fast index." << std::endl;
    return usage_index2();
  }

  build_phases = phase_timer();
  std::string header;
  auto ms = read_masked_superstring(fn, &header);
  if (ms.size() == 0) {
    std::cerr << "ERROR: The file '" << fn << "' does not contain a masked superstring." << std::endl;
    return usage_index2();
  }
  build_phases.lap("read");
  std::cerr << "Read masked superstring of length " << ms.size() << std::endl;
  if (k == 0) {
    k = infer_k(ms);
    std::cerr << "Inferred k from the masked case convention: " << k << std::endl;
  }

  size_t length = ms.size();
  // The mask along the superstring, and the rows of the sampled positions that the file then needs.
  fmsi2::text_mask along_text;
  if (compact_file) {
    along_text.runs = fmsi2::mask_runs(ms);
    along_text.characters = length;
    mbwt_sample_step = fmsi2::SAMPLE_STEP;
  }
  // The rows that start with the same k-mer tell the type of the mask, or check the one that is stated.
  mbwt_with_same_kmer = true;
  if (use_klcp && !fast && k > 65) {
    std::cerr << "ERROR: The kLCP array (-l) of the --small index is available for k <= 65." << std::endl;
    return 1;
  }
  mbwt m;
  if (fast && lean) {
    m = doubled::compute_mbwt_lean(std::move(ms), k, use_klcp, threads, index_path + ".tmp");
  } else if (fast) {
    m = doubled::compute_mbwt(std::move(ms), k, use_klcp, threads);
  } else {
    m = k <= 32 ? compute_mbwt<uint64_t>(ms, k, use_klcp, use_qsufsort, threads)
                : compute_mbwt<__uint128_t>(ms, k, use_klcp, use_qsufsort, threads);
  }
  std::cerr << "Constructed the masked Burrows-Wheeler transform" << (fast ? " of both strands" : "") << std::endl;
  std::string().swap(ms);
  // The index is written from the MBWT; its layout in memory is built when it is loaded.
  const fmsi2::mask_statistics statistics = fmsi2::decide_mask_type(fmsi2::classify_mask(m), mask_type, header_field(header, "mask"));
  if (!statistics.warning.empty()) std::cerr << "WARNING: " << statistics.warning << "." << std::endl;
  const fmsi2::description description =
      fmsi2::save_mbwt(m, k, preset, layout, statistics, index_path, threads, header, compact_file ? &along_text : nullptr);
  const double index_bytes = (double)std::filesystem::file_size(index_path);
  std::cerr << "Written " << index_path << ": " << (uint64_t)index_bytes << " bytes (" << 8 * index_bytes / length
            << " bits per character of the superstring), mask " << statistics.type << " (" << statistics.source << ")" << std::endl;
  const std::string klcp_path = index_path + fmsi2::KLCP_EXTENSION;
  if (!use_klcp && std::filesystem::exists(klcp_path)) {
    // The kLCP array of the index that this file replaces.
    std::filesystem::remove(klcp_path);
    std::cerr << "Removed " << klcp_path << ", the kLCP array of the previous index of this name" << std::endl;
  }
  if (use_klcp) {
    fmsi2::save_klcp_bits(m.klcp, k, description.uid(), klcp_path, threads);
    const double klcp_bytes = (double)std::filesystem::file_size(klcp_path);
    std::cerr << "Written " << klcp_path << ": " << (uint64_t)klcp_bytes << " bytes (" << 8 * klcp_bytes / length
              << " bits per character of the superstring)" << std::endl;
  }
  build_phases.lap("write");

  double total = 0;
  std::cerr << "[fmsi2-build]\tpreset=" << preset << "\tlean=" << lean << "\tstrands=" << description.strands
            << "\tthreads=" << threads << "\topenmp=" << SUFFIX_ARRAY_THREADED;
  for (const auto &phase : build_phases.phases) {
    std::cerr << "\t" << phase.first << "_s=" << phase.second;
    total += phase.second;
  }
  std::cerr << "\ttotal_s=" << total << std::endl;
  return 0;
}

template <bool hashes>
int run_query2(const fmsi2::index &index, std::string query_fn, unsigned batch_size, bool streaming, unsigned threads,
               size_t block_bases, fmsi2::query_stats &stats) {
  gzFile fp = OpenFile(query_fn);
  gzbuffer(fp, 1 << 18);
  bool ok = fmsi2::run_queries<hashes>(index, fp, stdout, streaming, batch_size, threads, block_bases, stats);
  gzclose(fp);
  if (!ok) {
    std::cerr << "ERROR: Query sequences of 4 Gbp and more are not supported." << std::endl;
    return 1;
  }
  return 0;
}

int ms_query2(int argc, char *argv[], bool output_orders) {
  bool usage = false;
  int c;
  int k = 0;
  std::string fn;

  if (argc > 1 && std::string(argv[argc - 1]) != "-h") {
    fn = argv[argc - 1];
    argc--;
  }

  std::string query_fn = "-";
  int batch_size = 32;
  // -1: the default of the index (fmsi2::default_prefix_length).
  int prefix_length = -1;
  bool streaming = false;
  bool use_klcp = true;
  int threads = 1;
  // Blocks of 128 kB keep the memory of the queries small and are as fast as larger ones.
  long long block_bases = 1 << 17;
  fmsi2::layout layout;
  enum { OPTION_LAYOUT_EXTRA_SMALL = 256, OPTION_LAYOUT_SMALL };
  static const struct option long_options[] = {
      {"extra-small", no_argument, nullptr, OPTION_LAYOUT_EXTRA_SMALL},
      {"small", no_argument, nullptr, OPTION_LAYOUT_SMALL},
      {nullptr, 0, nullptr, 0},
  };
  while ((c = getopt_long(argc, argv, "hq:k:b:p:SLt:B:c:m:O", long_options, nullptr)) >= 0) {
    switch (c) {
    case OPTION_LAYOUT_EXTRA_SMALL:
      layout = fmsi2::preset_layout("extra-small");
      break;
    case OPTION_LAYOUT_SMALL:
      layout = fmsi2::preset_layout("small");
      break;
    case 'h':
      usage = true;
      break;
    case 'q':
      query_fn = optarg;
      break;
    case 'k':
      k = atoi(optarg);
      break;
    case 'S':
      streaming = true;
      break;
    case 'L':
      use_klcp = false;
      break;
    case 't':
      threads = atoi(optarg);
      break;
    case 'B':
      block_bases = atoll(optarg);
      break;
    case 'b':
      batch_size = atoi(optarg);
      break;
    case 'p':
      prefix_length = atoi(optarg);
      break;
    case 'c':
    case 'm':
      if (!parse_layout_option(c, optarg, layout)) return usage_query2(output_orders);
      break;
    case 'O':
      // The option of FMSI 0.4 for max-one masks: not needed, as any mask is tested in constant time.
      break;
    default:
      return usage_query2(output_orders);
    }
  }
  if (usage) {
    usage_query2(output_orders);
    return 0;
  } else if (fn.empty()) {
    std::cerr << "ERROR: The index file is a required argument." << std::endl;
    return usage_query2(output_orders);
  } else if (batch_size < 1 || threads < 1 || block_bases < 1 || prefix_length < -1 || prefix_length > 14) {
    std::cerr << "ERROR: The batch size, the block size and the number of threads must be positive and the table pattern length between 0 and 14." << std::endl;
    return usage_query2(output_orders);
  }

  auto start = std::chrono::steady_clock::now();
  auto seconds_since_start = [&]() {
    return std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count();
  };

  const std::string index_path = fmsi2::find_index_file(fn);
  if (index_path.empty()) {
    std::cerr << "ERROR: No index file " << fn << " (nor " << fn << fmsi2::INDEX_EXTENSION << ")." << std::endl;
    return 1;
  }
  fmsi2::index index;
  fmsi2::load(index, index_path, threads, layout);
  if (output_orders && index.mask_type_name != "min-one" && index.mask_type_name != "min-one,max-one") {
    std::cerr << "Note: the mask of the index is " << (index.mask_type_name.empty() ? "of an unknown type" : index.mask_type_name)
              << ", not min-one: the hashes are unique but may not be minimal (a k-mer set more than once has several)." << std::endl;
  }
  if (k != 0 && k != (int)index.k) {
    std::cerr << "ERROR: Mismatch. Provided k (" << k << ") does not match the k of the index (" << index.k << ")." << std::endl;
    return usage_query2(output_orders);
  }
  if (streaming && use_klcp && std::filesystem::exists(index_path + fmsi2::KLCP_EXTENSION)) {
    fmsi2::load_klcp(index, index_path + fmsi2::KLCP_EXTENSION, threads);
  }
  double load_seconds = seconds_since_start();
  if (prefix_length < 0) prefix_length = fmsi2::default_prefix_length(index);
  fmsi2::build_prefix_table(index, prefix_length);
  double table_seconds = seconds_since_start() - load_seconds;
  // Resident memory once the index is ready, and the peak of loading it.
  const uint64_t loaded_bytes = fmsi2::resident_bytes(), load_peak_bytes = fmsi2::peak_resident_bytes();

  fmsi2::query_stats stats;
  int ret = output_orders ? run_query2<true>(index, query_fn, batch_size, streaming, threads, block_bases, stats)
                          : run_query2<false>(index, query_fn, batch_size, streaming, threads, block_bases, stats);
  fflush(stdout);

  // The times of the stages run by the searching threads are averages per thread.
  const fmsi2::lookup_stats &lookups = stats.lookups;
  std::cerr << "[fmsi2]"
            << "\tkmers=" << lookups.kmers << "\tfound=" << lookups.found
            << "\tsteps_per_kmer=" << (lookups.kmers ? (double)lookups.steps / lookups.kmers : 0)
            << "\tsearches=" << lookups.searches << "\trechecks=" << lookups.rechecks
            << "\tklcp=" << index.has_klcp() << "\tstrands=" << index.strands << "\tthreads=" << threads
            << "\tbwt_block_bytes=" << index.bwt.block_bytes() << "\tmask=" << fmsi2::mask_name(index.mask_type)
            << "\tindex_bytes_in_memory=" << index.bytes()
            << "\tloaded_rss_kb=" << loaded_bytes / 1000 << "\tload_peak_rss_kb=" << load_peak_bytes / 1000
            << "\tend_rss_kb=" << fmsi2::resident_bytes() / 1000 << "\tpeak_rss_kb=" << fmsi2::peak_resident_bytes() / 1000
            << "\ttable_length=" << index.prefix_length << "\ttable_bytes=" << index.prefix_table_bytes()
            << "\tload_s=" << load_seconds << "\ttable_s=" << table_seconds
            << "\tread_s=" << stats.read << "\tparse_s=" << stats.parse / threads
            << "\tscan_s=" << stats.scan / threads
            << "\tsearch_s=" << stats.search / threads
            << "\tformat_s=" << stats.format / threads << "\twrite_s=" << stats.write
            << "\tsearch_ns_per_kmer=" << (lookups.kmers ? 1e9 * stats.search / threads / lookups.kmers : 0)
            << "\ttotal_s=" << seconds_since_start() << std::endl;
  return ret;
}

int ms_merge(int argc, char *argv[]) {
  bool usage = false;
  int c;
  std::vector<std::string> fns;
  std::string result_fn;
  bool k_param_set = false;
  while ((c = getopt(argc, argv, "p:hr:k:")) >= 0) {
    switch (c) {
    case 'h':
      usage = true;
      break;
    case 'p':
      fns.emplace_back(optarg);
      break;
    case 'r':
      result_fn = optarg;
      break;
    case 'k':
        k_param_set = true;
        break;
    default:
      return usage_merge();
    }
  }
  if (usage) {
    usage_merge();
    return 0;
  }
  if (fns.size() < 2) {
    std::cerr << "At least two indices are required for merging." << std::endl;
    return usage_merge();
  }
  if (result_fn.empty()) {
      std::cerr << "Path to the result file is a required argument." << std::endl;
      return usage_merge();
  }
  if (k_param_set) {
      std::cerr << "WARNING: Parameter -k is ignored." << std::endl;
  }

  fms_index res = load_index(fns[0]);
  std::cerr << "Loaded index " << fns[0] << std::endl;

  for (size_t i = 1; i < fns.size(); ++i) {
      // TODO: investigate why this is needed (the ranks must have gotten set to nullptr somewhere).
      res.ac_gt_rank = sdsl::rank_support_v5<1>(&res.ac_gt);
      res.ac_rank = sdsl::rank_support_v5<1>(&res.ac);
      res.gt_rank = sdsl::rank_support_v5<1>(&res.gt);
      auto current = load_index(fns[i]);
      if (res.k != current.k) {
          std::cerr << "Mismatch. The k of the index " << fns[i] << " (" << current.k << ") does not match the k of the index " << fns[0] << "(" << res.k << ")." << std::endl;
          return usage_merge();
      }
      res = merge(res, current);
      std::cerr << "Loaded and merged index " << fns[i] << std::endl;
  }

  dump_index(res, result_fn);
  std::cerr << "Result written" << std::endl;
  return 0;
}

int ms_normalize(int argc, char *argv[]) {
  bool usage = false;
  int c;
  int k = 0;
  bool only_print = false;
  std::string fn;
  if (argc > 1 && std::string(argv[argc - 1]) != "-h") {
    fn = argv[argc - 1];
    argc--;
  }

  std::function<bool(size_t, size_t)> f = mask_function("or", true);
  while ((c = getopt(argc, argv, "hk:f:s")) >= 0) {
    switch (c) {
    case 'f':
      try {
        f = mask_function(optarg, true);
      } catch (std::invalid_argument &) {
        std::cerr << "Function '" << optarg << "' not recognized." << std::endl;
        return usage_query();
      }
      break;
    case 'h':
      usage = true;
      break;
    case 'p':
      fn = optarg;
      break;
    case 'k':
      k = atoi(optarg);
      break;
    case 's':
      only_print = true;
      break;
    default:
      return usage_normalize();
    }
  }
  if (usage) {
    usage_normalize();
    return 0;
  }
  if (fn.empty()) {
    std::cerr << "Path to the fasta file is a required argument." << std::endl;
    return usage_normalize();
  }

  std::cerr << "Starting " << fn << std::endl;
  fms_index index = load_index(fn);
    if (k != 0 && k != index.k) {
        std::cerr << "Provided k-mer size (" << k << ") does not match the one of the index (" << index.k << ")." << std::endl;
        return usage_normalize();
    }
  std::cerr << "Loaded index" << std::endl;
  auto ms = export_ms(index);
  ms = normalize(ms, index.k, f);
  std::cerr << "Compacted" << std::endl;
  if (only_print) {
    std::cout << ">exported f-masked superstring" << std::endl;
    std::cout << ms << std::endl;
    return 0;
  }
    if (index.k <= 32) index = construct<uint64_t>(ms, index.k, index.klcp.size() > 0);
    else index = construct<__uint128_t>(ms, index.k, index.klcp.size() > 0);
  dump_index(index, fn);
  std::cerr << "Written index" << std::endl;
  return 0;
}

int ms_export(int argc, char *argv[]) {
    bool usage = false;
    int c;
    std::string fn;
    if (argc > 1 && std::string(argv[argc - 1]) != "-h") {
      fn = argv[argc - 1];
      argc--;
    }

    while ((c = getopt(argc, argv, "h")) >= 0) {
        switch (c) {
            case 'h':
                usage = true;
                break;
            default:
                return usage_export();
        }
    }
    if (usage) {
        usage_export();
        return 0;
    }

    fms_index index = load_index(fn);

    if (index.sa_transformed_mask.size() == 0) {
      std::cerr << "ERROR: index not correctly loaded. Ensure that you correctly call `fmsi index` before." << std::endl;
      return usage_export();
    }

    std::cout << ">exported masked superstring" << std::endl;
    std::cout << export_ms(index) << std::endl;
    return 0;
}

int ms_op(int argc, char *argv[], std::string op) {
    bool usage = false;
    int c;
    int k = 0;
    std::vector<std::string> fns;
    std::string result_fn;
    while ((c = getopt(argc, argv, "p:hr:k:")) >= 0) {
        switch (c) {
            case 'h':
                usage = true;
                break;
            case 'p':
                fns.emplace_back(optarg);
                break;
            case 'r':
                result_fn = optarg;
                break;
            case 'k':
                k = atoi(optarg);
                break;
            default:
                return usage_merge();
        }
    }
    if (usage) {
        usage_op(op);
        return 0;
    }
    if (fns.size() < 2) {
        std::cerr << "At least two indices are required for merging." << std::endl;
        return usage_op(op);
    }
    if (result_fn.empty()) {
        std::cerr << "Path to the result file is a required argument." << std::endl;
        return usage_op(op);
    }
    fms_index res = load_index(fns[0]);
    std::cerr << "Loaded index " << fns[0] << std::endl;
    if (k != 0 && k != res.k) {
        std::cerr << "Provided k (." << k << ") does not match the k of the first index (" << res.k << ")." << std::endl;
        return usage_op(op);
    }


    for (size_t i = 1; i < fns.size(); ++i) {
        // TODO: investigate why this is needed (the ranks must have gotten set to nullptr somewhere).
        res.ac_gt_rank = sdsl::rank_support_v5<1>(&res.ac_gt);
        res.ac_rank = sdsl::rank_support_v5<1>(&res.ac);
        res.gt_rank = sdsl::rank_support_v5<1>(&res.gt);
        res = merge(res, load_index(fns[i]));
        if (op == "diff") {
            // TODO: investigate why this is needed (the ranks must have gotten set to nullptr somewhere).
            res.ac_gt_rank = sdsl::rank_support_v5<1>(&res.ac_gt);
            res.ac_rank = sdsl::rank_support_v5<1>(&res.ac);
            res.gt_rank = sdsl::rank_support_v5<1>(&res.gt);
            res = merge(res, load_index(fns[i]));
        }
        std::cerr << "Loaded and merged index " << fns[i] << std::endl;
    }


    // TODO: investigate why this is needed (the ranks must have gotten set to nullptr somewhere).
    res.ac_gt_rank = sdsl::rank_support_v5<1>(&res.ac_gt);
    res.ac_rank = sdsl::rank_support_v5<1>(&res.ac);
    res.gt_rank = sdsl::rank_support_v5<1>(&res.gt);

    auto ms = export_ms(res);

    demasking_function_t function = nullptr;
    if (op == "union") function =  mask_function("or", true);
    else if (op == "symdiff") function = mask_function("xor", true);
    else if (op == "diff") function = mask_function("1-1", true);
    else if (op == "inter") function = mask_function(std::to_string(fns.size()) + "-" + std::to_string(fns.size()), true);
    assert(function != nullptr);

    ms = normalize(ms, res.k, function);
    std::cerr << "Compacted result" << std::endl;

    if (res.k <= 32) res = construct<uint64_t>(ms, res.k, res.klcp.size() > 0);
    else res = construct<__uint128_t>(ms, res.k, res.klcp.size() > 0);

    dump_index(res, result_fn);
    std::cerr << "Result written" << std::endl;
    return 0;
}

int ms_clean(int argc, char *argv[]) {
  bool usage = false;
  int c;
  std::string fn;
  if (argc > 1 && std::string(argv[argc - 1]) != "-h") {
    fn = argv[argc - 1];
    argc--;
  }

  while ((c = getopt(argc, argv, "h")) >= 0) {
    switch (c) {
    case 'h':
      usage = true;
      break;
    default:
      return usage_clean();
    }
  }
  if (usage) {
    usage_clean();
    return 0;
  }
  if (fn == "") {
    std::cerr << "Path to the fasta file is a required argument." << std::endl;
    return usage_clean();
  }

  std::filesystem::remove(fn + ".fmsi.ac_gt");
  std::filesystem::remove(fn + ".fmsi.ac");
  std::filesystem::remove(fn + ".fmsi.gt");
  std::filesystem::remove(fn + ".fmsi.mask");
  std::filesystem::remove(fn + ".fmsi.misc");
  if (std::filesystem::exists(fn + ".fmsi.klcp")) {
    std::filesystem::remove(fn + ".fmsi.klcp");
  }
  std::cerr << "Cleaned " << fn << std::endl;
  return 0;
}

static int usage_access2() {
  std::cerr << std::endl;
  std::cerr << "Usage:   fmsi access [options] <index.fmsi>" << std::endl << std::endl;
  std::cerr << "Dictionary queries: prints the k-mer of every hash, as it occurs in the masked superstring; the inverse of" << std::endl;
  std::cerr << "`fmsi lookup` (for a min-one mask, of its minimal perfect hash). A value that is not a hash gives a line with `-`." << std::endl << std::endl;
  std::cerr << "Options:" << std::endl;
  std::cerr << "  -q FILE - Hashes, one non-negative integer per line; if `-` or not provided then standard input is used." << std::endl;
  std::cerr << "  -b INT  - Number of k-mers retrieved at the same time, with prefetching [default: 32]" << std::endl;
  std::cerr << std::endl;
  return 1;
}

int ms_access2(int argc, char *argv[]) {
  bool usage = false;
  int c;
  std::string fn;
  if (argc > 1 && std::string(argv[argc - 1]) != "-h") {
    fn = argv[argc - 1];
    argc--;
  }
  std::string query_fn = "-";
  int batch_size = 32;
  while ((c = getopt(argc, argv, "hq:b:")) >= 0) {
    switch (c) {
    case 'h':
      usage = true;
      break;
    case 'q':
      query_fn = optarg;
      break;
    case 'b':
      batch_size = atoi(optarg);
      break;
    default:
      return usage_access2();
    }
  }
  if (usage) {
    usage_access2();
    return 0;
  } else if (fn.empty()) {
    std::cerr << "ERROR: The index file is a required argument." << std::endl;
    return usage_access2();
  } else if (batch_size < 1) {
    std::cerr << "ERROR: The batch size must be positive." << std::endl;
    return usage_access2();
  }

  const std::string index_path = fmsi2::find_index_file(fn);
  if (index_path.empty()) {
    std::cerr << "ERROR: No index file " << fn << " (nor " << fn << fmsi2::INDEX_EXTENSION << ")." << std::endl;
    return 1;
  }
  fmsi2::index index;
  // Access selects on the mask, which only the plain encoding supports directly.
  fmsi2::layout layout;
  layout.mask = (int)fmsi2::mask_kind::plain;
  fmsi2::load(index, index_path, 1, layout);
  fmsi2::build_access_support(index);

  std::ifstream file;
  if (query_fn != "-") {
    file.open(query_fn);
    if (!file) {
      std::cerr << "ERROR: Cannot open " << query_fn << std::endl;
      return 1;
    }
  }
  std::istream &in = query_fn == "-" ? std::cin : file;

  // The hashes are processed in chunks; a value that is not a hash is replaced by a valid one and its line by `-`.
  const size_t chunk = 1 << 16, k = index.k;
  std::vector<uint64_t> hashes;
  std::vector<bool> valid;
  std::string kmers, output, line;
  uint64_t total = 0, invalid = 0;
  double seconds = 0;
  bool more = true;
  while (more) {
    hashes.clear();
    valid.clear();
    while (hashes.size() < chunk && (more = (bool)std::getline(in, line))) {
      char *end = nullptr;
      errno = 0;
      unsigned long long value = strtoull(line.c_str(), &end, 10);
      bool is_hash = !line.empty() && isdigit((unsigned char)line[0]) && errno == 0 && value < index.mask_ones &&
                     (*end == '\0' || isspace((unsigned char)*end));
      hashes.push_back(is_hash ? value : 0);
      valid.push_back(is_hash);
      invalid += !is_hash;
    }
    if (index.mask_ones == 0) std::fill(valid.begin(), valid.end(), false);
    kmers.assign(hashes.size() * k, 0);
    auto start = std::chrono::steady_clock::now();
    if (index.mask_ones > 0) fmsi2::access_batch(index, hashes.data(), hashes.size(), kmers.data(), batch_size);
    seconds += std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count();
    output.clear();
    for (size_t i = 0; i < hashes.size(); ++i) {
      if (valid[i]) {
        output.append(kmers, i * k, k);
      } else {
        output.push_back('-');
      }
      output.push_back('\n');
    }
    fwrite(output.data(), 1, output.size(), stdout);
    total += hashes.size();
  }
  fflush(stdout);
  std::cerr << "[fmsi2-access]"
            << "\thashes=" << total << "\tinvalid=" << invalid << "\tnumber_of_hashes=" << index.mask_ones
            << "\taccess_ns_per_kmer=" << (total ? 1e9 * seconds / total : 0) << std::endl;
  return 0;
}

/// Load an index of any supported format and write it as a .fmsi file (and its kLCP array, if it has one).
int ms_convert(int argc, char *argv[]) {
  int threads = 1;
  std::string output;
  int c;
  auto usage = []() {
    std::cerr << std::endl << "Usage:   fmsi convert [-t threads] [-o index.fmsi] <index>" << std::endl << std::endl;
    std::cerr << "Rewrites an index of an older format of FMSI 2 (<prefix>.fmsi2) as <prefix>.fmsi, with its kLCP array." << std::endl << std::endl;
    return 1;
  };
  while ((c = getopt(argc, argv, "t:o:h")) >= 0) {
    if (c == 't') {
      threads = atoi(optarg);
    } else if (c == 'o') {
      output = optarg;
    } else {
      return usage();
    }
  }
  if (optind != argc - 1 || threads < 1) return usage();
  const std::string name = argv[optind], path = fmsi2::find_index_file(name);
  if (path.empty()) {
    std::cerr << "ERROR: No index file " << name << "." << std::endl;
    return 1;
  }
  if (output.empty()) {
    const std::string old = fmsi2::OLD_INDEX_EXTENSION;
    const bool is_old = path.size() > old.size() && path.compare(path.size() - old.size(), old.size(), old) == 0;
    output = is_old ? path.substr(0, path.size() - old.size()) : path;
  }
  output = index_file_name(output);
  fmsi2::index index;
  fmsi2::load(index, path, threads);
  if (std::filesystem::exists(path + fmsi2::KLCP_EXTENSION)) fmsi2::load_klcp(index, path + fmsi2::KLCP_EXTENSION, threads);
  // The preset decides the layout that the new file recommends, as fmsi index would.
  const fmsi2::layout recommended = fmsi2::preset_layout(index.preset);
  const std::string temporary = output + ".convert";
  fmsi2::save(index, temporary, threads, &recommended);
  std::filesystem::rename(temporary, output);
  if (index.has_klcp()) {
    fmsi2::save_klcp(index, temporary, threads);
    std::filesystem::rename(temporary, output + fmsi2::KLCP_EXTENSION);
  }
  std::cerr << "Written " << output << (index.has_klcp() ? " and its kLCP array" : "") << std::endl;
  return 0;
}

/// Print the masked superstring of an index as FASTA.
int ms_export2(int argc, char *argv[]) {
  if (argc != 2 || std::string(argv[1]) == "-h") {
    std::cerr << std::endl << "Usage:   fmsi export <index.fmsi>" << std::endl << std::endl;
    std::cerr << "Prints the masked superstring of the index as FASTA: the superstring, in upper case where the mask is" << std::endl;
    std::cerr << "set, under the header of the masked superstring that was indexed. `fmsi index` gives the same index back." << std::endl << std::endl;
    return argc == 2 ? 0 : 1;
  }
  const std::string name = argv[1], path = fmsi2::find_index_file(name);
  if (path.empty()) {
    std::cerr << "ERROR: No index file " << name << "." << std::endl;
    return 1;
  }
  fmsi2::index index;
  fmsi2::load(index, path);
  const std::string superstring = fmsi2::export_masked_superstring(index);
  std::cout << ">" << (index.superstring_header.empty() ? "maskedsuperstring k=" + std::to_string(index.k) : index.superstring_header)
            << "\n" << superstring << "\n";
  return 0;
}

/// Print what an index file holds: its metadata, its sections, and the kLCP file next to it.
int ms_info(int argc, char *argv[]) {
  bool check = false;
  static const struct option long_options[] = {{"check", no_argument, nullptr, 'C'}, {nullptr, 0, nullptr, 0}};
  int c;
  auto usage = []() {
    std::cerr << std::endl << "Usage:   fmsi info [--check] <index.fmsi>" << std::endl << std::endl;
    std::cerr << "Prints what the index file holds and how: the description of the index, the sections of the file with" << std::endl;
    std::cerr << "their encodings and sizes, and the kLCP file next to it. --check verifies the checksums." << std::endl << std::endl;
    return 1;
  };
  while ((c = getopt_long(argc, argv, "h", long_options, nullptr)) >= 0) {
    if (c == 'C') {
      check = true;
    } else {
      return usage();
    }
  }
  if (optind != argc - 1) return usage();
  const std::string name = argv[optind], path = fmsi2::find_index_file(name);
  if (path.empty()) {
    std::cerr << "ERROR: No index file " << name << "." << std::endl;
    return 1;
  }
  if (fmsi2::file::file_version(path) != fmsi2::file::CONTAINER_VERSION) {
    std::cout << "file\t" << path << std::endl;
    std::cout << "format\t" << fmsi2::file::file_version(path) << " (the prototype of FMSI 2; `fmsi convert` rewrites it)" << std::endl;
    return 0;
  }
  auto print = [&](const std::string &file_path) {
    fmsi2::file::reader in(file_path);
    const fmsi2::file::metadata meta = in.meta();
    std::cout << "file\t" << file_path << std::endl;
    std::cout << "format\t" << fmsi2::file::CONTAINER_VERSION << std::endl;
    std::cout << "bytes\t" << in.file_bytes() << std::endl;
    for (const auto &item : meta.items()) std::cout << item.first << "\t" << item.second << std::endl;
    if (meta.has("rows") && meta.get("kind") == "index") {
      const double strands = (double)meta.number("strands", 1);
      // The superstring has (rows - strands) / strands characters; the rows include the sentinels.
      const double characters = (meta.number("rows") - strands) / strands;
      std::cout << "superstring_characters\t" << (uint64_t)characters << std::endl;
      std::cout << "bits_per_superstring_character\t" << 8.0 * in.file_bytes() / characters << std::endl;
      if (meta.get("mask_type").find("min-one") == 0) {
        std::cout << "kmers\t" << meta.number("mask_ones") / (uint64_t)strands << std::endl;
        std::cout << "bits_per_kmer\t" << 8.0 * in.file_bytes() * strands / meta.number("mask_ones") << std::endl;
      }
    }
    for (const auto &e : in.sections()) {
      std::cout << "section\t" << e.name() << "\tencoding=" << e.encoding << "\tflags=" << ((e.flags & fmsi2::file::REQUIRED) ? "required" : "optional")
                << ((e.flags & fmsi2::file::DEFLATE) ? ",deflate" : "") << "\titems=" << e.items << "\tbytes_in_file=" << e.stored_bytes
                << "\tbytes=" << e.raw_bytes;
      if (check) {
        in.read(e, 1, [](uint64_t, const uint8_t *, size_t) {}, true);
        std::cout << "\tchecksum=ok";
      }
      std::cout << std::endl;
    }
  };
  print(path);
  if (std::filesystem::exists(path + fmsi2::KLCP_EXTENSION) && fmsi2::file::file_version(path + fmsi2::KLCP_EXTENSION) != 0) {
    std::cout << std::endl;
    print(path + fmsi2::KLCP_EXTENSION);
  }
  return 0;
}

/// The commands of FMSI 0.4.
static int legacy(int argc, char *argv[]) {
  if (argc < 2) return usage_legacy();
  std::string op = argv[1];
  if (op == "index") return ms_index(argc - 1, argv + 1);
  if (op == "query") return ms_query(argc - 1, argv + 1, false);
  if (op == "lookup") return ms_query(argc - 1, argv + 1, true);
  if (op == "clean") return ms_clean(argc - 1, argv + 1);
  if (op == "merge") return ms_merge(argc - 1, argv + 1);
  // Recognize "normalize" for backwards compatibility
  if (op == "normalize" || op == "compact") return ms_normalize(argc - 1, argv + 1);
  if (op == "export") return ms_export(argc - 1, argv + 1);
  if (op == "union" || op == "inter" || op == "diff" || op == "symdiff") return ms_op(argc, argv, op);
  if (op == "-h") {
    usage_legacy();
    return 0;
  }
  return usage_legacy();
}

int main(int argc, char *argv[]) {
  if (argc < 2)
    return usage();
  std::string op = argv[1];
  try {
    // The names of the prototype (index2, ...) are kept as synonyms.
    if (op == "index" || op == "index2") return ms_index2(argc - 1, argv + 1);
    if (op == "query" || op == "query2") return ms_query2(argc - 1, argv + 1, false);
    if (op == "lookup" || op == "lookup2") return ms_query2(argc - 1, argv + 1, true);
    if (op == "access" || op == "access2") return ms_access2(argc - 1, argv + 1);
    if (op == "info") return ms_info(argc - 1, argv + 1);
    if (op == "export") return ms_export2(argc - 1, argv + 1);
    if (op == "convert") return ms_convert(argc - 1, argv + 1);
    if (op == "legacy") return legacy(argc - 1, argv + 1);
  } catch (const std::exception &error) {
    std::cerr << "ERROR: " << error.what() << std::endl;
    return 1;
  }
  if (op == "-v") return version();
  if (op == "-h") {
    usage();
    return 0;
  }
  return usage();
}
