#include <filesystem>
#include <sdsl/suffix_arrays.hpp>

#include "fms_index.h"
#include "fms_index2.h"
#include "doubled.h"
#include "query_pipeline.h"
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
  std::cerr << "Program: FMSI - a tool for space-efficient k-mer set indexing via masked superstrings." << std::endl;
  std::cerr << "Version: " << VERSION << std::endl;
  std::cerr << "Contact: Ondrej Sladky (ondra.sladky@gmail.com)" << std::endl << std::endl;
  std::cerr << "Usage:   fmsi <command> [options]" << std::endl << std::endl;
  std::cerr << "Command (stable):" << std::endl;
  std::cerr << "    index   - Creates a BWT based index of the given masked superstring." << std::endl;
  std::cerr << "    query   - Queries k-mers against an index." << std::endl;
  std::cerr << "    lookup  - Return unique hashes of present k-mers." << std::endl;
  std::cerr << "    export  - Print the underlying masked superstring to stdout." << std::endl << std::endl;
  std::cerr << "Command (FMSI 2 prototype):" << std::endl;
  std::cerr << "    index2  - Creates an index with a cache-friendly layout." << std::endl;
  std::cerr << "    query2  - Queries k-mers against an index2 index in batches." << std::endl;
  std::cerr << "    lookup2 - Return unique hashes of present k-mers using an index2 index." << std::endl;
  std::cerr << "    access2 - Return the k-mers of hashes using an index2 index." << std::endl << std::endl;
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
  std::cerr << "Usage:   fmsi index [options] <masked-superstring-input>" << std::endl << std::endl;
  std::cerr << "Options:" << std::endl;
  std::cerr << "    -k INT  - size of k-mers [recommended, default: number of mask trailing zeros - 1]"
            << std::endl;
  std::cerr << "    -x      - do not compute the kLCP array used for faster streaming queries."
            << std::endl;
  std::cerr << "    -Q      - construct the suffix array by QSufSort instead of libsais."
            << std::endl << std::endl;
  std::cerr << "Note: `fmsi index` accepts only masked superstrings - these can be computed e.g. by KmerCamel from any FASTA file." << std::endl << std::endl;
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
  std::cerr << "Usage:   fmsi index2 [options] <masked-superstring-input>" << std::endl << std::endl;
  std::cerr << "Presets:" << std::endl;
  std::cerr << "    --small - the masked superstring alone [default]" << std::endl;
  std::cerr << "    --fast  - the masked superstring and its reverse complement (twice the space): a k-mer is" << std::endl;
  std::cerr << "              searched as it is only, never as its reverse complement" << std::endl;
  std::cerr << "    Both use the plain mask and BWT blocks of 64 bytes, which -m and -c change." << std::endl << std::endl;
  std::cerr << "Options:" << std::endl;
  std::cerr << "    -k INT  - size of k-mers [default: number of mask trailing zeros - 1]" << std::endl;
  std::cerr << "    -m STR  - encoding of the mask: plain (the fastest), packed (smaller), or rrr (the smallest, slow) [default: plain]" << std::endl;
  std::cerr << "    -c INT  - bytes per block of the BWT: 64 (the fastest), 128 or 256 (the smallest) [default: 64]" << std::endl;
  std::cerr << "    -l      - also construct the kLCP array (1 bit per character), which speeds up streamed queries" << std::endl;
  std::cerr << "    -o STR  - prefix of the index file [default: path to the masked superstring]" << std::endl;
  std::cerr << "    -Q      - construct the suffix array by QSufSort instead of libsais (--small only)" << std::endl;
  std::cerr << "    -t INT  - number of threads [default: 1]; the suffix array is constructed on "
            << (SUFFIX_ARRAY_THREADED ? "all of them (compiled with OpenMP)" : "one thread (compiled without OpenMP)")
            << std::endl << std::endl;
  std::cerr << "The index is written to <prefix>.fmsi2 and the kLCP array to <prefix>.fmsi2.klcp." << std::endl << std::endl;
  return 1;
}

static int usage_query2(bool lookup) {
  std::cerr << std::endl;
  std::cerr << "Usage:   fmsi " << (lookup ? "lookup2" : "query2") << " [options] <index-prefix>" << std::endl << std::endl;
  std::cerr << "Options:" << std::endl;
  std::cerr << "  -q FILE - Path to FASTA/FASTQ with queries [default: stdin]" << std::endl;
  std::cerr << "  -k INT  - Size of k-mers [default: infer automatically from index]" << std::endl;
  std::cerr << "  -S      - Streamed queries: search the k-mers of a sequence by extending the match of the previous one;" << std::endl;
  std::cerr << "            uses the kLCP array if the index has one" << std::endl;
  std::cerr << "  -L      - Do not use the kLCP array" << std::endl;
  std::cerr << "  -b INT  - Number of k-mers (or sequences if -S) searched at the same time; 1 turns batching off [default: 32]" << std::endl;
  std::cerr << "  -p INT  - Length of the patterns with tabulated SA intervals; 0 turns the table off [default: 10]" << std::endl;
  std::cerr << "  -t INT  - Number of threads searching; with more than one, the queries are also read and the output printed" << std::endl;
  std::cerr << "            by a thread each [default: 1]" << std::endl;
  std::cerr << "  -B INT  - Size of a block of sequences handed to a thread, in bytes of the query file [default: 1048576]" << std::endl;
  std::cerr << std::endl;
  return 1;
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
  fmsi2::mask_kind mask_type = fmsi2::mask_kind::plain;
  int block_bytes = 64;
  bool use_klcp = false;
  int threads = 1;
  // Presets: --small (one strand) or --fast (both strands).
  bool fast = false;
  enum { OPTION_SMALL = 256, OPTION_FAST };
  static const struct option long_options[] = {
      {"small", no_argument, nullptr, OPTION_SMALL},
      {"fast", no_argument, nullptr, OPTION_FAST},
      {nullptr, 0, nullptr, 0},
  };
  while ((c = getopt_long(argc, argv, "hk:m:c:lo:Qt:", long_options, nullptr)) >= 0) {
    switch (c) {
    case OPTION_SMALL:
      fast = false;
      break;
    case OPTION_FAST:
      fast = true;
      break;
    case 'h':
      usage = true;
      break;
    case 'k':
      k = atoi(optarg);
      break;
    case 'c':
      block_bytes = atoi(optarg);
      if (!fmsi2::quad_sequence::valid_block_bytes(block_bytes)) {
        std::cerr << "ERROR: The blocks of the BWT can have 64, 128 or 256 bytes." << std::endl;
        return usage_index2();
      }
      break;
    case 'm':
      if (std::string(optarg) == "plain") {
        mask_type = fmsi2::mask_kind::plain;
      } else if (std::string(optarg) == "packed") {
        mask_type = fmsi2::mask_kind::packed;
      } else if (std::string(optarg) == "rrr") {
        mask_type = fmsi2::mask_kind::rrr;
      } else {
        std::cerr << "ERROR: Mask encoding '" << optarg << "' not recognized." << std::endl;
        return usage_index2();
      }
      break;
    case 'l':
      use_klcp = true;
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
  if (prefix.empty()) prefix = fn;
  if (fast && use_qsufsort) {
    std::cerr << "ERROR: QSufSort (-Q) constructs the --small index only." << std::endl;
    return usage_index2();
  }

  build_phases = phase_timer();
  auto ms = read_masked_superstring(fn);
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
  mbwt m;
  if (fast) {
    m = doubled::compute_mbwt(std::move(ms), k, use_klcp, threads);
  } else {
    m = k <= 32 ? compute_mbwt<uint64_t>(ms, k, use_klcp, use_qsufsort, threads)
                : compute_mbwt<__uint128_t>(ms, k, use_klcp, use_qsufsort, threads);
  }
  std::cerr << "Constructed the masked Burrows-Wheeler transform" << (fast ? " of both strands" : "") << std::endl;
  std::string().swap(ms);
  fmsi2::index index;
  fmsi2::build(index, m, k, mask_type, block_bytes, threads);
  fmsi2::save(index, prefix + ".fmsi2");
  std::cerr << "Written index: BWT " << index.bwt.bytes() << " bytes (" << 8.0 * index.bwt.bytes() / length
            << " bits per character of the superstring), mask " << index.mask_bytes() << " bytes ("
            << 8.0 * index.mask_bytes() / length << " bits per character of the superstring)" << std::endl;
  if (use_klcp) {
    fmsi2::save_klcp(index, prefix + ".fmsi2.klcp");
    std::cerr << "Written kLCP array: " << index.strands << " bit(s) per character of the superstring" << std::endl;
  }
  build_phases.lap("write");

  double total = 0;
  std::cerr << "[fmsi2-build]\tpreset=" << (fast ? "fast" : "small") << "\tstrands=" << index.strands
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
  int prefix_length = 10;
  bool streaming = false;
  bool use_klcp = true;
  int threads = 1;
  long long block_bases = 1 << 20;
  while ((c = getopt(argc, argv, "hq:k:b:p:SLt:B:")) >= 0) {
    switch (c) {
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
    default:
      return usage_query2(output_orders);
    }
  }
  if (usage) {
    usage_query2(output_orders);
    return 0;
  } else if (fn.empty()) {
    std::cerr << "ERROR: Prefix of the index is a required argument." << std::endl;
    return usage_query2(output_orders);
  } else if (batch_size < 1 || threads < 1 || block_bases < 1 || prefix_length < 0 || prefix_length > 14) {
    std::cerr << "ERROR: The batch size, the block size and the number of threads must be positive and the table pattern length between 0 and 14." << std::endl;
    return usage_query2(output_orders);
  }

  auto start = std::chrono::steady_clock::now();
  auto seconds_since_start = [&]() {
    return std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count();
  };

  fmsi2::index index;
  fmsi2::load(index, fn + ".fmsi2");
  if (k != 0 && k != (int)index.k) {
    std::cerr << "ERROR: Mismatch. Provided k (" << k << ") does not match the k of the index (" << index.k << ")." << std::endl;
    return usage_query2(output_orders);
  }
  if (streaming && use_klcp && std::filesystem::exists(fn + ".fmsi2.klcp")) {
    fmsi2::load_klcp(index, fn + ".fmsi2.klcp");
  }
  double load_seconds = seconds_since_start();
  fmsi2::build_prefix_table(index, prefix_length);
  double table_seconds = seconds_since_start() - load_seconds;

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
  std::cerr << "Usage:   fmsi access2 [options] <index-prefix>" << std::endl << std::endl;
  std::cerr << "Prints the k-mer of every hash, as it occurs in the masked superstring; the inverse of `fmsi lookup2`" << std::endl;
  std::cerr << "for a mask with the minimum number of ones. A value that is not a hash gives a line with `-`." << std::endl << std::endl;
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
    std::cerr << "ERROR: Prefix of the index is a required argument." << std::endl;
    return usage_access2();
  } else if (batch_size < 1) {
    std::cerr << "ERROR: The batch size must be positive." << std::endl;
    return usage_access2();
  }

  fmsi2::index index;
  fmsi2::load(index, fn + ".fmsi2");
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

int main(int argc, char *argv[]) {
  int ret;
  if (argc < 2)
    return usage();
  std::string op = argv[1];
  if (op == "index")
    ret = ms_index(argc - 1, argv + 1);
  else if (op == "query")
    ret = ms_query(argc - 1, argv + 1, false);
  else if (op == "lookup")
    ret = ms_query(argc - 1, argv + 1, true);
  else if (op == "index2")
    ret = ms_index2(argc - 1, argv + 1);
  else if (op == "query2")
    ret = ms_query2(argc - 1, argv + 1, false);
  else if (op == "lookup2")
    ret = ms_query2(argc - 1, argv + 1, true);
  else if (op == "access2")
    ret = ms_access2(argc - 1, argv + 1);
  else if (op == "clean")
    ret = ms_clean(argc - 1, argv + 1);
  else if (op == "merge")
    ret = ms_merge(argc - 1, argv + 1);
  // Recognize "normalize" for backwards compatibility
  else if (op == "normalize" || op == "compact")
    ret = ms_normalize(argc - 1, argv + 1);
  else if (op == "export")
    ret = ms_export(argc - 1, argv + 1);
  else if (op == "union" || op == "inter" || op == "diff" || op == "symdiff")
    ret = ms_op(argc, argv, op);
  else if (op == "-v")
    return version();
  else if (op == "-h") {
    usage();
    return 0;
  } else
    return usage();

  return ret;
}
