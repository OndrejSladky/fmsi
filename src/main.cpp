#include <inttypes.h>
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

#include "index.h"
#include "bwa.h"
#include "mask.h"

//#include "bwt_gen.c"

static int usage() {
  fprintf(stderr, "HERE WILL BE USAGE\n");
  return 1;
}

static int usage_index() {
  fprintf(stderr, "HERE WILL BE USAGE\n");
  return 1;
}

int bwa_fa2pac(int argc, char *argv[]);
int bwa_pac2bwt(int argc, char *argv[]);
int bwt_bwtgen_main(int argc, char *argv[]);
int bwa_bwtupdate(int argc, char *argv[]);
int bwa_bwt2sa(int argc, char *argv[]);

int ms_index(int argc, char *argv[]) {
  int c;
  int usage = 0;
  int k = 13;
  while ((c = getopt(argc, argv, "k:h")) >= 0) {
    switch (c) {
    case 'k':
        k = atoi(optarg);
      break;
    case 'h':
      usage = 1;
      break;
    default:
      usage_index();
      return 1;
    }
  }
  if (usage) {
    usage_index();
    return 0;
  }
  if (optind + 1 > argc) {
    usage_index();
    return 1;
  }
  char *prefix = (char*)malloc(strlen(argv[optind]) * sizeof(char));
  strcpy(prefix, argv[optind]);

  // Construct and dump the FM-index.
  char *arguments[3];
  arguments[0] = (char*)malloc(10 * sizeof(char));
  arguments[1] = (char*)malloc((strlen(prefix) + 10) * sizeof(char));
  arguments[2] = (char*)malloc((strlen(prefix) + 10) * sizeof(char));
  /*
  strcpy(arguments[0], "fa2pac");
  strcpy(arguments[1], prefix);
  strcpy(arguments[2], prefix);
  optind = 1;
  bwa_fa2pac(3, arguments);
  strcpy(arguments[0], "pac2bwt");
  strcat(arguments[1], ".pac");
  strcat(arguments[2], ".bwt");
  optind = 1;
  bwt_bwtgen_main(3, arguments);
  strcpy(arguments[0], "bwtupdate");
  strcpy(arguments[1], prefix);
  strcat(arguments[1], ".bwt");
  optind = 1;
  bwa_bwtupdate(2, arguments);*/

  // Construct and dump the BW-transformed mask.
  mask_t bwTransformedMask = construct_bw_transformed_mask(prefix, k);
  strcpy(arguments[0], prefix);
  strcat(arguments[0], ".mask");
  //mask_dump(arguments[0], bwTransformedMask);
  free(prefix);
  return 0;
}

int main(int argc, char *argv[]) {
  int ret = 0;
  if (argc < 2) {
    usage();
    return 0;
  } else if (strcmp(argv[1], "index") == 0)
    ret = ms_index(argc - 1, argv + 1);
  else
    return usage();

  return ret;
}
