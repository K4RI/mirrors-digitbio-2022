/* Open source implementation of the OPUS Miner algorithm which applies OPUS search for Filtered Top-k Association Discovery of Self-Sufficient Itemsets
** Copyright (C) 2012 Geoffrey I Webb
**
** This program is free software: you can redistribute it and/or modify
** it under the terms of the GNU General Public License as published by
** the Free Software Foundation, either version 3 of the License, or
** (at your option) any later version.
** 
** This program is distributed in the hope that it will be useful,
** but WITHOUT ANY WARRANTY; without even the implied warranty of
** MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
** GNU General Public License for more details.
** 
** You should have received a copy of the GNU General Public License
** along with this program.  If not, see <http://www.gnu.org/licenses/>.
**
** Please report any bugs to Geoff Webb <geoff.webb@monash.edu>
*/

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <limits>
#ifdef _WIN32
  #include <time.h>
#else
  #include <sys/times.h>
#endif

#ifdef _DEBUG
#define _CRTDBG_MAP_ALLOC
#include <crtdbg.h>
#endif

#include <unistd.h>
#include "opus_miner.h"
#include "globals.h"
#include "load_data.h"
#include "find_itemsets.h"
#include "print_itemsets.h"
#include "filter_itemsets.h"
#include "utils.h"
#include "fisher.h"

std::priority_queue<itemsetRec> itemsets;

void print_header(FILE *f, int argc, const char **argv) {
  fprintf(f,
    "OPUS Miner: Filtered Top-k Association Discovery of Self-Sufficient Itemsets\n"
    "Version 1.1\n"
    "Copyright (C) 2012-2013 Geoffrey I Webb\n"
    "This program comes with ABSOLUTELY NO WARRANTY. This is free software, \n"
    "and you are welcome to redistribute it under certain conditions.\n"
    "See the GNU General Public Licence <http://www.gnu.org/licenses/> for details.\n"
    "\n"
    "If you publish results obtained by using this software please cite\n"
    "  Webb, G.I. & Vreeken, J. (in-press) Efficient Discovery of the Most Interesting\n"
    "  Associations.  ACM Transactions on Knowledge Discovery from Data.\n"
    "\n"
    );
  for (int i = 1; i < argc; i++) {
    if (argv[i][0] == '-' && argv[i][1] == 'k' && i < argc - 1) {
      fprintf(f, "  %s %s\n", argv[i], argv[i+1]);
      ++i;
    }
    else fprintf(f, "  %s\n", argv[i]);
  }
}

int
#ifdef _WIN32
cdecl
#endif
 main(const int argc, const char **argv) {
#ifdef _DEBUG
  _CrtSetDbgFlag ( _CRTDBG_ALLOC_MEM_DF | _CRTDBG_LEAK_CHECK_DF );
#endif
  char const *inputFileName = NULL;
  char const *outputFileName = NULL;
  static const char *usageStr = "Usage: %s [-c] [-f] [-k <k>] [-l] [-r] <input file> <output file>\n";
  std::vector<itemsetRec> is;
  FILE *outf = NULL;

  if (argc < 3) {
    fprintf(stderr, usageStr, argv[0]);
    exit(1);
  }

  print_header(stdout, argc, argv);

  int i;

  for (i = 1; i < argc; i++) {
    if (argv[i][0] == '-') {
      switch (argv[i][1]) {
      case 'c':
        printClosures = true;
        break;
      case 'f':
        filter = false;
        break;
      case 'k':
        if (argv[i][2] == '\0') {
          k = getNum(argv[++i]);
        }
        else {
          k = getNum(argv[i]+2);
        }
        break;
      case 'l':
        searchByLift = true;
        break;
      case 'p':
        correctionForMultCompare = false; // do not correct alpha for the size of the search space
        break;
      case 'r':
        redundancyTests = false;
	      break;
      default: 
        fprintf(stderr, usageStr, argv[0]);
        exit(1);
      }
    }
    else if (inputFileName == NULL) {
      inputFileName = argv[i];
    }
    else if (outputFileName == NULL) {
      outputFileName = argv[i];
    }
    else {
      fprintf(stderr, usageStr, argv[0]);
      exit(1);
    }
  }

#ifdef _WIN32
  time_t start_t = time(NULL);
#else
  struct tms start_t;
  times(&start_t);
#endif // _WIN32

  try {
    printf("\nLoading data from %s\n", inputFileName);
    load_data(inputFileName);

    printf("%d transactions, %d items\n", noOfTransactions, noOfItems);

    outf = fopen(outputFileName, "w");

    if (outf == NULL) {
      fprintf(stderr, "Cannot open output file '%s'\n", outputFileName);
      exit(1);
    }

    print_header(outf, argc, argv);

    fprintf(outf, "\n%s: %ld items, %ld transactions\n", inputFileName, noOfItems, noOfTransactions);


    printf("Finding itemsets\n");
  #ifdef _WIN32
    time_t find_start_t = time(NULL);
  #else
    struct tms find_start_t;
    times(&find_start_t);
  #endif // _WIN32

    find_itemsets();

  #ifdef _WIN32
    time_t find_end_t = time(NULL);
  #else
    struct tms find_end_t;
    times(&find_end_t);
  #endif // _WIN32

    // extract the itemsets from the priority queue
    while (!itemsets.empty()) {
      is.push_back(itemsets.top());
      itemsets.pop();
    }

    if (filter) {
      printf("Filtering itemsets\n");
      filter_itemsets(is);
    }

  #ifdef _WIN32
    time_t print_start_t = time(NULL);

    const long tm = static_cast<long>(print_start_t-find_start_t);
  #else
    struct tms print_start_t;
    times(&print_start_t);

    const long tm = static_cast<long>(print_start_t.tms_utime-find_start_t.tms_utime) / sysconf(_SC_CLK_TCK);
  #endif // _WIN32

    fprintf(outf, "Found %1d non-redundant productive itemsets in %ld seconds\n", static_cast<long>(is.size()), tm);

    printf("Printing itemsets\n");
    print_itemsets(outf, is);

  #ifdef _WIN32
    time_t end_t = time(NULL);
  #else
    struct tms end_t;
    times(&end_t);
  #endif // _WIN32

  #ifdef _WIN32
    const long t =  static_cast<long>(end_t-start_t);

    printf("%ld seconds (%ld input, %ld search, %ld filter, %ld output)", t, static_cast<long>(find_start_t-start_t), static_cast<long>(find_end_t-find_start_t), static_cast<long>(print_start_t-find_end_t), static_cast<long>(end_t-print_start_t));
  #else
    const long ticks_sec = sysconf(_SC_CLK_TCK);
    const long t = static_cast<long>(end_t.tms_utime-start_t.tms_utime) / ticks_sec;

    printf("%ld seconds (%ld input, %ld search, %ld filter, %ld output)", t, static_cast<long>(find_start_t.tms_utime-start_t.tms_utime)/ticks_sec, static_cast<long>(find_end_t.tms_utime-find_start_t.tms_utime)/ticks_sec, static_cast<long>(print_start_t.tms_utime-find_end_t.tms_utime)/ticks_sec, static_cast<long>(end_t.tms_utime-print_start_t.tms_utime)/ticks_sec);
  #endif // _WIN32
    printf(" for %ld itemsets\n", static_cast<long>(is.size()));
  }
  catch (std::bad_alloc) {
    fprintf(stderr, "Out of memory\n");
  }
  catch (...) {
    fprintf(stderr, "Unhandled exception\n");
  }

  fclose(outf);
  return 0;
}
