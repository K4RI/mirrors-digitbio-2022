/* print_itemsets.cpp - a module of OPUS Miner providing print_itemsets, a procedure to print the top-k itemsets to a file.
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
*/

#include <stdio.h>
#include <vector>
#include <algorithm>
#include "globals.h"
#include "itemset.h"
#include "utils.h"
#include "find_closure.h"

bool valgt(itemsetRec i1, itemsetRec i2) {
  return i1.value > i2.value;
}

void print_itemset(FILE *f, const itemset &is) {
  itemset::const_iterator item_it;

  for (item_it = is.begin(); item_it != is.end(); item_it++) {
    if (item_it != is.begin()) {
      fputc(',', f);
    }

    fprintf(f, "%s", itemNames[*item_it].c_str());
  }
}

void print_itemsetRec(FILE *f, const itemsetRec &is) {
  print_itemset(f, is);

  fprintf(f, " [%d,%f", is.count, is.value);
  fprintf(f, " %g]", is.p);

  if (printClosures) {
    itemset closure;

    find_closure(is, closure);

    if (closure.size() > is.size()) {
      fprintf(f, " closure: ");
      print_itemset(f, closure);
    }
  }

  putc('\n', f);
}

void print_itemsets(FILE *f, std::vector<itemsetRec> &is) {
  int i;
  
  for (i = 2; i < alpha.size(); i++) {
      fprintf(f, "Alpha for size %d = %g\n", i, alpha[i]);
  }

  fprintf(f, "\nSELF-SUFFICIENT ITEMSETS:\n");

  std::sort(is.begin(), is.end(), valgt);

  std::vector<itemsetRec>::const_iterator it;

  int failed_count = 0;

  for (it = is.begin(); it != is.end(); it++) {
    if (!it->self_sufficient) {
      failed_count++;
    }
    else {
      print_itemsetRec(f, *it);
    }
  }

  if (failed_count) {
    fprintf(f, "\n%d itemsets failed test for self sufficiency\n", failed_count);
    for (it = is.begin(); it != is.end(); it++) {
      if (!it->self_sufficient) {
        print_itemsetRec(f, *it);
      }
    }
  }
}
