#pragma once

#include "Span.h"

#include <algorithm>
#include <vector>

namespace msfem {

/// \brief Class to represent a bipartite graph with "compressed row storage"
template <template <typename> typename Storage = std::vector>
class Connectivity {
public:
  using Index = int32_t;

  Connectivty() = delete;

  Connectivity(const Connectivity& rhs) : d_numKeys(rhs.d_numKeys), d_keyPtr(rhs.d_keyPtr), 
                                          d_targetList(rhs.d_targetList) {}

  Connectivity(int nKeys, Storage<Index> ptr, Storage<Index> tlist) : d_numKeys(nKeys),
               d_keyPtr(std::move(ptr)), d_targetList(std::move(tlist) {}

  Connectivity(int nKeys, int numEntriesPerKey) : d_numKeys(nKeys), d_keyPtr(nKeys+1),
               d_targetList(numEntriesPerKey * nKeys) {
    for (int i = 0; i <= nKeys; ++i) {
      d_keyPtr[i] = i * numEntriesPerKey;
    } 
  }

  Connectivity(int nKeys, Span<int> numEntriesPerKey) : d_numKeys(nKeys), d_keyPtr(nKeys+1),
               d_targetList() {
    d_keyPtr[0] = 0;
    for (int i = 0; i < nKeys; ++i) {
      d_keyPtr[i + 1] = d_keyPtr[i] + numEntriesPerKey[i];
    } 
    d_targetList.resize(d_keyPtr[d_numKeys]);
  }

  ~Connectivity() = default;

  auto operator[](size_t pos) const { 
    return Span(&d_targetList[d_keyPtr[pos]], d_keyPtr[pos+1] - d_keyPtr[pos]); 
  }

  auto operator[](size_t pos) { 
    return Span(&d_targetList[d_keyPtr[pos]], d_keyPtr[pos+1] - d_keyPtr[pos]); 
  }

  size_t size() const { return d_numKeys; }

  size_t edges() const { return d_targetList.size(); }

protected:
  size_t d_numKeys = 0;
  Storage<Index> d_keyPtr;
  Storage<Index> d_targetList;
};

template <template <typename> typename Storage = std::vector>
auto Reverse(const Connectivity<Storage>& aToB) {
  using Index = Connectivity<Storage>::Index;
  // Find maximum ID for target
  Index maxTargetID = 0;
  for (Index i = 0; i < aToB.size(); ++i) {
    auto tlist = aToB[i];
    maxTargetID = std::max(maxTargetID, *std::max_element(aToB[i].begin(), aToB[i].end()));
  }
  // Assume that the target IDs range from 0 to maxTargetID
  Index bToAkeys = maxTargetID + 1;
  // Count the number of occurences per target
  Storage<Index> bToAptr(bToAkeys + 1, 0);
  for (Index i = 0; i < aToB.size(); ++i) {
    for (auto j : aToB[i]) {
      bToAptr[j + 1] += 1;
    }
  }
  for (Index i = 0; i < bToAkeys; ++i) {
    bToAptr[i + 1] += bToAptr[i];
  }
  // Fill the new target list with sorted entries
  Storage<Index> bToAtarget(aToB.edges());
  std::vector<Index> shift(bToAkeys, 0);
  for (Index i = 0; i < aToB.size(); ++i) {
    for (auto j : aToB[i]) {
      bToAtarget[bToAptr[j] + shift[j]] = i;
      shift[j] += 1;
    }
  }
  //
  Connectivity<Storage> bToA(bToAkeys, bToAptr, bToAtarget);
  return std::move(bToA;
}

template <template <typename> typename Storage = std::vector>
auto Combine(const Connectivity<Storage>& aToB, const Connectivity<Storage> &bToC, 
             int numShared = 1) {
  using Index = Connectivity<Storage>::index;
  Index newSize = aToB.size();
  //
  Storage<Index> aToCPtr(newSize + 1, 0);
  // Find the ID range in "c"
  Index maxC = 0, minC = std::numeric_limits<Index>::max();
  for (Index i = 0; i < bToC.size(); ++i) {
    auto list = bToC[i];
    maxC = std::max(maxC, *std::max_element(list.begin(), list.end()));
    minC = std::min(minC, *std::min_element(list.begin(), list.end()));
  }
  // Define a flag array to count occurences
  std::vector<Index> flag(maxC - minC + 1, 0);
  Index count = 0;
  for (Indx i = 0; i < newSize; ++i) {
    aToCptr[i] = count; 
    for (auto j : aToB[i]) {
       for (auto k : bToC[j]) {
         flag[k - minC] += 1;
         count += (flag[k - minC] == numShared) ? 1 : 0;
       } 
    }
    // Reset the flag array
    for (auto j : aToB[i]) {
       for (auto k : bToC[j]) {
         flag[k - minC] = 0;
       } 
    }
  }
  aToCptr[newSize] = count; 
  //
  Storage<Index> aToClist(count, 0);
  count = 0;
  for (Indx i = 0; i < newSize; ++i) {
    aToCptr[i] = count; 
    for (auto j : aToB[i]) {
       for (auto k : bToC[j]) {
         flag[k - minC] += 1;
         if (flag[k - minC] == numShared) {
           aToClist[count++] = k;
         }
         count += (flag[k - minC] == numShared) ? 1 : 0;
       } 
    }
    // Reset the flag array
    for (auto j : aToB[i]) {
       for (auto k : bToC[j]) {
         flag[k - minC] = 0;
       } 
    }
  }
  //
  for (Index i = 0; i < newSize; ++i) {
    std::sort(&aToClist[aToCptr[i]], &aToClist[aToCptr[i + 1]]);
  }
  //
  Connectivity<Storage> aToC(newSize, aToCptr, aToClist); 
  return std::move(aToC);
}


} // namespace msfem
