#include "Geometry/Partition.h"

#include <cstring>
#include <string>


using namespace FECode;


#ifdef FEC_USE_METIS

extern "C" {
void METIS_PartGraphKway(int *n, int *xadj, int *adjncy, int *vwgt, int *adjwgt,
     int *wgtflag, int *numflag, int *nparts, int *options, int *edgecut, int *part);
void METIS_PartGraphRecursive(int *n, int *xadj, int *adjncy, int *vwgt, int *adjwgt,
     int *wgtflag, int *numflag, int *nparts, int *options, int *edgecut, int *part);
}

void Partitioner::METISCellPartition(int numPart, int numCell, int *cellToPart,
                                     int *xadj, int *adjncy) {

  // Call METIS to distribute cells

  int zero = 0;
  int options[5];
  options[0] = 0; options[1] = 0; options[2] = 0; options[3] = 0; options[4] = 0;
  int edgecut = 0;

  memset(cellToPart, 0, numCell*sizeof(int));
  if (numPart <= 1)
    return;

  if (numPart <= 8) {
    METIS_PartGraphRecursive(&numCell, xadj, adjncy, NULL, NULL, &zero, &zero,
                        &numPart, &options[0], &edgecut, cellToPart);
  }
  else {
    METIS_PartGraphKway(&numCell, xadj, adjncy, NULL, NULL, &zero, &zero,
                        &numPart, &options[0], &edgecut, cellToPart);
  }

}


void Partitioner::METISPartitionRecursive(int numPart, int numCell, int *cellToPart,
                                     int *xadj, int *adjncy) {

  // Call METIS to distribute cells
  // Note: This algorithm seems to give nested partitions for successive powers of 2

  int zero = 0;
  int options[5];
  options[0] = 0; options[1] = 0; options[2] = 0; options[3] = 0; options[4] = 0;
  int edgecut = 0;

  memset(cellToPart, 0, numCell*sizeof(int));
  if (numPart <= 1)
    return;

  METIS_PartGraphRecursive(&numCell, xadj, adjncy, NULL, NULL, &zero, &zero,
                      &numPart, &options[0], &edgecut, cellToPart);

}


void Partitioner::METISPartitionKWay(int numPart, int numCell, int *cellToPart,
                                     int *xadj, int *adjncy) {

  // Call METIS to distribute cells

  int zero = 0;
  int options[5];
  options[0] = 0; options[1] = 0; options[2] = 0; options[3] = 0; options[4] = 0;
  int edgecut = 0;

  memset(cellToPart, 0, numCell*sizeof(int));
  if (numPart <= 1)
    return;

  METIS_PartGraphKway(&numCell, xadj, adjncy, NULL, NULL, &zero, &zero,
                      &numPart, &options[0], &edgecut, cellToPart);

}
#endif


void Partitioner::UniformCellPartition(int numPart, int numCell, int *cellToPart) {

  // Define a uniform partition of cells

  int remainder = numCell % numPart;
  int ratio = numCell / numPart;
  int flag = remainder * (ratio + 1);
  for (int iE = 0; iE < numCell; ++iE) {
    if (iE < flag)
      cellToPart[iE] = iE / (ratio + 1);
    else
      cellToPart[iE] = remainder + (iE - flag) / ratio;
  }

}

