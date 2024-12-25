#ifndef FECODE_PARTITION_H
#define FECODE_PARTITION_H

namespace FECode {

  /// A class to implement partitioning strategies for meshes.

  /// This class implements partitioning strategies for meshes.
  /// It uses METIS on each processor when the code is linked to
  /// the METIS library.
  class Partitioner {
  
    public:
    
      #ifdef FEC_USE_METIS
      //! Routine to partition cells with METIS.
      /// \param[in] numPart: Number of partitions to be created.
      /// \param[in] numCell: Number of cells to partition.
      /// \param[out] cellToPart: Array whose length is the number of cells.
      /// On exit, the array stores the one-to-one map of a cell to a partition.
      /// \param[in] xadj: Array of integers...
      /// \param[in] adjacency: Array of integers...
      /// \note This routines requires the METIS library.
      /// \note When the number of partitions is smaller or equal to 8,
      /// the algorithm is METIS_PartGraphRecursive.
      /// \note When the number of partitions is greater than 8,
      /// the algorithm is METIS_PartGraphKWay.
      void METISCellPartition(int numPart, int numCell, int *cellToPart,
                              int *xadj, int *adjncy);

      //! Routine to partition cells with METIS.
      /// \param[in] numPart: Number of partitions to be created.
      /// \param[in] numCell: Number of cells to partition.
      /// \param[out] cellToPart: Array whose length is the number of cells.
      /// On exit, the array stores the one-to-one map of a cell to a partition.
      /// \param[in] xadj: ...
      /// \param[in] adjacency: ...
      /// \note This routines requires the METIS library.
      /// \note The algorithm is METIS_PartGraphRecursive.
      void METISPartitionRecursive(int numPart, int numCell, int *cellToPart,
                              int *xadj, int *adjncy);

      //! Routine to partition cells with METIS.
      /// \param[in] numPart: Number of partitions to be created.
      /// \param[in] numCell: Number of cells to partition.
      /// \param[out] cellToPart: Array whose length is the number of cells.
      /// On exit, the array stores the one-to-one map of a cell to a partition.
      /// \param[in] xadj: ...
      /// \param[in] adjacency: ...
      /// \note This routines requires the METIS library.
      /// \note The algorithm is METIS_PartGraphKWay.
      void METISPartitionKWay(int numPart, int numCell, int *cellToPart,
                              int *xadj, int *adjncy);
      #endif

      //! Routine to uniformly partition cells.
      /// \param[in] numPart: Number of partitions to be created.
      /// \param[in] numCell: Number of cells to partition.
      /// \param[out] cellToPart: Array whose length is the number of cells.
      /// On exit, the array stores the one-to-one map of a cell to a partition.
      /// \note The partition is based on the element numbering.
      /// It distributes uniformly the elements.
      /// When working on a mesh, there is no guarantee that the resulting
      /// subdomains are connected.
      void UniformCellPartition(int numPart, int numCell, int *cellToPart);

  };

}

#endif
