#!/bin/bash
#
# Helper script to set MFEM ratio in source code
# Usage: ./set_ratio.sh <ratio>
# Example: ./set_ratio.sh 32
#

if [ $# -ne 1 ]; then
    echo "Usage: $0 <ratio>"
    echo "  ratio: 8, 16, 32, or 64"
    echo ""
    echo "Example:"
    echo "  $0 32          # Set ratio to 32"
    echo "  cd build && make cuda_assembly openmp_assembly"
    echo "  ./run_cpp_cuda_r32.sh"
    exit 1
fi

RATIO=$1

# Validate ratio
if [[ ! "$RATIO" =~ ^(8|16|32|64)$ ]]; then
    echo "Error: Invalid ratio '$RATIO'"
    echo "Valid ratios: 8, 16, 32, 64"
    exit 1
fi

SOURCE_FILE="src/ScaledLaplacian.h"
BACKUP_FILE="${SOURCE_FILE}.bak"

# Check if source file exists
if [ ! -f "$SOURCE_FILE" ]; then
    echo "Error: $SOURCE_FILE not found"
    exit 1
fi

# Create backup
cp "$SOURCE_FILE" "$BACKUP_FILE"
echo "Created backup: $BACKUP_FILE"

# Update ratio on line 144
# The line should be: "  static constexpr int ratio = XX;"
sed -i.tmp "144s/ratio = [0-9]\+/ratio = $RATIO/" "$SOURCE_FILE"
rm -f "${SOURCE_FILE}.tmp"

# Verify the change
CURRENT_RATIO=$(sed -n '144p' "$SOURCE_FILE" | grep -oP 'ratio = \K[0-9]+')

if [ "$CURRENT_RATIO" = "$RATIO" ]; then
    echo "✓ Successfully set ratio to $RATIO in $SOURCE_FILE (line 144)"
    echo ""
    echo "Next steps:"
    echo "  1. Recompile: cd build && make cuda_assembly openmp_assembly && cd .."
    echo "  2. Run tests: ./run_cpp_cuda_r${RATIO}.sh"
else
    echo "✗ Failed to update ratio (current: $CURRENT_RATIO, expected: $RATIO)"
    echo "Restoring backup..."
    mv "$BACKUP_FILE" "$SOURCE_FILE"
    exit 1
fi
