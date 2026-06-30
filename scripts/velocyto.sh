#!/bin/bash
#SBATCH --job-name=velocyto
#SBATCH -o velocyto_%j.out
#SBATCH -e velocyto_aggr_%j.err
#SBATCH --ntasks=1
#SBATCH --mem=64G
#SBATCH --cpus-per-task=4
#SBATCH --time=24:00:00
#SBATCH --mail-type=end
#SBATCH --mail-user=david.wragg@abdn.ac.uk

# For more options see: https://slurm.schedmd.com/sbatch.html

BARCODES=
OUT=
RPT=
BAM=
GENES=

# Parsing named arguments
while [[ $# -gt 0 ]]; do
  key="$1"
  case $key in
    --barcodes)
      BARCODES="$2"
      shift
      shift
      ;;
    --genes)
      GENES="$2"
      shift
      shift
      ;;
    --repeats)
      RPT="$2"
      shift
      shift
      ;;
    --bam)
      BAM="$2"
      shift
      shift
      ;;
    --out)
      OUT="$2"
      shift
      shift
      ;;
    *)
      echo "Unknown option: $1"
      exit 1
      ;;
  esac
done

MEM_GB=$((SLURM_MEM_PER_NODE / 1024))
echo "Mem GB: ${MEM_GB}"

mkdir -p ${OUT}
velocyto run \
  --samtools-threads ${SLURM_CPUS_PER_TASK} \
  --samtools-memory ${MEM_GB} \
  -b ${BARCODES} \
  -o ${OUT} \
  -m ${RPT} \
  ${BAM} \
  ${GENES}


echo "Finished"
