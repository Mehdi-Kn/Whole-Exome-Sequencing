#!/usr/bin/env bash
###############################################################################
# parse_args.sh
# Parse command-line arguments (--fastq, --out, --threads).
# Usage: source parse_args.sh "$@"
###############################################################################

usage() {
  echo "Usage: $0 [--fastq <FASTQ_DIR>] [--out <OUTPUT_DIR>] [--threads <N>]"
  exit 1
}

# Initialize empty variables (in case they are not specified)
FASTQ_DIR=""
OUTPUT_DIR=""

while [[ $# -gt 0 ]]; do
  key="$1"
  case $key in
    --fastq)
      FASTQ_DIR="$2"
      shift; shift
      ;;
    --out)
      OUTPUT_DIR="$2"
      shift; shift
      ;;
    --threads)
      THREADS="$2"
      shift; shift
      ;;
    *)
      usage
      ;;
  esac
done

# Check for mandatory arguments
if [[ -z "${FASTQ_DIR}" || -z "${OUTPUT_DIR}" ]]; then
  usage
fi
