#!/bin/bash

# Check if a job name was provided
if [ -z "$1" ]; then
  echo "Error: No job name provided."
  echo "Usage: ./run_abaqus.sh <job_name>"
  exit 1
fi

JOBNAME=$1

echo "Starting Abaqus job: $JOBNAME with uel.f90..."

# Run the abaqus command with the provided job name
abaqus job=$JOBNAME user=uel.f90 ask_delete=no interactive > ${JOBNAME}_output.txt

echo "Job $JOBNAME finished! Output saved to ${JOBNAME}_output.txt"
