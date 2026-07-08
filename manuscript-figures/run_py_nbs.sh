#!/bin/bash

for dir in figure-* supplementary-*; do
  if [ -d "$dir" ]; then
    echo "=========================================="
    echo "Entering directory: $dir"
    echo "=========================================="
    for nb in "$dir"/*.ipynb; do
      if [ -f "$nb" ]; then
        kernel=$(python3 -c "import json; nb=json.load(open('$nb')); print(nb.get('metadata',{}).get('kernelspec',{}).get('language','unknown'))")
        if [ "$kernel" != "python" ]; then
          echo "[SKIP]  Non-Python notebook ($kernel), skipping: $nb"
          echo ""
          continue
        fi

        echo "[START] Running: $nb"
        jupyter nbconvert --to notebook --execute --inplace \
          --ExecutePreprocessor.timeout=-1 \
          --ExecutePreprocessor.kernel_name=ndmm-scrna \
          --ExecutePreprocessor.shutdown_kernel=immediate \
          "$nb"
        if [ $? -eq 0 ]; then
          echo "[DONE]  Successfully finished: $nb"
        else
          echo "[FAIL]  Error running: $nb"
        fi
        echo ""
      fi
    done
    echo "Finished all notebooks in: $dir"
    echo ""
  fi
done

echo "=========================================="
echo "All directories processed!"
echo "=========================================="