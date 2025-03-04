#!/bin/bash

find . -type f \( -name "*.out" -o -name "*.err" \) -exec rm -v {} +

echo "Deletion completed."
