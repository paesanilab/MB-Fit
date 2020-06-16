#!/bin/bash

if [ -z "$MBFIT_HOME" ]; then
    echo "Error: MBFIT_HOME not set"
else
    if [ -z "$PYTHONPATH" ]; then
        export PYTHONPATH="${MBFIT_HOME}"
    else
	export PYTHONPATH="$PYTHONPATH:${MBFIT_HOME}"
    fi
fi
