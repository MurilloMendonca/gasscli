#!/bin/bash

# Create the cache folder if it does not exist
if [ ! -d "$HOME/.gasscli" ]; then
    mkdir "$HOME/.gasscli"
fi
