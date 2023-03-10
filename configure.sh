#!/bin/bash

# Create the cache folder if it does not exist
if [ ! -d "$HOME/.gasscli" ]; then
    mkdir "$HOME/.gasscli"
    mkdir "$HOME/.gasscli/cache"
    cp config.json "$HOME/.gasscli/config.json"
fi
