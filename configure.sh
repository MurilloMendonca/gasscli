#!/bin/bash

# Create the cache folder if it does not exist
if [ ! -d "$HOME/.gasscli" ]; then
    mkdir "$HOME/.gasscli"
    mkdir "$HOME/.gasscli/cache"
    cd "confs"
    for file in *; do
        cp "$file" "$HOME/.gasscli/$file"
        sed -i "s|/home/username|$HOME|g" "$HOME/.gasscli/$file" # Replace the username with the current user's home directory
    done
fi
    
