# Build

```bash
$ make
```

# Dependencies

* [boost](http://www.boost.org/)
* [libcurl](http://curl.haxx.se/libcurl/)

# Usage

## Download proteins to a cache folder
```bash
$ ./gasscli download -c <cache folder> -p <protein list> -f <protein list file name>
```
Example:
```bash
$ ./gasscli download -c cache/ -p 3NOS 5CNV P29474
```

## Run GASS

```bash
$ ./gasscli run -c <cache folder> -t <target protein> -s <template site> -r <reference protein> --mutations <mutation> 
```
Where:

* `<cache folder>` is the folder where the proteins are stored
* `<target protein>` is the target protein code in pdb or uniprot format, if not on cache folder it will be downloaded
* `<template site>` is the template site in the format 'RESIDUE,NUMBER,CHAIN;' Ex: 'GLU,100,A;'
* `<reference protein>` is the reference protein code in pdb or uniprot format, if not on cache folder it will be downloaded
* `<mutation>` is the mutation in the format 'RESIDUE,RESIDUE;' Ex: 'GLU,ALA;'

Example:
```bash
$ ./gasscli run -r 3NOS -s "CYS,99,B;CYS,99,A;CYS,94,A;CYS,94,B;" -t P29474
```