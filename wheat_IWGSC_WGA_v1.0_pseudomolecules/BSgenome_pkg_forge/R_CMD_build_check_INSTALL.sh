#!/bin/bash

# Build, check and install forged BSgenome package

# Usage:
# ./R_CMD_build_check_INSTALL.sh BSgenome.Taestivum.Cambridge.iwgsc_refseqv1_0 '1.50.0'

pkgName=$1
pkgVersion=$2

R CMD build ${pkgName} 
R CMD check ${pkgName}"_"${pkgVersion}".tar.gz"
R CMD INSTALL ${pkgName}"_"${pkgVersion}".tar.gz"
