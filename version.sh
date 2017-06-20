#!/bin/sh

# Master version, for use in tarballs or non-git source copies
VERSION=1.4.1

# If we have a git clone, then check against the current tag
if [ -e .git ]
then
    # If we ever get to 10.x this will need to be more liberal
    VERSION=`git describe --match '[0-9].[0-9]*' --dirty`
fi

# Numeric version is for use in .dylib or .so libraries
#
# Follows the same logic from the Makefile commit c2e93911
# as non-numeric versions get bumped to patch level 255 to indicate
# an unknown value.
if [ "$1" = "numeric" ]
then
    v1=`expr "$VERSION" : '\([0-9]*\)'`
    v2=`expr "$VERSION" : '[0-9]*.\([0-9]*\)'`
    v3=`expr "$VERSION" : '[0-9]*.[0-9]*.\([0-9]*\)'`
    if [ -z "`expr "$VERSION" : '^\([0-9.]*\)$'`" ]
    then
	VERSION="$v1.$v2.255"
    else
	VERSION="$v1.$v2${v3:+.}$v3"
    fi
fi

echo $VERSION
