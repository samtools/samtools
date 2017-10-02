#!/bin/sh

# Master version, for use in tarballs or non-git source copies
VERSION=1.6

# If we have a git clone, then check against the current tag
if [ -e .git ]
then
    # If we ever get to 10.x this will need to be more liberal
    VERSION=`git describe --match '[0-9].[0-9]*' --dirty`
fi

echo $VERSION
