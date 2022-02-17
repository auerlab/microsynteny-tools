#!/bin/sh -e

: ${PREFIX:=../local}
: ${LOCALBASE:=/usr/local}

for file in *.c; do
    cproto -I$PREFIX/include -I$LOCALBASE/include $file
done
