#!/bin/sh -e

case $(uname) in
FreeBSD)
    pv=38
    pkg install -y python3 py$pv-matplotlib py$pv-dna-features-viewer
    ;;

*)
    cat << EOM

$(uname) is not yet supported.  Please consider helping out by adding
a case for $(uname) and submitting it to the project on Github.

EOM
    ;;

esac
