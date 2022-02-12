#!/bin/sh -e

pyversion=38

pkg install -y python3 py$pyversion-matplotlib
# Replace this with pkg install once the port is well-tested
wip-reinstall-port -u -r py-dna-features-viewer
