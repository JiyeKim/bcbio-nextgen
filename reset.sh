#!/bin/sh

echo "remove old packages..."
rm -rf /incogwas/igsite/storage/bcbio_dev/env/lib/python2.7/site-packages/bcbio*

echo "setup..."
python setup.py install

echo "testing..."
./tests/run_tests.sh dgseq
