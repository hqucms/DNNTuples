#!/usr/bin/env python
import os
import sys
import re
import argparse
import logging

def create(inputdir):
    train_val_files = []
    test_files = []

    for dp, dn, filenames in os.walk(inputdir):
        if 'failed' in dp or 'ignore' in dp:
            continue
        for f in filenames:
            if not f.endswith('.root'):
                continue
            fullpath = os.path.join(dp, f)
            try:
                filesize = os.path.getsize(fullpath)
                if filesize > 1000:
                    relpath = os.path.relpath(fullpath, start=inputdir)
                    if 'test_sample' in dp:
                        test_files.append(relpath)
                    else:
                        train_val_files.append(relpath)
                else:
                    logging.warning('Ignore file %s: size=%d' % (fullpath, filesize))
            except OSError:
                logging.warning('Ignore file %s: IO Error' % fullpath)

    with open(os.path.join(inputdir, 'train_val_samples.txt'), 'w') as f:
        f.write('\n'.join(train_val_files))
        print 'Write training/validation files to %s' % os.path.join(inputdir, 'train_val_samples.txt')

    with open(os.path.join(inputdir, 'test_samples.txt'), 'w') as f:
        f.write('\n'.join(test_files))
        print 'Write testing files to %s' % os.path.join(inputdir, 'test_samples.txt')

def main():
    parser = argparse.ArgumentParser(description='Create filelist for merging.\nFiles reserved for testing should be placed in *test_sample* directory!')
    parser.add_argument("inputdir", help="Input (parent) directory for a specific sample.")
    args = parser.parse_args()

    subdirs = os.listdir(args.inputdir)
    for s in subdirs:
        p = os.path.join(args.inputdir, s)
        if os.path.isdir(p):
            print p
            create(p)

if __name__ == '__main__':
    main()
