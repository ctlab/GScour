#!/usr/bin/env python
# -*- coding: utf-8 -*-
import multiprocessing
import os
import logging

logging.basicConfig(format='%(asctime)s : %(levelname)s : %(message)s', level=logging.INFO)

directory_in = '/home/alina_grf/progprojects/data/outprank2/'
directory_out = '/home/alina_grf/progprojects/data/outgblocks/'


def parse_dir():
    for infile in os.listdir(directory_in):
        yield infile


def get_data_to_gblocks(infile):
    #print(infile)
    os.chdir(directory_out)
    launch = '/home/alina_grf/biotools/Gblocks_0.91b/Gblocks {0} -b1=3 -b2=3' \
             ' -b3=7 -b4=3 -b5=h -b6=Yes -p=Yes'.format(os.path.join(directory_in, infile))
    os.system(launch)


if __name__ == '__main__':
    multiprocessing.log_to_stderr()
    logger = multiprocessing.get_logger()
    logger.setLevel(logging.INFO)
#with Pool(2) as p:
    pool = multiprocessing.Pool(3)
    inputs = list(parse_dir())
    pool.map(get_data_to_gblocks, inputs)
