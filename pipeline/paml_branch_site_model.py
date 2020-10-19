#!/usr/bin/env python
import argparse
import sys
from Bio.Phylo.PAML import codeml
import logging
import os
import logging

logging.basicConfig(format='%(asctime)s : %(levelname)s : %(message)s', level=logging.INFO)
WRITTEN_FILES = 0
EXCEPTION_NUMBER = 0
FILES_NUMBER = 0
"""There are two hypothesis:
H0: model = 2, NSsites = 2 (branch-site model),
    fix_kappa = 0   * 1: kappa fixed, 0: kappa to be estimated
    kappa = 2   * initial or fixed kappa
    fix_omega = 1   * 1: omega or omega_1 fixed, 0: estimate
    omega = 1   * initial or fixed omega, for codons or codon-based AAs
H1: model = 2, NSsites = 2 (branch-site model),
    fix_kappa = 0   * 1: kappa fixed, 0: kappa to be estimated
    kappa = 2   * initial or fixed kappa
    fix_omega = 0   * 1: omega or omega_1 fixed, 0: estimate
    omega = 1   * initial or fixed omega, for codons or codon-based AAs
    
in analyse_pamlout.py:
For analysis performs: ln0, np0 from H0; ln1, np1 from H1; 
                       ΔLRT = 2×(lnL1 - lnL0)
                       n = np1 - np0
                       p_val = 1 - stats.chi2.cdf(ΔLRT, n)

"""

"""
BROKEN_LIST = ['1002', '1006', '1010', '1014', '1037', '1041', '1044', '1051', '1053', '1059', '106', '1077', '107',
               '1089', '1099', '10', '1102', '1106', '1107', '1111', '1113', '1115', '1116', '1117', '1120', '1121', 
               '1124', '1125', '1127', '1139', '1152', '1154', '1158', '1159', '116', '1178', '1190', '11', '1203', 
               '1213', '1234', '1239', '1252', '1254', '1255', '1257', '1265', '1273', '127', '1287', '1293', '1299',
                '129', '12', '1306', '1344', '1373', '1384', '138', '1392', '1394', '1399', '1409', '1413', '1433', 
                '1437', '1441', '1444', '1458', '1466', '1467', '1477', '1482', '1483', '1490', '1493', '1500', '1501',
                 '1506', '1511', '1519', '1524', '1527', '1528', '1530', '1535', '1538', '1539', '1546', '1547', '154',
                  '1559', '1571', '1573', '1578', '1579', '1582', '1584', '1596', '1598', '1602', '1603', '1605', '1609',
                   '1618', '1630', '1654', '1658', '1661', '1662', '1668', '1669', '1672', '1687', '1699', '1706', '1709',
                    '1710', '1714', '1718', '171', '1725', '1727', '1728', '1736', '1738', '1739', '1741', '1742', '1759',
                     '1760', '1762', '1764', '1772', '1773', '177', '1782', '1792', '1793', '1797', '1798', '1805', '1807',
                      '180', '1811', '1814', '181', '1833', '1835', '1837', '1840', '1841', '1842', '1846', '1853', '1854', 
                      '1858', '1860', '1865', '1868', '186', '1871', '1875', '1887', '1888', '1892', '1893', '1895', '189',
                       '18', '1904', '1912', '1923', '1948', '1962', '1965', '1974', '1977', '1978', '1979', '1982', '1989',
                        '1990', '1997', '19', '1', '2004', '2021', '2034', '2054', '2055', '2061', '2064', '2065', '206',
                         '2076', '2079', '2089', '208', '2091', '2093', '2097', '20', '2101', '2103', '2105', '2109', '2110', '2113', '2115', '2116', '2120', '2121', '2126', '2135', '2142', '2146', '2152', '2154', '2160', '2168', '2169', '216', '2173', '2174', '2175', '217', '2181', '2189', '2196', '219', '21', '2200', '2202', '2204', '2212', '2223', '2227', '2229', '2231', '2249', '2263', '2276', '2278', '2291', '2292', '2294', '22', '2304', '2308', '2309', '2314', '2315', '2316', '2318', '231', '2327', '2328', '2331', '2336', '233', '2346', '2348', '2360', '2372', '2374', '2378', '2379', '2385', '2397', '2399', '2401', '2404', '2429', '2438', '2443', '2445', '2450', '2452', '2461', '2464', '2478', '2481', '2482', '2487', '2491', '2498', '24', '2505', '2506', '2512', '2519', '2521', '2525', '2542', '2543', '2551', '2553', '2558', '2564', '2586', '2594', '2596', '2603', '2614', '2616', '2617', '2628', '2634', '2638', '2639', '2642', '2649', '2657', '2658', '2661', '2670', '2679', '2680', '2681', '2696', '2700', '2704', '270', '2713', '2716', '2759', '2761', '2766', '2767', '2769', '2772', '2780', '2785', '2787', '278', '2791', '2798', '27', '2816', '2835', '2841', '2843', '2862', '2864', '2865', '2871', '2873', '2876', '2882', '2886', '2888', '2895', '2897', '2905', '2924', '2943', '2947', '2957', '2958', '2970', '2972', '2987', '29', '2', '3006', '3017', '3028', '3030', '3031', '3035', '3039', '3040', '3043', '305', '3065', '3070', '3075', '3079', '307', '30', '3111', '3119', '3123', '3130', '3139', '3141', '3152', '3157', '3171', '3187', '3190', '3215', '3218', '3226', '323', '3243', '3247', '3264', '3287', '3297', '32', '3318', '331', '3322', '3323', '3326', '3328', '3341', '3342', '3353', '3358', '3362', '3368', '3375', '3378', '3386', '3390', '3391', '3410', '3414', '3419', '3424', '3433', '3435', '3436', '3442', '3446', '3448', '3455', '3459', '3463', '3474', '3489', '3493', '3495', '3500', '3506', '3508', '3516', '3523', '3525', '3527', '3532', '3533', '3537', '3545', '3546', '3549', '3550', '3553', '3555', '3556', '3560', '35', '360', '366', '36', '370', '371', '373', '37', '389', '39', '3', '408', '40', '420', '432', '436', '43', '46', '480', '484', '488', '489', '496', '498', '499', '49', '4', '501', '504', '50', '510', '51', '535', '537', '53', '542', '559', '561', '565', '566', '567', '56', '576', '59', '605', '61', '627', '62', '630', '640', '644', '647', '66', '672', '675', '680', '68', '690', '6', '700', '710', '718', '725', '72', '732', '734', '741', '743', '752', '753', '75', '762', '763', '76', '771', '772', '773', '775', '779', '78', '790', '791', '7', '815', '816', '81', '821', '830', '831', '837', '84', '858', '85', '872', '874', '875', '877', '885', '88', '892', '898', '903', '91', '926', '931', '935', '947', '956', '95', '96', '986', '987', '98', '99', '9', '1173', '1245', '1315', '1396', '160', '1644', '1924', '2122', '2321', '2604', '2808', '3076', '343', '3444', '616', '638', '891']
"""


def parse_dir(infolder):
    global FILES_NUMBER
    for personal_folder in os.scandir(infolder):
        if os.path.isdir(personal_folder):
            for infile in os.listdir(personal_folder):
                if infile.split('.')[-1] == 'phy':
                    FILES_NUMBER += 1
                    yield os.path.join(infolder, personal_folder, infile)


def set_alternative_hypothesis(infile, tree, personal_dir):
    file_out_path = infile.replace('.phy', '_alter1.out')
    cml = codeml.Codeml(
        alignment=infile,
        tree=tree,
        out_file=file_out_path,
        working_dir=personal_dir,
        )
    cml.set_options(noisy=9)
    cml.set_options(verbose=0)
    cml.set_options(runmode=0)
    cml.set_options(seqtype=1)
    cml.set_options(CodonFreq=2)
    cml.set_options(clock=0)
    cml.set_options(aaDist=0)
    cml.set_options(model=2)
    cml.set_options(NSsites=[2])
    cml.set_options(icode=0)
    cml.set_options(Mgene=0)
    cml.set_options(fix_kappa=0)
    cml.set_options(kappa=2)
    cml.set_options(fix_omega=0)
    cml.set_options(omega=1)
    cml.set_options(getSE=0)
    cml.set_options(RateAncestor=0)
    cml.set_options(Small_Diff=.45e-6)
    cml.set_options(cleandata=1)
    cml.set_options(fix_blength=1)
    return cml, file_out_path


def set_null_hypothesis(infile, tree, personal_dir):
    file_out_path = infile.replace('.phy', '_null1.out')
    cml = codeml.Codeml(
        alignment=infile,
        tree=tree,
        out_file=file_out_path,
        working_dir=personal_dir,
        )
    cml.set_options(noisy=9)
    cml.set_options(verbose=0)
    cml.set_options(runmode=0)
    cml.set_options(seqtype=1)
    cml.set_options(CodonFreq=2)
    cml.set_options(clock=0)
    cml.set_options(aaDist=0)
    cml.set_options(model=2)
    cml.set_options(NSsites=[2])
    cml.set_options(icode=0)
    cml.set_options(Mgene=0)
    cml.set_options(fix_kappa=0)
    cml.set_options(kappa=2)
    cml.set_options(fix_omega=1)
    cml.set_options(omega=1)
    cml.set_options(getSE=0)
    cml.set_options(RateAncestor=0)
    cml.set_options(Small_Diff=.45e-6)
    cml.set_options(cleandata=1)
    cml.set_options(fix_blength=1)
    return cml, file_out_path


def run_paml(infile, tree):
    global WRITTEN_FILES, EXCEPTION_NUMBER
    personal_dir = os.path.split(infile)[0]
    cml, file_out_path = set_null_hypothesis(infile, tree, personal_dir)
    try:
        cml.run(command="codeml", verbose=True) #/home/alina_grf/BIOTOOLS/paml4.9j/bin/
        logging.info("paml out file {} has been written".format(file_out_path))
        WRITTEN_FILES += 1
    except:
        logging.exception("null, infile {}, sys.exc_info() {}".format(infile, sys.exc_info()))
        EXCEPTION_NUMBER += 1

    cml, file_out_path = set_alternative_hypothesis(infile, tree, personal_dir)
    try:
        cml.run(command="codeml", verbose=True) #/home/alina_grf/BIOTOOLS/paml4.9j/bin/
        logging.info("paml out file {} has been written".format(file_out_path))
        WRITTEN_FILES += 1
    except:
        logging.exception("alter, infile {},  sys.exc_info() {}".format(infile, sys.exc_info()))
        EXCEPTION_NUMBER += 1


def main(infolder, tree):
    for infile in parse_dir(infolder):
        run_paml(infile, tree)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--infolder', help='The full path to the folder with input files for paml', nargs='?')
    parser.add_argument('--tree', help='Path to the tree for paml', nargs='?')
    args = parser.parse_args()
    try:
        main(args.infolder, args.tree)
        logging.info("Number of files have been analyzed: {}".format(FILES_NUMBER))
        logging.info("Number of written files: {}".format(WRITTEN_FILES))
        logging.info("Number of exceptions: {}".format(EXCEPTION_NUMBER))
    except:
        logging.exception("Unexpected error")

    logging.info("The work has been completed")
