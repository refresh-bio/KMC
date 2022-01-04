#!/usr/bin/env python3
import sys
import py_kmc_api as pka
import argparse
import textwrap

VER = "3.2.1"
DATE = "2022-01-04"

def GeneralHelp():
    print('KMC dump ver. {} ({})'.format(VER, DATE))

class VersionHelpAction(argparse.Action):
    def __call__(self, parser, namespace, values, option_string=None):        
        GeneralHelp()
        parser.print_help()        
        parser.exit()

class MyParamsParser(argparse.ArgumentParser):
    def error(self, message):
        sys.stderr.write('error: %s\n' % message)
        GeneralHelp()
        self.print_help()
        sys.exit(1)

parser = MyParamsParser(add_help=False)
parser.register('action', 'version_help_action', VersionHelpAction)


group1 = parser.add_argument_group('help and version')
group1.add_argument("--version", action='version_help_action', nargs=0, help="print help and version")
group1.add_argument("-h", "--help", action='version_help_action', nargs=0, help="print help and version")

group2 = parser.add_argument_group('normal run')
group2.add_argument("kmc_database", help="kmc database")
group2.add_argument("output_file", help="output file")
group2.add_argument("-ci", "--cutoff_min", type=int, help="exclude k-mers occurring less than CI times", default=0)
group2.add_argument("-cx", "--cutoff_max", type=int, help="exclude k-mers occurring more of than CX times", default=0)

args = parser.parse_args()

kmer_data_base = pka.KMCFile()
if not kmer_data_base.OpenForListing(args.kmc_database):
    print("Error: cannot open kmc database")
    sys.exit(1)

info = kmer_data_base.Info()
kmer_object = pka.KmerAPI(info.kmer_length)

if args.cutoff_min > 0:
    if not kmer_data_base.SetMinCount(args.cutoff_min):
        print("Error: cannot set cutoff min")
        sys.exit(1)

if args.cutoff_max > 0:
    if not kmer_data_base.SetMaxCount(args.cutoff_max):
        print("Error: cannot set cutoff max")
        sys.exit(1)

output_file = open(args.output_file, 'w')

counter = pka.Count()
while kmer_data_base.ReadNextKmer(kmer_object, counter):
    output_file.write("{}\t{}\n".format(kmer_object, counter.value))
output_file.close()