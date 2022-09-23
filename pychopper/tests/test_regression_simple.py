import unittest
import os
from os import path
import subprocess


class TestIntegration(unittest.TestCase):

    def testIntegration(self):
        """ Integration test. """
        base = path.dirname(__file__)
        test_base = path.join(base, 'data')

        barcodes = path.join(test_base, 'barcodes.fas')
        input_fasta = path.join(test_base, 'ref.fq')
        output_fasta = path.join(test_base, 'test_output.fq')
        expected_output = path.join(test_base, 'expected_output.fas')

        subprocess.call("{} {} {} {} {}".format('pychopper', "-Y 0 -B 3 -q 0.5 -m edlib -b", barcodes, input_fasta, output_fasta), shell=True)
        retval = subprocess.call(['cmp', output_fasta, expected_output])
        self.assertEqual(retval, 0)
        os.remove(output_fasta)

    def testIntegration_umi(self):
        """ Integration test. """
        base = path.dirname(__file__)
        test_base = path.join(base, 'data')

        input_fasta = path.join(test_base, 'PCS111_umi_test_reads.fastq')
        output_fasta = path.join(test_base, 'test_output_umi.fq')
        expected_output = path.join(test_base, 'PCS111_umi_test_reads_expected.fastq')

        subprocess.call("{} {} {} {}".format('pychopper', "-U -m edlib -k PCS111", input_fasta, output_fasta), shell=True)
        retval = subprocess.call(['cmp', output_fasta, expected_output])
        self.assertEqual(retval, 0)
        os.remove(output_fasta)

