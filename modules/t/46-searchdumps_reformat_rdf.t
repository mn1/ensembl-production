#!perl
# Copyright [2009-2017] EMBL-European Bioinformatics Institute
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#      http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
use warnings;
use strict; 
use Test::More;
use File::Slurp;
use FindBin qw( $Bin );
use File::Spec;
use Bio::EnsEMBL::Test::MultiTestDB;

BEGIN {
	use_ok('Bio::EnsEMBL::Production::Search::RDFFormatter');
}

diag("Testing ensembl-production Bio::EnsEMBL::Production::Search, Perl $], $^X"
);
my $test     = Bio::EnsEMBL::Test::MultiTestDB->new('homo_sapiens_dump');
my $core_dba = $test->get_DBAdaptor('core');

my $test_onto = Bio::EnsEMBL::Test::MultiTestDB->new('multi');
my $onto_dba  = $test_onto->get_DBAdaptor('ontology');
my $genome_in_file = File::Spec->catfile( $Bin, "genome_test.json" );
my $genes_in_file = File::Spec->catfile( $Bin, "genes_test.json" );

my $formatter = Bio::EnsEMBL::Production::Search::RDFFormatter->new(
			-CONFIG_FILE => File::Spec->catfile( $Bin, "xref_LOD_mapping.json" ),
			-ONTOLOGY_DBA    => $onto_dba
);
subtest "RDF core", sub {
	my $out_file = File::Spec->catfile( $Bin, "rdf_core_test.ttl" );
	$formatter->reformat_core($genome_in_file, $genes_in_file, $core_dba, $out_file );
	ok(-e $out_file, "File exists");
	unlink $out_file;
};
subtest "RDF xrefs", sub {
	my $out_file = File::Spec->catfile( $Bin, "rdf_xrefs_test.ttl" );
	$formatter->reformat_xrefs($genome_in_file, $genes_in_file, $core_dba, $out_file );
	ok(-e $out_file, "File exists");
	unlink $out_file;
};

done_testing;
