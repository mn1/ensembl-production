
=head1 LICENSE

Copyright [2009-2016] EMBL-European Bioinformatics Institute

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

     http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.

=cut

package Bio::EnsEMBL::Production::Search::RDFFormatter;

use warnings;
use strict;
use Bio::EnsEMBL::Utils::Exception qw(throw);
use Bio::EnsEMBL::Utils::Argument qw(rearrange);
use Data::Dumper;
use Log::Log4perl qw/get_logger/;
use Bio::EnsEMBL::IO::Object::RDF;
use Bio::EnsEMBL::IO::Translator::Slice;
use Bio::EnsEMBL::IO::Translator::BulkFetcherFeature;
use Bio::EnsEMBL::IO::Writer::RDF;
use Bio::EnsEMBL::IO::Writer::RDF::XRefs;

use JSON;
use Carp;
use File::Slurp;

use Exporter 'import';
our @EXPORT = qw(_array_nonempty _id_ver _base);

sub new {
	my ( $class, @args ) = @_;
	my $self = bless( {}, ref($class) || $class );
	( $self->{config_file}, $self->{onto_dba} ) =
	  rearrange( [qw/config_file ontology_dba/], @args );
	  
	$self->{log} = get_logger();
	return $self;
}

sub log {
	my ($self) = @_;
	return $self->{log};
}

sub reformat_core {
	my ( $self, $genome_file, $genes_file, $dba, $outfile, $slices ) = @_;

	if ( !defined $slices ) {
		$slices = [ sort { $b->length() <=> $a->length() }
					  @{$dba->get_SliceAdaptor()
						  ->fetch_all( 'toplevel', undef, 1, undef, undef ) }
		];
	}

	my $genome = decode_json( read_file($genome_file) );
	my $genes  = decode_json( read_file($genes_file) );

	### Dump core RDF ###

	my $meta_adaptor = $dba->get_MetaContainer();

	$genome->{dbname} =~ m/.*_([a-z]+)_([0-9]+)_([0-9]+)(_([0-9]+))?/;
	my $type    = $1;
	my $release = $2;

	# start writing out: namespaces and species info
	my $fh = IO::File->new( $outfile, "w" ) || die "$! $outfile";
	my $core_writer = Bio::EnsEMBL::IO::Writer::RDF->new();
	$core_writer->open($fh);
	$core_writer->write( Bio::EnsEMBL::IO::Object::RDF->namespaces() );
	$core_writer->write(
					Bio::EnsEMBL::IO::Object::RDF->species(
						taxon_id        => $genome->{organism}{taxonomy_id},
						scientific_name => $genome->{organism}{scientific_name},
						common_name     => $genome->{organism}{display_name} )
	);

	# write sequence regions
	my $slice_trans =
	  Bio::EnsEMBL::IO::Translator::Slice->new( version      => $release,
												meta_adaptor => $meta_adaptor );
	map { $core_writer->write( $_, $slice_trans ) } @{$slices};

	# write BulkFetcher 'features'
	my $feature_trans = Bio::EnsEMBL::IO::Translator::BulkFetcherFeature->new(
		version           => $release,
		xref_mapping_file => $self->{config_file},
		,    # required for mapping Ensembl things to RDF
		ontology_adaptor => $self->{onto_dba}->get_OntologyTermAdaptor(),
		meta_adaptor     => $meta_adaptor );
	map { $core_writer->write( $_, $feature_trans ) } @{$genes};

	# finally write connecting triple to master RDF file
	$core_writer->write(
					  Bio::EnsEMBL::IO::Object::RDF->dataset(
						  version => $release,
						  project => $meta_adaptor->get_division() || 'Ensembl',
						  production_name => $genome->{organism}{name} ) );
	$core_writer->close();
	return;
} ## end sub reformat_core

sub reformat_xrefs {
	my ( $self, $genome_file, $genes_file, $dba, $outfile ) = @_;

	my $genome = decode_json( read_file($genome_file) );
	my $genes  = decode_json( read_file($genes_file) );

	$genome->{dbname} =~ m/.*_([a-z]+)_([0-9]+)_([0-9]+)(_([0-9]+))?/;
	my $type    = $1;
	my $release = $2;
	### Xrefs RDF ###
	#
	my $fh = IO::File->new( $outfile, "w" ) || die "$! $outfile";
	my $feature_trans = Bio::EnsEMBL::IO::Translator::BulkFetcherFeature->new(
		version => $release,
		xref_mapping_file =>
		  $self->{config_file},     # required for mapping Ensembl things to RDF
		ontology_adaptor => $self->{onto_dba}->get_OntologyTermAdaptor(),
		meta_adaptor     => $dba->get_MetaContainer() );

	my $xrefs_writer =
	  Bio::EnsEMBL::IO::Writer::RDF::XRefs->new($feature_trans);
	$xrefs_writer->open($fh);

	# write namespaces
	$xrefs_writer->write( Bio::EnsEMBL::IO::Object::RDF->namespaces() );

	# then dump feature xrefs
	map { $xrefs_writer->write($_) } @{$genes};

	$xrefs_writer->close();
	return;
} ## end sub reformat_xrefs

1;
