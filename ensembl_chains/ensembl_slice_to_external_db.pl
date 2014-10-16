#!/usr/bin/env perl
# Copyright [1999-2014] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
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

=pod

=head1 CHAIN FILES

This code writes chain files from Ensembl assembly records. Chain files are defined at

  http://genome.ucsc.edu/goldenPath/help/chain.html

They can be used with liftover (http://genome.ucsc.edu/cgi-bin/hgLiftOver) and crossmap (http://crossmap.sourceforge.net). The file format says the following:

  chain score tName tSize tStrand tStart tEnd qName qSize qStrand qStart qEnd id 
  size dt dq
  size

Majority of fields are self explanatory expect for the declaration of target & query sequences. In a mapping of NCBI36 (hg18) to GRCh37 (hg19) the following roles are designated:

  - target == NCBI36
  - query  == GRCh37

You can see this in the the example chain files for human mappings http://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg18.over.chain.gz and http://hgdownload.soe.ucsc.edu/goldenPath/hg18/liftOver/hg18ToHg19.over.chain.gz.

dt & dq are meant to be the differences between the end of this current block and the start of the next in the target and query. The final line is just a size to indicate no further offset.

tSize and qSize are the total sizes of the target and reference sequences.

=cut

use strict;
use warnings;
use Getopt::Long;
use Bio::EnsEMBL::Registry;
use Bio::EnsEMBL::Utils::IO qw/work_with_file gz_work_with_file/;
use feature qw/say/;
use File::Path qw/mkpath/;
use File::Spec;

my $chain_id = 1;
my %global_mapings;

sub get_options {
  my ($db_name, $db_host, $db_user, $db_pass, $db_port, $help, @species, $group, $release, $dir, $compress, @external_db);
  # my ($asm, $cmp) = ('GRCh37', 'NCBI36');
  $db_port = 3306;
  # $species = 'human';
  $group = 'core';
  $dir = File::Spec->curdir();

  GetOptions(
    "db_name|dbname|database=s"       => \$db_name,
    "db_host|dbhost|host=s"           => \$db_host,
    "db_user|dbuser|user|username=s"  => \$db_user,
    "db_pass|dbpass|pass|password=s"  => \$db_pass,
    "db_port|dbport|port=s"           => \$db_port,
    "species=s@"                      => \@species,
    "dir=s"                           => \$dir,
    "compress!"                       => \$compress,
    "version|release=i"               => \$release,
    "external_db=s@"                  => \@external_db,
    "h!"                              => \$help,
    "help!"                           => \$help,
  );

  die "No -external_db given" unless @external_db;
  
  Bio::EnsEMBL::Registry->load_registry_from_db(
    -HOST => $db_host, -PORT => $db_port, 
    -USER => $db_user, -PASS => $db_pass,
    -DB_VERSION => $release
  );

  my @dbas;
    
  if(@species) {
    say "Working against a restricted species list";
    foreach my $s (@species) {
      my $dba = Bio::EnsEMBL::Registry->get_DBAdaptor($s, $group);
      die "Cannot find a DBAdaptor for the species ${s}" unless $dba;
      push(@dbas, $dba);
    }
  }
  else {
    say "Dumping chain file for all available species";
    @dbas = sort {$a->dbc->dbname cmp $b->dbc->dbname} @{Bio::EnsEMBL::Registry->get_all_DBAdaptors(-GROUP => 'core')};
  }

  return {
    dbas => \@dbas,
    external_dbs => \@external_db,
    dir => $dir,
    compress => $compress
  };
}

run();

sub run {
  my $opts = get_options();
  foreach my $dba (@{$opts->{dbas}}) {
    run_for_dba($dba, $opts);
    $dba->dbc->disconnect_if_idle;
  }
  return;
}

sub run_for_dba {
  my ($core_dba, $opts) = @_;
  my $prod_name = $core_dba->get_MetaContainer->get_production_name();

  say "Working with ${prod_name}";
  say "\tFetching slices";
  my $slices = slices($core_dba);
  my $data;
  foreach my $external_db (@{$opts->{external_dbs}}) {
    say "\tProcessing ${external_db}";

    if(! should_run($core_dba, $external_db)) {
      say "\tNo synonyms found. Nothing to convert";
      next;
    }
    $data = 1;
    write_mappings($opts, $core_dba, $slices, $external_db, 1);
    write_mappings($opts, $core_dba, $slices, $external_db, 0);
    say "\tDone";
    say q{};
  }
  write_tab_lookup($opts, $core_dba, $prod_name) if $data;
}

sub write_tab_lookup {
  my ($opts, $core_dba, $prod_name) = @_;

  my @external_dbs = @{$opts->{external_dbs}};
  my $assembly = $core_dba->get_GenomeContainer->get_version();
  my $file = "${assembly}_region_lookup.tsv";
  my $dir = File::Spec->catdir($opts->{dir}, $prod_name);
  mkpath($dir);
  my $path = File::Spec->catfile($dir, $file);

  say "\tWriting tab lookup to $path";
  write_to_file($path, $opts, sub {
    my ($fh) = @_;
    say $fh '#'.join("\t", ('Ensembl', @external_dbs));
    foreach my $slice (@{slices($core_dba)}) {
      my $seq_region_name = $slice->seq_region_name();
      next unless exists $global_mapings{$prod_name}{$seq_region_name};
      say $fh join("\t", $seq_region_name, map { defined($_) ? $_ : '.' } map {$global_mapings{$prod_name}{$seq_region_name}{$_}} @external_dbs);
    }
  });
  return;
}

sub should_run {
  my ($core_dba, $external_db) = @_;
  my $sql = <<'SQL';
select count(*) 
from seq_region 
join seq_region_synonym using(seq_region_id) 
join external_db using (external_db_id) 
join coord_system using (coord_system_id)
where db_name = ?
and species_id =?
SQL
  my $params = [$external_db, $core_dba->species_id()];
  return $core_dba->dbc->sql_helper->execute_single_result(-SQL => $sql, -PARAMS => $params);
}

sub write_mappings {
  my ($opts, $core_dba, $slices, $external_db, $to_ensembl) = @_;
  
  my $assembly = $core_dba->get_GenomeContainer->get_version();
  my $prod_name = $core_dba->get_MetaContainer->get_production_name();
  my $direction = ($to_ensembl) ? "from ${external_db} to Ensembl" : "from Ensembl to ${external_db}";
  say "\tBuilding chain mappings $direction";
  my $chains = build_chains($slices, $external_db, $to_ensembl, $prod_name);

  # Setup path
  my $file = ($to_ensembl) ? "${external_db}ToEnsembl.chain" : "EnsemblTo${external_db}.chain";
  my $dir = File::Spec->catdir($opts->{dir}, $prod_name, $assembly);
  mkpath($dir);
  my $path = File::Spec->catfile($dir, $file);

  say "\tWriting mappings to $path";
  write_to_file($path, $opts, sub {
    my ($fh) = @_;
    print_chains($fh, $chains);
  });
  return;
}

sub slices {
  my ($core_dba) = @_;
  my %hash; 
  my $slices = $core_dba->get_SliceAdaptor->fetch_all('toplevel', undef, 'nonref');
  while(my $slice = shift @{$slices}) {
    $hash{$slice->seq_region_name()} = $slice;
  }
  return [
    map { $_->[0] }
    sort { $a->[1] <=> $b->[1] || $a->[2] cmp $b->[2] } 
    map { my $rank = $_->karyotype_rank(); $rank = 1e8 if $rank == 0; [$_, $rank, $_->seq_region_name()] } 
    values %hash
  ];
}

sub build_chains {
  my ($slices, $external_db, $to_ensembl, $prod_name) = @_;
  my @chains;
  foreach my $slice (@{$slices}) {
    my $chain = build_chain($slice, $external_db, $to_ensembl, $prod_name);
    next unless defined $chain;
    push(@chains, $chain);
  }
  return \@chains;
}

# trg is the source data & q is the target i.e. converting from UCSC to Ensembl means trg == chr1 & q == 1
sub build_chain {
  my ($slice, $external_db, $to_ensembl, $prod_name) = @_;
  my ($synonym) = @{$slice->get_all_synonyms($external_db)};
  my $seq_region_name = $slice->seq_region_name();
  return unless $synonym;
  #Write into a global hash so we know what the synonyms were for the TSV file
  $global_mapings{$prod_name}{$seq_region_name}{$external_db} = $synonym->name();
  my ($t_name, $q_name) = ($to_ensembl) ? ($synonym->name(), $seq_region_name) : ($seq_region_name, $synonym->name());
  my $size = $slice->seq_region_length();
  return {
    # header => ['chain', $chain_score, $t_name, $t_size, $t_strand, $t_start, $t_end, $q_name, $q_size, $q_strand, $q_start, $q_end, $chain_id],
    header => ['chain', 1, $t_name, $size, '+', 0, $size, $q_name, $size, '+', 0, $size, $chain_id++],
    # gaps => [[match t_gap_to_next_block q_gap_to_next_block]]
    gaps => [[$size]]
  };
}

sub print_chains {
  my ($fh, $chains) = @_;
  while(my $chain_def = shift @{$chains}) {
    my $header = $chain_def->{header};
    my $gaps = $chain_def->{gaps};
    print $fh join(q{ }, @{$header});
    print $fh "\n";
    foreach my $gap (@{$gaps}) {
      print $fh join(q{ }, @{$gap});
      print $fh "\n";
    }
    print $fh "\n";
  }
}

sub write_to_file {
  my ($path, $opts, $writer) = @_;
  if($opts->{compress}) {
    $path .= '.gz';
    gz_work_with_file($path, 'w', $writer);
  }
  else {
    work_with_file($path, 'w', $writer);
  }
  return;
}
