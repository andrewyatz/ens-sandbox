#!/usr/bin/env perl
# Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
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

use strict;
use warnings;

sub basic_TO_JSON {
  my ($obj) = @_;
  my %shallow_clone = %{$obj};
  delete $shallow_clone{adaptor} if exists $shallow_clone{adaptor};
  return \%shallow_clone;
}

sub booleanise {
  my ($obj, $key) = @_;
  if(exists $obj->{$key}) {
    $obj->{$key} = ($obj->{$key}) ? $Types::Serialiser::true : $Types::Serialiser::false;
  }
  return;
}

sub numerify {
  my ($obj, $key) = @_;
  $obj->{$key} = $obj->{$key}+0 if exists $obj->{$key};
  return;
}

sub Bio::EnsEMBL::Storable::TO_JSON {
  my $clone = basic_TO_JSON(@_);
  numerify($clone, $_) for qw/dbID start end strand/;
  booleanise($clone, $_) for qw/is_current is_constitutive is_canonical/;
  return $clone;
}

sub Bio::EnsEMBL::Slice::TO_JSON {
  my ($slice) = @_;
  return {
    name => $slice->seq_region_name(),
    coord_system => $slice->coord_system(),
  };
}

sub Bio::EnsEMBL::CoordSystem::TO_JSON {
  my ($coord_system) = @_;
  return {
    name => $coord_system->name(),
    version => $coord_system->version()
  };
}

sub Bio::EnsEMBL::Attribute::TO_JSON {
  my ($attrib) = @_;
  return {
    name => $attrib->name(),
    code => $attrib->code(),
    value => $attrib->value()
  };
}


sub Bio::EnsEMBL::Analysis::TO_JSON {
  my ($self) = @_;
  return {
    logic_name => $self->logic_name(),
  };
}

sub Bio::EnsEMBL::Transcript::TO_JSON {
  my ($self) = @_;
  my $shallow = Bio::EnsEMBL::Storable::TO_JSON($self);
  $shallow->{cdna} = $self->spliced_seq();
  $shallow->{cds} = $self->translateable_seq();
  return $shallow;
}

sub Bio::EnsEMBL::Translation::TO_JSON {
  my ($self) = @_;
  my $shallow = Bio::EnsEMBL::Storable::TO_JSON($self);
  delete $shallow->{transcript} if exists $shallow->{transcript};
  return $shallow;
}

sub Bio::EnsEMBL::Gene::TO_JSON {
  my ($self) = @_;
  my $shallow = Bio::EnsEMBL::Storable::TO_JSON($self);
  $shallow->{canonical_transcript} = $shallow->{canonical_transcript}->stable_id() if exists $shallow->{canonical_transcript};
  return $shallow;
}

sub Bio::EnsEMBL::Exon::TO_JSON {
  my ($self) = @_;
  my $shallow = Bio::EnsEMBL::Storable::TO_JSON($self);
  delete $shallow->{_seq_cache} if exists $shallow->{_seq_cache};
  return $shallow; 
}

package Script;

use strict;
use warnings;

use Bio::EnsEMBL::Registry;
use Getopt::Long qw/:config no_ignore_case auto_version bundling_override/;
use Pod::Usage;
use JSON::XS;
use File::Spec::Functions qw/:ALL/;
use File::Path qw/make_path/;
use Bio::EnsEMBL::Utils::IO qw/work_with_file/;
use PerlIO::gzip;

sub run {
  my ($class) = @_;
  my $self = bless({}, $class);
  $self->args();
  $self->check();
  $self->setup();
  $self->process();
  return;
}

sub args {
  my ($self) = @_;
  my $opts = {
    port => 3306
  };
  GetOptions(
    $opts, qw/
      release|version=i
      host|hostname|h=s
      port|P=i
      username|user|u=s
      password|pass|p=s
      dir=s
      species=s
      gzip
      pretty
      verbose
      help
      man
      /
  ) or pod2usage(-verbose => 1, -exitval => 1);
  pod2usage(-verbose => 1, -exitval => 0) if $opts->{help};
  pod2usage(-verbose => 2, -exitval => 0) if $opts->{man};
  $self->{opts} = $opts;
  return;
}

sub opts {
  my ($self) = @_;
  return $self->{'opts'};
}

sub check {
  my ($self) = @_;
  my $o = $self->opts();

  my @required_params = qw/host username dir/;

  foreach my $r (@required_params) {
    if (!$o->{$r}) {
      pod2usage(
        -message => "-${r} has not been given at the command line but is a required parameter",
        -verbose => 1,
        -exitval => 1
      );
    }
  }
  
  if(-d $o->{dir}) {
    pod2usage(
      -message => "-dir exists. Choosing not to continue because I could overwrite existing data",
      -verbose => 1,
      -exitval => 1
    );
  }
  else {
    make_path($o->{dir});
    $o->{dir} = rel2abs($o->{dir});
  }
  
  return;
}

sub setup {
  my ($self) = @_;
  my $o = $self->opts();
  
  ##SETTING UP REGISTRY
  my %args = (
    -HOST => $o->{host}, 
    -PORT => $o->{port}, 
    -USER => $o->{username},
  );
  $args{-VERBOSE} = 1 if $o->{verbose};
  $args{-DB_VERSION} = $o->{release};
  $args{-PASS} = $o->{password} if $o->{password};
  Bio::EnsEMBL::Registry->no_version_check(1);
  my $loaded = Bio::EnsEMBL::Registry->load_registry_from_db(%args);
  
  return;
}

sub process {
  my ($self) = @_;
  my $dbas = $self->_get_core_dbs();
  while (my $dba = shift @{$dbas}) {
    $self->_process_dba($dba);
    $dba->dbc->disconnect_if_idle();
    last;
  }
  return;
}

sub _process_dba {
  my ($self, $dba) = @_;
  my $o = $self->opts();
  my ($ext, $mode) = ($o->{gzip}) ? ('json.gz', '>:gzip') : ('json', '>');
  my $species = $dba->get_MetaContainer()->get_production_name();
  my $release = $dba->get_MetaContainer()->get_schema_version();
  printf("Writing genes for %s\n", $species);
  my $json = JSON::XS->new()->convert_blessed();
  $json->pretty(1) if $o->{pretty};
  my $genes = $dba->get_GeneAdaptor()->fetch_all();
  while(my $gene = shift @{$genes}) {
    $gene->load(1);
    my $stable_id = $gene->stable_id();
    my ($extra_dir) = $stable_id =~ /(\w{2})$/;
    my $target_dir = catdir($o->{dir}, $release, $species, $extra_dir);
    make_path($target_dir);
    my $path = catfile($target_dir, sprintf("%s.%s", $gene->stable_id(), $ext));
    work_with_file($path, $mode, sub {
      my ($fh) = @_;
      print $fh $json->encode($gene);
      return;
    });
  }
  
  return;
}

sub _get_core_dbs {
  my ($self) = @_;
  my $dbas;
  if($self->opts->{species}) {
    $dbas = Bio::EnsEMBL::Registry->get_all_DBAdaptors(-SPECIES => $self->opts->{species});
  }
  else {
     $dbas = Bio::EnsEMBL::Registry->get_all_DBAdaptors(-GROUP => 'core');
  }
  my @final_dbas;
  while(my $dba = shift @{$dbas}) {
    next if $dba->species() eq 'multi';
    next if lc($dba->species()) eq 'ancestral sequences';
    next if $dba->dbc()->dbname() =~ /^.+_userdata$/xms;
    push(@final_dbas, $dba);
  }
  return [sort { $a->species() cmp $b->species() } @final_dbas];
}

Script->run();

1;
__END__

=pod

=head1 NAME

dump_ensembl_gene_json.pl

=head1 DESCRIPTION

Dump JSON file versions of all genes. Does EVERYTHING EVER.

=head1 OPTIONS

=over 8

=item B<--username | --user | -u>

REQUIRED. Username of the connecting account

=item B<--password | -pass | -p>

Password of the connecting user.

=item B<--release | --version>

REQUIRED. Indicates the release of Ensembl to process

=item B<--host | --host | -h>

REQUIRED. Host name of the database to connect to

=item B<--port | -P>

Optional integer of the database port. Defaults to 3306.

=item B<--dir>

  -dir /target/dir

REQUIRED. Target to write data to

=item B<--species>

Specify to run over a single species

=item B<--help>

Help message

=item B<--man>

Man page

=back

=head1 REQUIREMENTS

=over 8

=item Perl 5.8+

=item Bio::EnsEMBL

=back

=end

