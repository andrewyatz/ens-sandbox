package CrossMap;

# Process
# - Get binary, input file & type (gff, bed, vcf)
# - See if it has UCSC names; if so convert the names 1st before doing anything else (get the converting chain for the source assembly)
#   - Run CrossMap on this (write to a temp file)
# - Get source and target assembly - convert to chain file location
#   - Run CrossMap on this
# - Return the location of the file

use base qw/Bio::EnsEMBL::Hive::Process/;

sub param_defaults {
  return {
    ucsc_names => 0, # input has UCSC names and must be converted before anything else
    crossmap_binary => 'CrossMap', #location of CrossMap
    # species => '', #name of the species to process
    # input => '', #input file
    # output => '', #expected output file
    # chain_dir => '', #location of the chain files on disk
    # source_assembly => '', #assembly to convert from
    # target_assembly => '', #assembly to convert to
  };
}

sub fetch_input {
  my ($self) = @_;
  my $type = $self->param_required('type');
  my $type_count = grep { $type eq $_ } qw/bam bed gff gtf vcf wig bigwig/;
  $self->throw("Do not understand the type '${type}'") if ! $type_count;

  $self->param('type', 'gff') if $type eq 'gtf'; # type is just gff not gtf
  return 1;
}

sub run {
  my ($self) = @_;
  my $input = $self->param_required('input');
  my $new_input = $self->convert_UCSC($input);
  my $output = $self->convert($new_input);
  $self->param('output', $output);
  return 1;
}

sub convert_UCSC {
  my ($self, $input) = @_;
  return $input unless $self->param('ucsc_names');
  my $output = "${input}.ucsc";
  $self->run_CrossMap($self->chain_file('ucsc'), $input, $output);
  return $output;
}

sub convert {
  my ($self, $input) = @_;
  my $target_assembly = $self->param_required('target_assembly');
  my $output = ($self->param_is_defined('output')) ? $self->param('output') : "${input}.${target_assembly}";
  $self->run_CrossMap($self->chain_file($target_assembly), $input, $output);
  return $output;
}

# Assumes paths like
#   /path/to/chain_dir/species/sourceassembly/targetassembly.chain
#   /nfs/chains/homo_sapiens/GRCh37/GRCh38.chain
#   /nfs/chains/homo_sapiens/GRCh38/GRCh37.chain
#   /nfs/chains/homo_sapiens/GRCh37/ucsc.chain
sub chain_file {
  my ($self, $target_assembly) = @_;
  my $species = $self->param_required('species');
  my $source_assembly = $self->param_required('source_assembly');
  my $dir = File::Spec->catdir($self->param_required('chain_dir'), $species, $source_assembly);
  my $file = "${target_assembly}.chain";
  return File::Spec->catfile($dir, $file);
}

sub run_CrossMap {
  my ($self, $chain_file, $input_file, $output_file) = @_;
  $self->throw("Cannot do mapping. Chain file '${chain_file}' does not exist") if ! -f $chain_file;
  $self->throw("Cannot do mapping. Input file '${input_file}' does not exist") if ! -f $input_file;
  my $bin = $self->param_required('crossmap_binary');
  my $type = $self->param_required('type');
  # Command: CrossMap type chain input output
  my $cmd = "${bin} ${type} ${chain_file} ${input_file} ${output_file}";
  return $self->run_cmd($cmd);
}

sub run_cmd {
  my ($self, $cmd) = @_;
  warn("About to run '${cmd}'");
  my $output = `$cmd 2>&1`;
  my $rc = $? >> 8;
  $self->throw("Cannot run program '$cmd'. Return code was ${rc}. Program output was $output") if $rc;
  return ($rc, $output);
}

1;