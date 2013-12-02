#!/usr/bin/perl -w
#
# Andrei Gaponenko, 2013

use File::Basename;
use File::Copy;
use Getopt::Long;

my $numEvents = 1000;
my $g4OutFile = 'hist_g4beam.root';
my $clark = '/home/andr/Clark/Clark';

sub usage() {
    return "
Usage: stopdist.pl g4template.fcl clark_config {flat|gauss} momentum widthParameter

";
}

GetOptions('numEvents=i' => \$numEvents) or die "Bad command line\n";

die usage() unless($#ARGV == 4);
my ($template, $clarkconf, $shape, $momentum, $width) = @ARGV;

my $workdir = basename($template, ('.fcl')) . '_' . $shape . sprintf("_m%.3f", $momentum) . sprintf("_w%.4f", $width);
mkdir $workdir or die "Can't create directory $workdir: $!\n";


#----------------------------------------------------------------
# Prepare G4 inputs
my $fcl = $workdir . '/' . basename($template);
copy($template, $fcl) or die "Can't copy $template => $fcl: $!\n";

open(my $fh, '>>', $fcl) or die "Can't append to $fcl: $!\n";
print $fh "physics.producers.generate.energySpec.spectrum: $shape\n";
if($shape eq "flat") {
    print $fh "physics.producers.generate.energySpec.center: $momentum\n";
    print $fh "physics.producers.generate.energySpec.halfWidth: $width\n";
}
elsif($shape eq 'gauss') {
    print $fh "physics.producers.generate.energySpec.mean: $momentum\n";
    print $fh "physics.producers.generate.energySpec.sigma: $width\n";
}
else {
    die "Unknown spectrum $shape\n";
}

#----------------------------------------------------------------
# Prepare Clark inputs

copy($clarkconf, $workdir ) or die "Can't copy $clarkconf => $workdir: $!\n";

#----------------------------------------------------------------
# Run the jobs

chdir $workdir or die "Can't cd to $workdir: $!\n";

system("/usr/bin/time mu2e -n $numEvents -c " . basename($fcl) . " > mucap.log 2>&1") == 0
    or die "Error running the G4 job\n";

system("/usr/bin/time $clark -l clark.out -o clark.root ". basename($clarkconf) ." $g4OutFile > clark.log 2>&1") == 0
    or die "Error running the Clark job\n";

#----------------------------------------------------------------
