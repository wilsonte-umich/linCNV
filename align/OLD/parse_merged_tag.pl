use strict;
use warnings;

while (my $line = <STDIN>) {
    if($line =~ m/^\@/){
        chomp $line;
        my @f = split(" ", $line);
        print join(":", $f[0], $f[$#f] =~ m/^merged/ ? 'mrg' : 'unmrg'), "\n";
    } else {
        print $line;
    } 
}
    
    