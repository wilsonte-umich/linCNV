use strict;
use warnings;

# helper subs that support the launcher, broken out for code clarity

# working variables
use vars qw($pipelineDir @args %commands
            %options %optionValues %longOptions %shortOptions);

#------------------------------------------------------------------------------
# helper subs
#------------------------------------------------------------------------------
# load options for a specific command
sub loadCommandOptions {
    my ($command) = @_;
    if($command and $commands{$command}{require}){
        foreach my $codeGroup(@{$commands{$command}{require}}){ # add pipeline-specific options to universal options from launcher_options
            my $optionsFile = "$pipelineDir/$codeGroup/options.pl";
            -e $optionsFile and require $optionsFile;
        }        
    }
    foreach my $optionGroup(keys %options){ # create option lookups
        foreach my $option(@{$options{$optionGroup}}){
            $longOptions {$$option{long}}  = $option;
            $shortOptions{$$option{short} || ""} = $option; # variables can have only long-names
        }
    }    
}
# get command line option specifications
sub getCommandLineOptions {
    while (defined (my $optionList = shift @args)){
        defined $optionList or last; # no more options to process
        unless($optionList =~ m/^\-./){ # next item is a value, not an option
            unshift @args, $optionList;
            last;
        }
        if($optionList =~ m/^\-\-(.+)/){ # long option formatted request
            my $longOption = $1;
            defined $longOptions{$longOption} or showOptionHelp("'$longOption' is not a recognized option");
            setOption($longOptions{$longOption});
        } elsif ($optionList =~ m/^\-(.+)/){ # short option formatted request
            foreach my $shortOption(split('', $1)){
                defined $shortOptions{$shortOption} or showOptionHelp("'$shortOption' is not a recognized option");
                setOption($shortOptions{$shortOption});
            }   
        } else {
            showOptionHelp("malformed option list");
        }
    }     
}
# check and set option request
sub setOption { 
    my ($option) = @_;
    my $value = $$option{type} ? shift @args : 1; # type is false for boolean flag options
    (!defined $value or $value =~ m/^\-/) and showOptionHelp("missing value for option '$$option{long}'");
    $$option{type} ne 'int' or ($value =~ m|\D| and showOptionHelp("'$$option{name}' must be an integer"));
    $optionValues{$$option{long}} = $value;
    setEnvVariable($$option{long}, $value); # need here too for restricted commands    
    $$option{type} eq 'str' and $value =~ m/\s/ and $value = '"'.$value.'"';
    $ENV{CALL_VARIABLES_STRING} .= "--$$option{long} $value "; # used by workflow.sh for log reporting
}

# add option value to environment
sub setEnvVariable {
    my ($optionLong, $value) = @_;
    my $varName = uc($optionLong); # reformat option names (xxx-yyy) to variable names (XXX_YYY)
    $varName =~ s/-/_/g;
    $ENV{$varName} = $value;
}

# confirm dangerous actions
sub confirmAction {
    my ($message) = @_;
    $optionValues{force} and return 1;
    print "$message\n";
    print "Agree to continue? (y|n) ";
    my $response = <STDIN>;
    chomp $response;
    $response = (uc(substr($response, 0, 1)) eq "Y");
    $response or print "aborting with no action taken\n";
    return $response;
}

# convert string to integer RAM
sub getIntRam {
    my ($ramStr) = @_;
    my ($ram, $scale) = ($ramStr =~ m/^(\d+)(\w*)/);
    my %ramScales = (
        B => 1,
        K => 1e3,
        M => 1e6,
        G => 1e9
    );    
    $scale = $ramScales{uc($scale)};
    if($scale){
        $ram *= $scale
    } else {
        showOptionHelp("malformed RAM specification: $ramStr");
        exit;
    }
    return $ram;
}
sub getStrRam {
    my ($ramInt) = @_;
    if($ramInt < 1e3){ return $ramInt."B" }
    elsif($ramInt < 1e6){ return int($ramInt/1e3+0.5)."K" }
    elsif($ramInt < 1e9){ return int($ramInt/1e6+0.5)."M" }
    else{ return int($ramInt/1e9+0.5)."G" }
}

# check if a variable value is a valid directory or file
sub checkIsDirectory {
    my ($optionLong) = @_;
    $longOptions{$optionLong} or return;
    -d $optionValues{$optionLong} or showOptionHelp("'$optionLong' is not a directory: $optionValues{$optionLong}"); 
}
sub checkIsFile {
    my ($optionLong, $fileName) = @_;
    $optionLong and ($longOptions{$optionLong} or return);
    $fileName or $fileName = $optionValues{$optionLong};
    -e $fileName or showOptionHelp("file not found: $fileName"); 
}

1;

