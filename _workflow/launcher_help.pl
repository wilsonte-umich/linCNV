use strict;
use warnings;

# functions that provide executable help feedback on the command line

# working variables
use vars qw($pipelineName $pipelineMessage @commands %workflowCommands @dependencies
            %options @topLevelOptionGroups %optionValues);
my $commandTabLength = 15;
my $optionTabLength = 25;
our $leftPad = (" ") x 4;
my $usage = "$pipelineMessage\n\n".
             "usage\n$leftPad"."$pipelineName <command> <options>\n".
            "$leftPad"."$pipelineName <command> --help\n".
            "$leftPad"."$pipelineName --help\n";

sub showCommandHelp {
    my ($error, $exit) = @_;
    $error and print "\n".("!" x 60)."\n$error\n".("!" x 60)."\n";
    print  "\n$usage\n";
    addCommandGroup("pipeline-specific", 0);    
    addCommandGroup("general workflow",  1);
    $exit and exit;
    sub addCommandGroup {
        my ($group, $toggle) = @_;
        print "$group commands\n";
        foreach my $command(@commands){
            ($toggle xor $workflowCommands{$$command{name}}) and next;
            my $commandLength = length($$command{name});
            my $spaces = (" ") x ($commandTabLength - $commandLength);
            print  "$leftPad"."$$command{name}$spaces$$command{message}\n";
        }
        print  "\n";        
    }
}

sub showOptionHelp {
    my ($error, $useValues, $suppressExit) = @_;
    $useValues or showCommandHelp($error);
    my %groupSeen;
    foreach my $optionGroup(@topLevelOptionGroups, keys %options){
        $groupSeen{$optionGroup} and next;
        $groupSeen{$optionGroup} = 1;
        print  "$optionGroup options\n";
        foreach my $option(@{$options{$optionGroup}}){
            my $left = "-$$option{short},--$$option{long}";
            my $leftLength = length($left);
            my $nSpaces = $optionTabLength - $leftLength;
            my $spaces = (" ") x ($nSpaces > 0 ? $nSpaces : 0);
            if($useValues){
                defined $optionValues{$$option{long}} and 
                    print  "$leftPad"."$left$spaces$optionValues{$$option{long}}\n";            
            } else {
                my $type = $$option{type} ? "<$$option{type}> " : "";
                my $required = $$option{required} ? "*REQUIRED*" : ($$option{default} ? "[$$option{default}]" : "");                
                my $right = "$type$$option{message} $required";
                print  "$leftPad"."$left$spaces$right\n";                
            }
        }
        $useValues or print  "\n";
    }
    $suppressExit or exit;
}

sub showDependencyHelp {
    my ($suppressExit) = @_;
    foreach my $dependency(@dependencies){
        my $commandLength = length($$dependency{name});
        my $spaces = (" ") x ($commandTabLength - $commandLength);        
        my $path = qx(bash -c "command -v $$dependency{name}");
        print "$$dependency{name}$spaces".($path ? $path : "!!! MISSING !!!\n");
    }
    $suppressExit or exit;
}

1;

