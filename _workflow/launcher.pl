use strict;
use warnings;

# generic code invoked by a pipeline-specific wrapper script
# configures the environment and launches the master shell script

# working variables
use vars qw($pipelineDir $pipelineName $pipelineMessage @commands
            %workflowCommands %noParameterCommands $leftPad);
our ($command, @args)       = @ARGV;
$ENV{PIPELINE_DIR}          = $pipelineDir;
$ENV{PIPELINE_NAME}         = $pipelineName;
$ENV{PIPELINE_COMMAND}      = $command;
$ENV{PIPELINE_EXECUTABLE}   = $0;
$ENV{SLURP} = "$pipelineDir/_workflow/slurp";
our (%options, %longOptions, %shortOptions, %optionValues);

# load required launcher scripts
map { require "$pipelineDir/_workflow/launcher_$_.pl" } qw(help options utilities); 

# set the commands list and working command
our %commands = map { $$_{name} => $_ } @commands;
(!$command or $command eq '-h' or $command eq '--help') and showCommandHelp(undef, 1);
!$commands{$command} and showCommandHelp("unknown command: $command", 1);

# act on restricted commands (all end here)
if($command eq 'status'){
    loadCommandOptions();    # no command-specific options are needed
    getCommandLineOptions(); # but need data path
    exec "bash -c 'source $pipelineDir/_workflow/workflow.sh; showWorkflowStatus'";
} elsif($command eq 'rollback'){
    my ($targetCommand, $statusLevel);
    ($targetCommand, $statusLevel, @args) = @args;
    if(!defined $statusLevel){
        print  "\nusage: $pipelineName rollback <COMMAND> <LAST_SUCCESSFUL_STEP> <options>\n\n";
        exit;
    }
    loadCommandOptions($targetCommand);
    getCommandLineOptions();
    confirmAction("Pipeline status will be permanently reset.") or exit;
    $ENV{PIPELINE_COMMAND} = $targetCommand;
    $ENV{LAST_SUCCESSFUL_STEP} = $statusLevel;
    exec "bash -c 'source $pipelineDir/_workflow/workflow.sh; resetWorkflowStatus'";
} elsif($command eq 'options'){
    my ($targetCommand, $required) = @args;
    if(!defined $targetCommand){
        print  "\nusage: $pipelineName options <COMMAND> [required]\n\n";
        exit;
    }    
    loadCommandOptions($targetCommand); # need options but no values    
    my @optionsOut = sort { lc($$a{short}) cmp lc($$b{short}) or
                               $$a{short}  cmp    $$b{short} or
                               $$a{long}   cmp    $$b{long} } values %longOptions;
    foreach my $option(@optionsOut){
        if(!$required or $$option{required}){
            my $required = $$option{required} ? "*REQUIRED*" : "";         
            my $shortOut = $$option{short} ? "-$$option{short}" : "";
            print join("\t", $shortOut, "--$$option{long}", $required), "\n";
        }   
    }
    exit;
} elsif($command eq 'dependencies'){
    showDependencyHelp(); # no target command, no options or values
}

# set the options list for a pipeline action command
loadCommandOptions($command);
!$noParameterCommands{$command} and (!$args[0] or $args[0] eq '-h' or $args[0] eq '--help') and showOptionHelp();
getCommandLineOptions();

# check required options and fill defaults
foreach my $optionGroup(keys %options){
    foreach my $option(@{$options{$optionGroup}}){
        my $valueExists = defined $optionValues{$$option{long}};
        if($$option{required}){
            $valueExists or showOptionHelp("option '$$option{long}' is required");
        } elsif(!$valueExists and defined $$option{default}){ # options can carry 0 or zero-length strings
            $optionValues{$$option{long}} = $$option{default};
        }
    }
}

# load environment variables with provided values
foreach my $optionLong(keys %optionValues){
    setEnvVariable($optionLong, $optionValues{$optionLong});
}

# check for/create valid output paths
-d $ENV{OUTPUT_DIR} or showOptionHelp("directory does not exist: $ENV{OUTPUT_DIR}");
$ENV{LOGS_DIR} = "$ENV{OUTPUT_DIR}/$ENV{PIPELINE_NAME}_logs";
$ENV{LOGS_PREFIX} = "$ENV{LOGS_DIR}/$ENV{DATA_NAME}";
$ENV{DATA_PREFIX} = "$ENV{OUTPUT_DIR}/$ENV{DATA_NAME}";

# check memory requirements
$ENV{REQUIRED_RAM_INT} = getIntRam($commands{$command}{ram});
my $ramPerCpu = getIntRam($optionValues{'ram-per-cpu'});
$ENV{TOTAL_RAM_INT} = $ramPerCpu * $ENV{N_CPU};
if($ENV{TOTAL_RAM_INT} < $ENV{REQUIRED_RAM_INT}){
    showOptionHelp("insufficent net RAM, '$pipelineName $command' requires $commands{$command}{ram}");
}
$ENV{TOTAL_RAM} = getStrRam($ENV{TOTAL_RAM_INT});

# provide options and dependency feeback, for log file
unless($optionValues{quiet}){
    print "$pipelineName $command\n\n";    
    showOptionHelp(undef, 1, 1);
    print "\n";
    showDependencyHelp(1);
    print "\n";
    system "bash -c 'source $pipelineDir/_workflow/workflow.sh; showWorkflowStatus'";
    print "\n";
}

# launch pipeline inline and never return
unless($ENV{DRY_RUN}){
    &setDerivativeVariables and setDerivativeVariables();
    -d $ENV{LOGS_DIR} or mkdir $ENV{LOGS_DIR};
    $ENV{SCRIPT_TARGET} = $commands{$command}{script};
    exec "bash $pipelineDir/$ENV{SCRIPT_TARGET}";
}


