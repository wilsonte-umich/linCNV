use strict;
use warnings;

# universal command and option definitions that support all pipeline executables

use vars qw(@commands %options);

# non-specific pipeline action commands
push @commands, (
    {
        name    => "status",
        script  => "",
        message => "print the pipeline status file for the input data"
    },
    {
        name    => "rollback",
        script  => "",
        message => "revert the pipeline to an earlier step for the input data"
    },    
    {
        name    => "options",
        script  => "",
        message => "show a concise, alphabetically sorted option list for a command"
    },
    {
        name    => "dependencies",
        script  => "",
        message => "list and check non-standard applications called by the pipeline"
    },
);
our %workflowCommands = map { $_ => 1 } qw(status rollback options dependencies);
our %noParameterCommands = %workflowCommands;

# universal expected options
%options = ( 
    help => [
        {           
            long        => 'help',
            short       => 'h',
            type        => '',            
            required    => 0,
            default     => undef,
            message     => "show pipeline help",
        },
        {           
            long        => 'dry-run',
            short       => 'd',
            type        => '',            
            required    => 0,
            default     => undef,
            message     => "only show parsed variable values, do not execute the command",
        },
        {           
            long        => 'quiet',
            short       => 'q',
            type        => '',            
            required    => 0,
            default     => undef,
            message     => "suppress the variable and dependency feedback",
        },
    ],
    main => [
        {           
            long        => 'data-name',
            short       => 'n',
            type        => 'str',            
            required    => 1,
            default     => undef,
            message     => "simple name for the data (e.g. sample) being analyzed (no spaces)",
        },
        {           
            long        => 'output-dir',
            short       => 'o',
            type        => 'str',            
            required    => 1,
            default     => undef,
            message     => "the directory where output files will be placed; must already exist",
        },
        {           
            long        => 'force',
            short       => 'f',
            type        => '',            
            required    => 0,
            default     => undef,
            message     => "force various actions and outcomes (used by rollback and maybe pipeline)",
        }, 
    ],
    compute => [
        {           
            long        => 'n-cpu',
            short       => 'p',
            type        => 'int',            
            required    => 0,
            default     => 1,
            message     => "number of CPUs used for parallel processing",
        },      
        {           
            long        => 'ram-per-cpu',
            short       => 'r',
            type        => 'str',            
            required    => 0,
            default     => '4G',
            message     => "RAM allocated per CPU (e.g. 500M, 4G)",
        },
        {           
            long        => 'tmp-dir',
            short       => 't',
            type        => 'str',            
            required    => 0,
            default     => '/tmp',
            message     => "directory used for small temporary files (recommend SSD)",
        },
        {           
            long        => 'tmp-dir-large',
            short       => 'T',
            type        => 'str',            
            required    => 0,
            default     => '/tmp',
            message     => "directory used for large temporary files (generally >10GB)",
        },
    ],
);
our @topLevelOptionGroups = qw(help compute main);

