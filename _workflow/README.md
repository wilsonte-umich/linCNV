
.../_workflow

carries a library of scripts that support knitting a series
of worker scripts together into a single coherent pipeline that
can be called and shared as an executable.

launcher*.pl    'launcher' scripts are called by executables;
                they parse pipeline definitions and user options,
                load environment variables, and run pipeline scripts

workflow.*      'workflow' scripts are called by running pipeline scripts
                to support various common actions taken by pipelines

See scripts for more detail.

