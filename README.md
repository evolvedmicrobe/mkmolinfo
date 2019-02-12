## mkmolinfo

    mkmolinfo 1.0
    Nigel Delaney <nigel.delaney@10xgenomics.com>
    Takes a path to a CellRanger Output directory and creates an older (CellRanger < 3.0) style molecule_info.h5 file.    

    USAGE:
        mkmolinfo [OPTIONS] <OUT_DIRECTORY>    

    FLAGS:
        -h, --help       Prints help information
        -V, --version    Prints version information    

    OPTIONS:
        -o, --output <FILE>    Specify an output file name (default molecule_info.h5)    

    ARGS:
        <OUT_DIRECTORY>    Specify the output directory produced by CellRanger