## mkmolinfo

    mkmolinfo 1.0
    Nigel Delaney <nigel.delaney@10xgenomics.com>
    Takes a path to a Cell Ranger Output directory and creates an older (Cell Ranger < 3.0) style molecule_info.h5 file.    

    USAGE:
        mkmolinfo [FLAGS] [OPTIONS] <OUT_DIRECTORY>    

    FLAGS:
        -h, --help       Prints help information
        -v               Print progress as program parses through BAM file.
        -V, --version    Prints version information    

    OPTIONS:
        -o, --output <FILE>    Specify an output file name (default molecule_info_new.h5)    

    ARGS:
        <OUT_DIRECTORY>    Specify the output directory produced by Cell Ranger