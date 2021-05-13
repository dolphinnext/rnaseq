$HOSTNAME = ""
params.outdir = 'results'  

// enable required indexes to build them
params.use_Bowtie2_Index = (params.run_Sequential_Mapping == "yes" || params.run_Tophat == "yes") ? "yes" : ""
params.use_Bowtie_Index  = (params.run_Sequential_Mapping == "yes") ? "yes" : ""
params.use_STAR_Index    = (params.run_Sequential_Mapping == "yes" || params.run_STAR == "yes") ? "yes" : ""
params.use_Hisat2_Index  = (params.run_HISAT2 == "yes") ? "yes" : ""
params.use_RSEM_Index = (params.run_RSEM == "yes") ? "yes" : ""
params.use_Kallisto_Index = (params.run_Kallisto == "yes") ? "yes" : ""
params.nucleicAcidType = "rna"
if (!params.run_FeatureCounts_after_STAR){params.run_FeatureCounts_after_STAR = ""} 
if (!params.run_FeatureCounts_after_Hisat2){params.run_FeatureCounts_after_Hisat2 = ""} 
if (!params.run_FeatureCounts_after_Tophat2){params.run_FeatureCounts_after_Tophat2 = ""} 
if (!params.run_FeatureCounts_after_RSEM){params.run_FeatureCounts_after_RSEM = ""} 
if (!params.run_FeatureCounts_after_Kallisto){params.run_FeatureCounts_after_Kallisto = ""} 
if (!params.mate){params.mate = ""} 
if (!params.reads){params.reads = ""} 

Channel.value(params.run_FeatureCounts_after_STAR).set{g_179_run_featureCounts_g221_125}
Channel.value(params.run_FeatureCounts_after_Hisat2).set{g_188_run_featureCounts_g220_125}
Channel.value(params.run_FeatureCounts_after_Tophat2).set{g_189_run_featureCounts_g222_125}
Channel.value(params.run_FeatureCounts_after_RSEM).set{g_203_run_featureCounts_g219_125}
Channel.value(params.run_FeatureCounts_after_Kallisto).set{g_225_run_featureCounts_g224_125}
Channel.value(params.mate).into{g_229_mate_g213_3;g_229_mate_g213_11;g_229_mate_g213_16;g_229_mate_g213_18;g_229_mate_g213_19;g_229_mate_g213_20;g_229_mate_g213_21;g_229_mate_g216_11;g_229_mate_g217_16;g_229_mate_g218_3;g_229_mate_g218_11;g_229_mate_g219_82;g_229_mate_g219_95;g_229_mate_g219_123;g_229_mate_g219_126;g_229_mate_g220_82;g_229_mate_g220_95;g_229_mate_g220_123;g_229_mate_g220_126;g_229_mate_g221_82;g_229_mate_g221_95;g_229_mate_g221_123;g_229_mate_g221_126;g_229_mate_g222_82;g_229_mate_g222_95;g_229_mate_g222_123;g_229_mate_g222_126;g_229_mate_g224_82;g_229_mate_g224_95;g_229_mate_g224_123;g_229_mate_g224_126;g_229_mate_g240_21;g_229_mate_g241_26;g_229_mate_g241_30;g_229_mate_g241_34;g_229_mate_g242_22}
Channel
	.fromFilePairs( params.reads , size: params.mate == "single" ? 1 : params.mate == "pair" ? 2 : params.mate == "triple" ? 3 : params.mate == "quadruple" ? 4 : -1 )
	.ifEmpty { error "Cannot find any reads matching: ${params.reads}" }
	.into{g_230_reads_g213_3;g_230_reads_g213_18}


//* params.run_Adapter_Removal =   "no"   //* @dropdown @options:"yes","no" @show_settings:"Adapter_Removal"
//* @style @multicolumn:{seed_mismatches, palindrome_clip_threshold, simple_clip_threshold} @condition:{Tool_for_Adapter_Removal="trimmomatic", seed_mismatches, palindrome_clip_threshold, simple_clip_threshold}, {Tool_for_Adapter_Removal="fastx_clipper", discard_non_clipped}

//* autofill
//* platform
if ($HOSTNAME == "ghpcc06.umassrc.org"){
    $TIME = 1000
    $CPU  = 1
    $MEMORY = 24
    $QUEUE = "long"
}
//* platform
//* autofill
if (!((params.run_Adapter_Removal && (params.run_Adapter_Removal == "yes")) || !params.run_Adapter_Removal)){
g_230_reads_g213_18.into{g213_18_reads_g213_19}
g213_18_log_file_g213_11 = Channel.empty()
} else {


process Adapter_Trimmer_Quality_Module_Adapter_Removal {

input:
 set val(name), file(reads) from g_230_reads_g213_18
 val mate from g_229_mate_g213_18

output:
 set val(name), file("reads/*.fastq")  into g213_18_reads_g213_19
 file "*.{fastx,trimmomatic}.log"  into g213_18_log_file_g213_11

errorStrategy 'retry'

when:
(params.run_Adapter_Removal && (params.run_Adapter_Removal == "yes")) || !params.run_Adapter_Removal

shell:
phred = params.Adapter_Trimmer_Quality_Module_Adapter_Removal.phred
Tool_for_Adapter_Removal = params.Adapter_Trimmer_Quality_Module_Adapter_Removal.Tool_for_Adapter_Removal
Adapter_Sequence = params.Adapter_Trimmer_Quality_Module_Adapter_Removal.Adapter_Sequence
//trimmomatic_inputs
min_length = params.Adapter_Trimmer_Quality_Module_Adapter_Removal.min_length
seed_mismatches = params.Adapter_Trimmer_Quality_Module_Adapter_Removal.seed_mismatches
palindrome_clip_threshold = params.Adapter_Trimmer_Quality_Module_Adapter_Removal.palindrome_clip_threshold
simple_clip_threshold = params.Adapter_Trimmer_Quality_Module_Adapter_Removal.simple_clip_threshold

//fastx_clipper_inputs
discard_non_clipped = params.Adapter_Trimmer_Quality_Module_Adapter_Removal.discard_non_clipped
    
remove_previous_reads = params.Adapter_Trimmer_Quality_Module_Adapter_Removal.remove_previous_reads
discard_non_clipped_text = ""
if (discard_non_clipped == "yes") {discard_non_clipped_text = "-c"}
nameAll = reads.toString()
nameArray = nameAll.split(' ')
file2 = ""
if (nameAll.contains('.gz')) {
    newName =  nameArray[0] - ~/(\.fastq.gz)?(\.fq.gz)?$/
    file1 =  nameArray[0] - '.gz' 
    if (mate == "pair") {file2 =  nameArray[1] - '.gz'}
    runGzip = "ls *.gz | xargs -i echo gzip -df {} | sh"
} else {
    newName =  nameArray[0] - ~/(\.fastq)?(\.fq)?$/
    file1 =  nameArray[0]
    if (mate == "pair") {file2 =  nameArray[1]}
    runGzip = ''
}
'''
#!/usr/bin/env perl
 use List::Util qw[min max];
 use strict;
 use File::Basename;
 use Getopt::Long;
 use Pod::Usage;
 use Cwd qw();
 
runCmd("mkdir reads adapter unpaired");

open(OUT, ">adapter/adapter.fa");
my @adaps=split(/\n/,"!{Adapter_Sequence}");
my $i=1;
foreach my $adap (@adaps)
{
 print OUT ">adapter$i\\n$adap\\n";
 $i++;
}
close(OUT);

system("!{runGzip}");
my $quality="!{phred}";
print "fastq quality: $quality\\n";
print "tool: !{Tool_for_Adapter_Removal}\\n";

if ("!{mate}" eq "pair") {
    if ("!{Tool_for_Adapter_Removal}" eq "trimmomatic") {
        runCmd("trimmomatic PE -threads 1 -phred${quality} !{file1} !{file2} reads/!{name}.1.fastq unpaired/!{name}.1.fastq.unpaired reads/!{name}.2.fastq unpaired/!{name}.2.fastq.unpaired ILLUMINACLIP:adapter/adapter.fa:!{seed_mismatches}:!{palindrome_clip_threshold}:!{simple_clip_threshold} MINLEN:!{min_length} 2> !{name}.trimmomatic.log");
    } elsif ("!{Tool_for_Adapter_Removal}" eq "fastx_clipper") {
        print "Fastx_clipper is not suitable for paired reads.";
    }
} else {
    if ("!{Tool_for_Adapter_Removal}" eq "trimmomatic") {
        runCmd("trimmomatic SE -threads 1  -phred${quality} !{file1} reads/!{name}.fastq ILLUMINACLIP:adapter/adapter.fa:!{seed_mismatches}:!{palindrome_clip_threshold}:!{simple_clip_threshold} MINLEN:!{min_length} 2> !{name}.trimmomatic.log");
    } elsif ("!{Tool_for_Adapter_Removal}" eq "fastx_clipper") {
        runCmd("fastx_clipper  -Q $quality -a !{Adapter_Sequence} -l !{min_length} !{discard_non_clipped_text} -v -i !{file1} -o reads/!{name}.fastq > !{name}.fastx.log");
    }
}
if ("!{remove_previous_reads}" eq "true") {
    my $currpath = Cwd::cwd();
    my @paths = (split '/', $currpath);
    splice(@paths, -2);
    my $workdir= join '/', @paths;
    splice(@paths, -1);
    my $inputsdir = join '/', @paths;
    $inputsdir .= "/work";
    print "INFO: inputs reads will be removed if they are located in the $workdir $inputsdir\\n";
    my @listOfFiles = `readlink -e !{file1} !{file2}`;
    foreach my $targetFile (@listOfFiles){
        if (index($targetFile, $workdir) != -1 || index($targetFile, $inputsdir) != -1) {
            system("rm -f $targetFile");
            print "INFO: $targetFile deleted.\\n";
        }
    }
}


##Subroutines
sub runCmd {
    my ($com) = @_;
    my $error = system($com);
    if   ($error) { die "Command failed: $error $com\\n"; }
    else          { print "Command successful: $com\\n"; }
}
'''

}
}


//* params.run_Trimmer =   "no"   //* @dropdown @options:"yes","no" @show_settings:"Trimmer"
//* @style @multicolumn:{trim_length_5prime,trim_length_3prime}, {trim_length_5prime_R1,trim_length_3prime_R1}, {trim_length_5prime_R2,trim_length_3prime_R2} @condition:{single_or_paired_end_reads="single", trim_length_5prime,trim_length_3prime}, {single_or_paired_end_reads="pair", trim_length_5prime_R1,trim_length_3prime_R1,trim_length_5prime_R2,trim_length_3prime_R2}

//* autofill
//* platform
if ($HOSTNAME == "ghpcc06.umassrc.org"){
    $TIME = 500
    $CPU  = 1
    $MEMORY = 8
    $QUEUE = "long"
}
//* platform
//* autofill
if (!((params.run_Trimmer && (params.run_Trimmer == "yes")) || !params.run_Trimmer)){
g213_18_reads_g213_19.into{g213_19_reads_g213_20}
g213_19_log_file_g213_21 = Channel.empty()
} else {


process Adapter_Trimmer_Quality_Module_Trimmer {

input:
 set val(name), file(reads) from g213_18_reads_g213_19
 val mate from g_229_mate_g213_19

output:
 set val(name), file("reads/*q")  into g213_19_reads_g213_20
 file "*.log" optional true  into g213_19_log_file_g213_21

errorStrategy 'retry'

when:
(params.run_Trimmer && (params.run_Trimmer == "yes")) || !params.run_Trimmer

shell:
phred = params.Adapter_Trimmer_Quality_Module_Trimmer.phred
single_or_paired_end_reads = params.Adapter_Trimmer_Quality_Module_Trimmer.single_or_paired_end_reads
trim_length_5prime = params.Adapter_Trimmer_Quality_Module_Trimmer.trim_length_5prime
trim_length_3prime = params.Adapter_Trimmer_Quality_Module_Trimmer.trim_length_3prime
trim_length_5prime_R1 = params.Adapter_Trimmer_Quality_Module_Trimmer.trim_length_5prime_R1
trim_length_3prime_R1 = params.Adapter_Trimmer_Quality_Module_Trimmer.trim_length_3prime_R1
trim_length_5prime_R2 = params.Adapter_Trimmer_Quality_Module_Trimmer.trim_length_5prime_R2
trim_length_3prime_R2 = params.Adapter_Trimmer_Quality_Module_Trimmer.trim_length_3prime_R2
remove_previous_reads = params.Adapter_Trimmer_Quality_Module_Trimmer.remove_previous_reads

nameAll = reads.toString()
nameArray = nameAll.split(' ')
file2 = ""
if (nameAll.contains('.gz')) {
    newName =  nameArray[0] - ~/(\.fastq.gz)?(\.fq.gz)?$/
    file1 =  nameArray[0] - '.gz' 
    if (mate == "pair") {file2 =  nameArray[1] - '.gz'}
    runGzip = "ls *.gz | xargs -i echo gzip -df {} | sh"
} else {
    newName =  nameArray[0] - ~/(\.fastq)?(\.fq)?$/
    file1 =  nameArray[0]
    if (mate == "pair") {file2 =  nameArray[1]}
    runGzip = ''
}
'''
#!/usr/bin/env perl
 use List::Util qw[min max];
 use strict;
 use File::Basename;
 use Getopt::Long;
 use Pod::Usage; 
 use Cwd qw();
 
system("mkdir reads");
system("!{runGzip}");
my $file1 = "";
my $file2 = "";
if ("!{mate}" eq "pair") {
    $file1 = "!{file1}";
    $file2 = "!{file2}";
    my $trim1 = "!{trim_length_5prime_R1}:!{trim_length_3prime_R1}";
    my $trim2 = "!{trim_length_5prime_R2}:!{trim_length_3prime_R2}";
    my $len=getLength($file1);
    print "length of $file1: $len\\n";
    trimFiles($file1, $trim1, $len);
    my $len=getLength($file2);
    print "INFO: length of $file2: $len\\n";
    trimFiles($file2, $trim2, $len);
} else {
    $file1 = "!{file1}";
    my $trim1 = "!{trim_length_5prime}:!{trim_length_3prime}";
    my $len=getLength($file1);
    print "INFO: length of file1: $len\\n";
    trimFiles($file1, $trim1, $len);
}
if ("!{remove_previous_reads}" eq "true") {
    my $currpath = Cwd::cwd();
    my @paths = (split '/', $currpath);
    splice(@paths, -2);
    my $workdir= join '/', @paths;
    splice(@paths, -1);
    my $inputsdir= join '/', @paths;
    $inputsdir .= "/inputs";
    print "INFO: inputs reads will be removed if they are located in the workdir inputsdir\\n";
    my @listOfFiles = `readlink -e !{file1} !{file2}`;
    foreach my $targetFile (@listOfFiles){
        if (index($targetFile, $workdir) != -1 || index($targetFile, $inputsdir) != -1) {
            system("rm -f $targetFile");
            print "INFO: $targetFile deleted.\\n";
        }
    }
}



sub trimFiles
{
  my ($file, $trim, $len)=@_;
    my @nts=split(/[,:\\s\\t]+/,$trim);
    my $inpfile="";
    my $com="";
    my $i=1;
    my $outfile="";
    my $param="";
    my $quality="-Q!{phred}";

    if (scalar(@nts)==2)
    {
      $param = "-f ".($nts[0]+1) if (exists($nts[0]) && $nts[0] >= 0 );
      $param .= " -l ".($len-$nts[1]) if (exists($nts[0]) && $nts[1] > 0 );
      $outfile="reads/$file";  
      $com="fastx_trimmer $quality -v $param -o $outfile -i $file > !{name}.fastx_trimmer.log" if ((exists($nts[0]) && $nts[0] > 0) || (exists($nts[0]) && $nts[1] > 0 ));
      print "INFO: $com\\n";
      if ($com eq ""){
          print "INFO: Trimmer skipped for $file \\n";
          system("mv $file reads/.");
      } else {
          runCmd("$com");
          print "INFO: Trimmer executed for $file \\n";
      }
    }

    
}


sub getLength
{
   my ($filename)=@_;
   open (IN, $filename);
   my $j=1;
   my $len=0;
   while(my $line=<IN>)
   {
     chomp($line);
     if ($j >50) { last;}
     if ($j%4==0)
     {
        $len=length($line);
     }
     $j++;
   }
   close(IN);
   return $len;
}

sub runCmd {
    my ($com) = @_;
    my $error = system($com);
    if   ($error) { die "Command failed: $error $com\\n"; }
    else          { print "Command successful: $com\\n"; }
}

'''

}
}



process Adapter_Trimmer_Quality_Module_Trimmer_Removal_Summary {

input:
 file logfile from g213_19_log_file_g213_21.collect()
 val mate from g_229_mate_g213_21

output:
 file "trimmer_summary.tsv"  into g213_21_outputFileTSV_g_198

shell:
'''
#!/usr/bin/env perl
use List::Util qw[min max];
use strict;
use File::Basename;
use Getopt::Long;
use Pod::Usage;
use Data::Dumper;

my @header;
my %all_files;
my %tsv;
my %tsvDetail;
my %headerHash;
my %headerText;
my %headerTextDetail;

my $i = 0;
chomp( my $contents = `ls *.log` );

my @files = split( /[\\n]+/, $contents );
foreach my $file (@files) {
    $i++;
    my $mapOrder = "1";
    if ($file =~ /(.*)\\.fastx_trimmer\\.log/){
        $file =~ /(.*)\\.fastx_trimmer\\.log/;
        my $mapper   = "fastx_trimmer";
        my $name = $1;    ##sample name
        push( @header, $mapper );
        my $in;
        my $out;
        chomp( $in =`cat $file | grep 'Input:' | awk '{sum+=\\$2} END {print sum}'` );
        chomp( $out =`cat $file | grep 'Output:' | awk '{sum+=\\$2} END {print sum}'` );

        $tsv{$name}{$mapper} = [ $in, $out ];
        $headerHash{$mapOrder} = $mapper;
        $headerText{$mapOrder} = [ "Total Reads", "Reads After Trimmer" ];
    }
}

my @mapOrderArray = ( keys %headerHash );
my @sortedOrderArray = sort { $a <=> $b } @mapOrderArray;

my $summary          = "trimmer_summary.tsv";
writeFile( $summary,          \\%headerText,       \\%tsv );

sub writeFile {
    my $summary    = $_[0];
    my %headerText = %{ $_[1] };
    my %tsv        = %{ $_[2] };
    open( OUT, ">$summary" );
    print OUT "Sample\\t";
    my @headArr = ();
    for my $mapOrder (@sortedOrderArray) {
        push( @headArr, @{ $headerText{$mapOrder} } );
    }
    my $headArrAll = join( "\\t", @headArr );
    print OUT "$headArrAll\\n";

    foreach my $name ( keys %tsv ) {
        my @rowArr = ();
        for my $mapOrder (@sortedOrderArray) {
            push( @rowArr, @{ $tsv{$name}{ $headerHash{$mapOrder} } } );
        }
        my $rowArrAll = join( "\\t", @rowArr );
        print OUT "$name\\t$rowArrAll\\n";
    }
    close(OUT);
}

'''
}

//* params.run_Quality_Filtering =   "no"   //* @dropdown @options:"yes","no" @show_settings:"Quality_Filtering"
//* @style @multicolumn:{window_size,required_quality}, {leading,trailing,minlen}, {minQuality,minPercent} @condition:{tool="trimmomatic", minlen, trailing, leading, required_quality_for_window_trimming, window_size}, {tool="fastx", minQuality, minPercent}

//* autofill
//* platform
if ($HOSTNAME == "ghpcc06.umassrc.org"){
    $TIME = 240
    $CPU  = 1
    $MEMORY = 8
    $QUEUE = "short"
}
//* platform
//* autofill
if (!((params.run_Quality_Filtering && (params.run_Quality_Filtering == "yes")) || !params.run_Quality_Filtering)){
g213_19_reads_g213_20.into{g213_20_reads_g241_34}
g213_20_log_file_g213_16 = Channel.empty()
} else {


process Adapter_Trimmer_Quality_Module_Quality_Filtering {

input:
 set val(name), file(reads) from g213_19_reads_g213_20
 val mate from g_229_mate_g213_20

output:
 set val(name), file("reads/*q")  into g213_20_reads_g241_34
 file "*.{fastx,trimmomatic}_quality.log" optional true  into g213_20_log_file_g213_16

errorStrategy 'retry'

when:
(params.run_Quality_Filtering && (params.run_Quality_Filtering == "yes")) || !params.run_Quality_Filtering    

shell:
tool = params.Adapter_Trimmer_Quality_Module_Quality_Filtering.tool
phred = params.Adapter_Trimmer_Quality_Module_Quality_Filtering.phred
window_size = params.Adapter_Trimmer_Quality_Module_Quality_Filtering.window_size
required_quality_for_window_trimming = params.Adapter_Trimmer_Quality_Module_Quality_Filtering.required_quality_for_window_trimming
leading = params.Adapter_Trimmer_Quality_Module_Quality_Filtering.leading
trailing = params.Adapter_Trimmer_Quality_Module_Quality_Filtering.trailing
minlen = params.Adapter_Trimmer_Quality_Module_Quality_Filtering.minlen


// fastx parameters
minQuality = params.Adapter_Trimmer_Quality_Module_Quality_Filtering.minQuality
minPercent = params.Adapter_Trimmer_Quality_Module_Quality_Filtering.minPercent

remove_previous_reads = params.Adapter_Trimmer_Quality_Module_Quality_Filtering.remove_previous_reads

nameAll = reads.toString()
nameArray = nameAll.split(' ')
file2 ="";
if (nameAll.contains('.gz')) {
    newName =  nameArray[0] - ~/(\.fastq.gz)?(\.fq.gz)?$/
    file1 =  nameArray[0] - '.gz' 
    if (mate == "pair") {file2 =  nameArray[1] - '.gz'}
    runGzip = "ls *.gz | xargs -i echo gzip -df {} | sh"
} else {
    newName =  nameArray[0] - ~/(\.fastq)?(\.fq)?$/
    file1 =  nameArray[0]
    if (mate == "pair") {file2 =  nameArray[1]}
    runGzip = ''
}
'''
#!/usr/bin/env perl
 use List::Util qw[min max];
 use strict;
 use File::Basename;
 use Getopt::Long;
 use Pod::Usage; 
 use Cwd qw();
 
system("mkdir reads unpaired");
system("!{runGzip}");
my $param = "SLIDINGWINDOW:"."!{window_size}".":"."!{required_quality_for_window_trimming}";
$param.=" LEADING:"."!{leading}";
$param.=" TRAILING:"."!{trailing}";
$param.=" MINLEN:"."!{minlen}";

my $quality="!{phred}";

print "INFO: fastq quality: $quality\\n";
     
if ("!{tool}" eq "trimmomatic") {
    if ("!{mate}" eq "pair") {
        runCmd("trimmomatic PE -phred${quality} !{file1} !{file2} reads/!{name}.1.fastq unpaired/!{name}.1.fastq.unpaired reads/!{name}.2.fastq unpaired/!{name}.1.fastq.unpaired $param 2> !{name}.trimmomatic_quality.log");
    } else {
        runCmd("trimmomatic SE -phred${quality} !{file1} reads/!{name}.fastq $param 2> !{name}.trimmomatic_quality.log");
    }
} elsif ("!{tool}" eq "fastx") {
    if ("!{mate}" eq "pair") {
        print("WARNING: Fastx option is not suitable for paired reads. This step will be skipped.");
        system("mv !{file1} !{file2} reads/.");
    } else {
        runCmd("fastq_quality_filter  -Q $quality -q !{minQuality} -p !{minPercent} -v -i !{file1} -o reads/!{name}.fastq > !{name}.fastx_quality.log");
    }
}
if ("!{remove_previous_reads}" eq "true") {
    my $currpath = Cwd::cwd();
    my @paths = (split '/', $currpath);
    splice(@paths, -2);
    my $workdir= join '/', @paths;
    splice(@paths, -1);
    my $inputsdir= join '/', @paths;
    $inputsdir .= "/inputs";
    print "INFO: inputs reads will be removed if they are located in the workdir inputsdir\\n";
    my @listOfFiles = `readlink -e !{file1} !{file2}`;
    foreach my $targetFile (@listOfFiles){
        if (index($targetFile, $workdir) != -1 || index($targetFile, $inputsdir) != -1) {
            system("rm -f $targetFile");
            print "INFO: $targetFile deleted.\\n";
        }
    }
}

##Subroutines
sub runCmd {
    my ($com) = @_;
    my $error = system($com);
    if   ($error) { die "Command failed: $error $com\\n"; }
    else          { print "Command successful: $com\\n"; }
}


'''

}
}



process Adapter_Trimmer_Quality_Module_Quality_Filtering_Summary {

input:
 file logfile from g213_20_log_file_g213_16.collect()
 val mate from g_229_mate_g213_16

output:
 file "quality_filter_summary.tsv"  into g213_16_outputFileTSV_g_198

shell:
'''
#!/usr/bin/env perl
use List::Util qw[min max];
use strict;
use File::Basename;
use Getopt::Long;
use Pod::Usage;
use Data::Dumper;

my @header;
my %all_files;
my %tsv;
my %headerHash;
my %headerText;

my $i = 0;
chomp( my $contents = `ls *.log` );
my @files = split( /[\\n]+/, $contents );
foreach my $file (@files) {
    $i++;
    my $mapper   = "";
    my $mapOrder = "1";
    if ($file =~ /(.*)\\.fastx_quality\\.log/){
        $mapper   = "fastx";
        $file =~ /(.*)\\.fastx_quality\\.log/;
        my $name = $1;    ##sample name
        push( @header, $mapper );
        my $in;
        my $out;
        chomp( $in =`cat $file | grep 'Input:' | awk '{sum+=\\$2} END {print sum}'` );
        chomp( $out =`cat $file | grep 'Output:' | awk '{sum+=\\$2} END {print sum}'` );
        $tsv{$name}{$mapper} = [ $in, $out ];
        $headerHash{$mapOrder} = $mapper;
        $headerText{$mapOrder} = [ "Total Reads", "Reads After Quality Filtering" ];
    } elsif ($file =~ /(.*)\\.trimmomatic_quality\\.log/){
        $mapper   = "trimmomatic";
        $file =~ /(.*)\\.trimmomatic_quality\\.log/;
        my $name = $1;    ##sample name
        push( @header, $mapper );
        my $in;
        my $out;
        if ( "!{mate}" eq "pair"){
            chomp( $in =`cat $file | grep 'Input Read Pairs:' | awk '{sum+=\\$4} END {print sum}'` );
            chomp( $out =`cat $file | grep 'Input Read Pairs:' | awk '{sum+=\\$7} END {print sum}'` );
        } else {
            chomp( $in =`cat $file | grep 'Input Reads:' | awk '{sum+=\\$3} END {print sum}'` );
            chomp( $out =`cat $file | grep 'Input Reads:' | awk '{sum+=\\$5} END {print sum}'` );
        }
        $tsv{$name}{$mapper} = [ $in, $out ];
        $headerHash{$mapOrder} = $mapper;
        $headerText{$mapOrder} = [ "Total Reads", "Reads After Quality Filtering" ];
    }
    
}

my @mapOrderArray = ( keys %headerHash );
my @sortedOrderArray = sort { $a <=> $b } @mapOrderArray;

my $summary          = "quality_filter_summary.tsv";
writeFile( $summary,          \\%headerText,       \\%tsv );

sub writeFile {
    my $summary    = $_[0];
    my %headerText = %{ $_[1] };
    my %tsv        = %{ $_[2] };
    open( OUT, ">$summary" );
    print OUT "Sample\\t";
    my @headArr = ();
    for my $mapOrder (@sortedOrderArray) {
        push( @headArr, @{ $headerText{$mapOrder} } );
    }
    my $headArrAll = join( "\\t", @headArr );
    print OUT "$headArrAll\\n";

    foreach my $name ( keys %tsv ) {
        my @rowArr = ();
        for my $mapOrder (@sortedOrderArray) {
            push( @rowArr, @{ $tsv{$name}{ $headerHash{$mapOrder} } } );
        }
        my $rowArrAll = join( "\\t", @rowArr );
        print OUT "$name\\t$rowArrAll\\n";
    }
    close(OUT);
}

'''
}


process Adapter_Trimmer_Quality_Module_Adapter_Removal_Summary {

publishDir params.outdir, overwrite: true, mode: 'copy',
	saveAs: {filename ->
	if (filename =~ /adapter_removal_detailed_summary.tsv$/) "adapter_removal_detailed_summary/$filename"
}

input:
 file logfile from g213_18_log_file_g213_11.collect()
 val mate from g_229_mate_g213_11

output:
 file "adapter_removal_summary.tsv"  into g213_11_outputFileTSV_g_198
 file "adapter_removal_detailed_summary.tsv" optional true  into g213_11_outputFile

shell:
'''
#!/usr/bin/env perl
use List::Util qw[min max];
use strict;
use File::Basename;
use Getopt::Long;
use Pod::Usage;
use Data::Dumper;

my @header;
my %all_files;
my %tsv;
my %tsvDetail;
my %headerHash;
my %headerText;
my %headerTextDetail;

my $i = 0;
chomp( my $contents = `ls *.log` );

my @files = split( /[\\n]+/, $contents );
foreach my $file (@files) {
    $i++;
    my $mapOrder = "1";
    if ($file =~ /(.*)\\.fastx\\.log/){
        $file =~ /(.*)\\.fastx\\.log/;
        my $mapper   = "fastx";
        my $name = $1;    ##sample name
        push( @header, $mapper );

        my $in;
        my $out;
        my $tooshort;
        my $adapteronly;
        my $noncliped;
        my $Nreads;

        chomp( $in =`cat $file | grep 'Input:' | awk '{sum+=\\$2} END {print sum}'` );
        chomp( $out =`cat $file | grep 'Output:' | awk '{sum+=\\$2} END {print sum}'` );
        chomp( $tooshort =`cat $file | grep 'too-short reads' | awk '{sum+=\\$2} END {print sum}'`);
        chomp( $adapteronly =`cat $file | grep 'adapter-only reads' | awk '{sum+=\\$2} END {print sum}'`);
        chomp( $noncliped =`cat $file | grep 'non-clipped reads.' | awk '{sum+=\\$2} END {print sum}'`);
        chomp( $Nreads =`cat $file | grep 'N reads.' | awk '{sum+=\\$2} END {print sum}'` );

        $tsv{$name}{$mapper} = [ $in, $out ];
        $headerHash{$mapOrder} = $mapper;
        $headerText{$mapOrder} = [ "Total Reads", "Reads After Adapter Removal" ];
        $tsvDetail{$name}{$mapper} = [ $in, $tooshort, $adapteronly, $noncliped, $Nreads, $out ];
        $headerTextDetail{$mapOrder} = ["Total Reads","Too-short reads","Adapter-only reads","Non-clipped reads","N reads","Reads After Adapter Removal"];
    } elsif ($file =~ /(.*)\\.trimmomatic\\.log/){
        $file =~ /(.*)\\.trimmomatic\\.log/;
        my $mapper   = "trimmomatic";
        my $name = $1;    ##sample name
        push( @header, $mapper );
        
        my $in;
        my $out;

        if ( "!{mate}" eq "pair"){
            chomp( $in =`cat $file | grep 'Input Read Pairs:' | awk '{sum+=\\$4} END {print sum}'` );
            chomp( $out =`cat $file | grep 'Input Read Pairs:' | awk '{sum+=\\$7} END {print sum}'` );
        } else {
            chomp( $in =`cat $file | grep 'Input Reads:' | awk '{sum+=\\$3} END {print sum}'` );
            chomp( $out =`cat $file | grep 'Input Reads:' | awk '{sum+=\\$5} END {print sum}'` );
        }
        


        $tsv{$name}{$mapper} = [ $in, $out ];
        $headerHash{$mapOrder} = $mapper;
        $headerText{$mapOrder} = [ "Total Reads", "Reads After Adapter Removal" ];
        
    }
    
}

my @mapOrderArray = ( keys %headerHash );
my @sortedOrderArray = sort { $a <=> $b } @mapOrderArray;

my $summary          = "adapter_removal_summary.tsv";
my $detailed_summary = "adapter_removal_detailed_summary.tsv";
writeFile( $summary,          \\%headerText,       \\%tsv );
if (%headerTextDetail){
    writeFile( $detailed_summary, \\%headerTextDetail, \\%tsvDetail );  
}

sub writeFile {
    my $summary    = $_[0];
    my %headerText = %{ $_[1] };
    my %tsv        = %{ $_[2] };
    open( OUT, ">$summary" );
    print OUT "Sample\\t";
    my @headArr = ();
    for my $mapOrder (@sortedOrderArray) {
        push( @headArr, @{ $headerText{$mapOrder} } );
    }
    my $headArrAll = join( "\\t", @headArr );
    print OUT "$headArrAll\\n";

    foreach my $name ( keys %tsv ) {
        my @rowArr = ();
        for my $mapOrder (@sortedOrderArray) {
            push( @rowArr, @{ $tsv{$name}{ $headerHash{$mapOrder} } } );
        }
        my $rowArrAll = join( "\\t", @rowArr );
        print OUT "$name\\t$rowArrAll\\n";
    }
    close(OUT);
}

'''
}

//* @style @array:{run_name,run_parameters} @multicolumn:{run_name,run_parameters}

process BAM_Analysis_Hisat2_featureCounts_Prep {

input:
 val run_featureCounts from g_188_run_featureCounts_g220_125

output:
 val run_params  into g220_125_run_parameters_g220_126

when:
run_featureCounts == "yes"

script:
run_name = params.BAM_Analysis_Hisat2_featureCounts_Prep.run_name
run_parameters = params.BAM_Analysis_Hisat2_featureCounts_Prep.run_parameters
sense_antisense = params.BAM_Analysis_Hisat2_featureCounts_Prep.sense_antisense

//define run_name and run_parameters in map item and push into run_params array
run_params = []
for (i = 0; i < run_parameters.size(); i++) {
   map = [:]
   map["run_name"] = run_name[i].replaceAll(" ","_").replaceAll(",","_").replaceAll(";","_").replaceAll("'","_").replaceAll('"',"_")
   map["run_parameters"] = run_parameters[i]
   run_params[i] = map
}
templateRunParams = run_parameters[0] ? run_parameters[0] : "-g gene_id -s 0 -Q 20 -T 2 -B -d 50 -D 1000 -C --fracOverlap 0 --minOverlap 1"
if (sense_antisense == "Yes"){
   map = [:]
   map["run_name"] = "gene_id_forward_expression"
   map["run_parameters"] = templateRunParams.replaceAll("transcript_id","gene_id").replaceAll("-s 0","-s 1").replaceAll("-s 2","-s 1")
   run_params.push(map)
   map = [:]
   map["run_name"] = "gene_id_reverse_expression"
   map["run_parameters"] = templateRunParams.replaceAll("transcript_id","gene_id").replaceAll("-s 0","-s 2").replaceAll("-s 1","-s 2")
   run_params.push(map)
   map = [:]
   map["run_name"] = "transcript_id_forward_expression"
   map["run_parameters"] = templateRunParams.replaceAll("gene_id","transcript_id").replaceAll("-s 0","-s 1").replaceAll("-s 2","-s 1")
   run_params.push(map)
   map = [:]
   map["run_name"] = "transcript_id_reverse_expression"
   map["run_parameters"] = templateRunParams.replaceAll("gene_id","transcript_id").replaceAll("-s 0","-s 2").replaceAll("-s 1","-s 2")
   run_params.push(map)
}
"""
"""

}

//* @style @array:{run_name,run_parameters} @multicolumn:{run_name,run_parameters}

process BAM_Analysis_STAR_featureCounts_Prep {

input:
 val run_featureCounts from g_179_run_featureCounts_g221_125

output:
 val run_params  into g221_125_run_parameters_g221_126

when:
run_featureCounts == "yes"

script:
run_name = params.BAM_Analysis_STAR_featureCounts_Prep.run_name
run_parameters = params.BAM_Analysis_STAR_featureCounts_Prep.run_parameters
sense_antisense = params.BAM_Analysis_STAR_featureCounts_Prep.sense_antisense

//define run_name and run_parameters in map item and push into run_params array
run_params = []
for (i = 0; i < run_parameters.size(); i++) {
   map = [:]
   map["run_name"] = run_name[i].replaceAll(" ","_").replaceAll(",","_").replaceAll(";","_").replaceAll("'","_").replaceAll('"',"_")
   map["run_parameters"] = run_parameters[i]
   run_params[i] = map
}
templateRunParams = run_parameters[0] ? run_parameters[0] : "-g gene_id -s 0 -Q 20 -T 2 -B -d 50 -D 1000 -C --fracOverlap 0 --minOverlap 1"
if (sense_antisense == "Yes"){
   map = [:]
   map["run_name"] = "gene_id_forward_expression"
   map["run_parameters"] = templateRunParams.replaceAll("transcript_id","gene_id").replaceAll("-s 0","-s 1").replaceAll("-s 2","-s 1")
   run_params.push(map)
   map = [:]
   map["run_name"] = "gene_id_reverse_expression"
   map["run_parameters"] = templateRunParams.replaceAll("transcript_id","gene_id").replaceAll("-s 0","-s 2").replaceAll("-s 1","-s 2")
   run_params.push(map)
   map = [:]
   map["run_name"] = "transcript_id_forward_expression"
   map["run_parameters"] = templateRunParams.replaceAll("gene_id","transcript_id").replaceAll("-s 0","-s 1").replaceAll("-s 2","-s 1")
   run_params.push(map)
   map = [:]
   map["run_name"] = "transcript_id_reverse_expression"
   map["run_parameters"] = templateRunParams.replaceAll("gene_id","transcript_id").replaceAll("-s 0","-s 2").replaceAll("-s 1","-s 2")
   run_params.push(map)
}
"""
"""

}

//* params.run_FastQC =  "no"  //* @dropdown @options:"yes","no" @description:"FastQC provides quality control checks on raw sequence data."



process Adapter_Trimmer_Quality_Module_FastQC {

publishDir params.outdir, overwrite: true, mode: 'copy',
	saveAs: {filename ->
	if (filename =~ /.*.(html|zip)$/) "fastqc/$filename"
}

input:
 val mate from g_229_mate_g213_3
 set val(name), file(reads) from g_230_reads_g213_3

output:
 file '*.{html,zip}'  into g213_3_FastQCout_g_177

errorStrategy 'retry'
maxRetries 3

when:
(params.run_FastQC && (params.run_FastQC == "yes"))

script:
nameAll = reads.toString()
if (nameAll.contains('.gz')) {
    file =  nameAll - '.gz' - '.gz'
    runGzip = "ls *.gz | xargs -i echo gzip -df {} | sh"
} else {
    file =  nameAll 
    runGzip = ''
}
"""
${runGzip}
fastqc ${file} 
"""
}

//* @style @array:{run_name,run_parameters} @multicolumn:{run_name,run_parameters}

process BAM_Analysis_Tophat2_featureCounts_Prep {

input:
 val run_featureCounts from g_189_run_featureCounts_g222_125

output:
 val run_params  into g222_125_run_parameters_g222_126

when:
run_featureCounts == "yes"

script:
run_name = params.BAM_Analysis_Tophat2_featureCounts_Prep.run_name
run_parameters = params.BAM_Analysis_Tophat2_featureCounts_Prep.run_parameters
sense_antisense = params.BAM_Analysis_Tophat2_featureCounts_Prep.sense_antisense

//define run_name and run_parameters in map item and push into run_params array
run_params = []
for (i = 0; i < run_parameters.size(); i++) {
   map = [:]
   map["run_name"] = run_name[i].replaceAll(" ","_").replaceAll(",","_").replaceAll(";","_").replaceAll("'","_").replaceAll('"',"_")
   map["run_parameters"] = run_parameters[i]
   run_params[i] = map
}
templateRunParams = run_parameters[0] ? run_parameters[0] : "-g gene_id -s 0 -Q 20 -T 2 -B -d 50 -D 1000 -C --fracOverlap 0 --minOverlap 1"
if (sense_antisense == "Yes"){
   map = [:]
   map["run_name"] = "gene_id_forward_expression"
   map["run_parameters"] = templateRunParams.replaceAll("transcript_id","gene_id").replaceAll("-s 0","-s 1").replaceAll("-s 2","-s 1")
   run_params.push(map)
   map = [:]
   map["run_name"] = "gene_id_reverse_expression"
   map["run_parameters"] = templateRunParams.replaceAll("transcript_id","gene_id").replaceAll("-s 0","-s 2").replaceAll("-s 1","-s 2")
   run_params.push(map)
   map = [:]
   map["run_name"] = "transcript_id_forward_expression"
   map["run_parameters"] = templateRunParams.replaceAll("gene_id","transcript_id").replaceAll("-s 0","-s 1").replaceAll("-s 2","-s 1")
   run_params.push(map)
   map = [:]
   map["run_name"] = "transcript_id_reverse_expression"
   map["run_parameters"] = templateRunParams.replaceAll("gene_id","transcript_id").replaceAll("-s 0","-s 2").replaceAll("-s 1","-s 2")
   run_params.push(map)
}
"""
"""

}

//* @style @array:{run_name,run_parameters} @multicolumn:{run_name,run_parameters}

process BAM_Analysis_RSEM_featureCounts_Prep {

input:
 val run_featureCounts from g_203_run_featureCounts_g219_125

output:
 val run_params  into g219_125_run_parameters_g219_126

when:
run_featureCounts == "yes"

script:
run_name = params.BAM_Analysis_RSEM_featureCounts_Prep.run_name
run_parameters = params.BAM_Analysis_RSEM_featureCounts_Prep.run_parameters
sense_antisense = params.BAM_Analysis_RSEM_featureCounts_Prep.sense_antisense

//define run_name and run_parameters in map item and push into run_params array
run_params = []
for (i = 0; i < run_parameters.size(); i++) {
   map = [:]
   map["run_name"] = run_name[i].replaceAll(" ","_").replaceAll(",","_").replaceAll(";","_").replaceAll("'","_").replaceAll('"',"_")
   map["run_parameters"] = run_parameters[i]
   run_params[i] = map
}
templateRunParams = run_parameters[0] ? run_parameters[0] : "-g gene_id -s 0 -Q 20 -T 2 -B -d 50 -D 1000 -C --fracOverlap 0 --minOverlap 1"
if (sense_antisense == "Yes"){
   map = [:]
   map["run_name"] = "gene_id_forward_expression"
   map["run_parameters"] = templateRunParams.replaceAll("transcript_id","gene_id").replaceAll("-s 0","-s 1").replaceAll("-s 2","-s 1")
   run_params.push(map)
   map = [:]
   map["run_name"] = "gene_id_reverse_expression"
   map["run_parameters"] = templateRunParams.replaceAll("transcript_id","gene_id").replaceAll("-s 0","-s 2").replaceAll("-s 1","-s 2")
   run_params.push(map)
   map = [:]
   map["run_name"] = "transcript_id_forward_expression"
   map["run_parameters"] = templateRunParams.replaceAll("gene_id","transcript_id").replaceAll("-s 0","-s 1").replaceAll("-s 2","-s 1")
   run_params.push(map)
   map = [:]
   map["run_name"] = "transcript_id_reverse_expression"
   map["run_parameters"] = templateRunParams.replaceAll("gene_id","transcript_id").replaceAll("-s 0","-s 2").replaceAll("-s 1","-s 2")
   run_params.push(map)
}
"""
"""

}

//* @style @array:{run_name,run_parameters} @multicolumn:{run_name,run_parameters}

process BAM_Analysis_Module_Kallisto_featureCounts_Prep {

input:
 val run_featureCounts from g_225_run_featureCounts_g224_125

output:
 val run_params  into g224_125_run_parameters_g224_126

when:
run_featureCounts == "yes"

script:
run_name = params.BAM_Analysis_Module_Kallisto_featureCounts_Prep.run_name
run_parameters = params.BAM_Analysis_Module_Kallisto_featureCounts_Prep.run_parameters
sense_antisense = params.BAM_Analysis_Module_Kallisto_featureCounts_Prep.sense_antisense

//define run_name and run_parameters in map item and push into run_params array
run_params = []
for (i = 0; i < run_parameters.size(); i++) {
   map = [:]
   map["run_name"] = run_name[i].replaceAll(" ","_").replaceAll(",","_").replaceAll(";","_").replaceAll("'","_").replaceAll('"',"_")
   map["run_parameters"] = run_parameters[i]
   run_params[i] = map
}
templateRunParams = run_parameters[0] ? run_parameters[0] : "-g gene_id -s 0 -Q 20 -T 2 -B -d 50 -D 1000 -C --fracOverlap 0 --minOverlap 1"
if (sense_antisense == "Yes"){
   map = [:]
   map["run_name"] = "gene_id_forward_expression"
   map["run_parameters"] = templateRunParams.replaceAll("transcript_id","gene_id").replaceAll("-s 0","-s 1").replaceAll("-s 2","-s 1")
   run_params.push(map)
   map = [:]
   map["run_name"] = "gene_id_reverse_expression"
   map["run_parameters"] = templateRunParams.replaceAll("transcript_id","gene_id").replaceAll("-s 0","-s 2").replaceAll("-s 1","-s 2")
   run_params.push(map)
   map = [:]
   map["run_name"] = "transcript_id_forward_expression"
   map["run_parameters"] = templateRunParams.replaceAll("gene_id","transcript_id").replaceAll("-s 0","-s 1").replaceAll("-s 2","-s 1")
   run_params.push(map)
   map = [:]
   map["run_name"] = "transcript_id_reverse_expression"
   map["run_parameters"] = templateRunParams.replaceAll("gene_id","transcript_id").replaceAll("-s 0","-s 2").replaceAll("-s 1","-s 2")
   run_params.push(map)
}
"""
"""

}

//* params.gtf =  ""  //* @input
//* params.genome =  ""  //* @input
//* params.commondb =  ""  //* @input
//* params.genome_url =  ""  //* @input
//* params.gtf_url =  ""  //* @input
//* params.commondb_url =  ""  //* @input

def downFile(path){
    if (path.take(1).indexOf("/") == 0){
      target=path
    } else {
      a=file(path)
      fname = a.getName().toString()
      target = "${baseDir}/work/${fname}"
      a.copyTo(target) 
    }
    return target
}

process Check_and_Build_Module_Check_Genome_GTF {


output:
 val "${params.genome}"  into g227_15_genomePath_g227_0, g227_15_genomePath_g227_6, g227_15_genomePath_g227_8, g227_15_genomePath_g227_10, g227_15_genomePath_g227_5, g227_15_genomePath_g227_13, g227_15_genomePath_g227_19
 val "${params.gtf}"  into g227_15_gtfPath_g227_0, g227_15_gtfPath_g227_6, g227_15_gtfPath_g227_8, g227_15_gtfPath_g227_10, g227_15_gtfPath_g227_4, g227_15_gtfPath_g227_13, g227_15_gtfPath_g227_19
 val "${params.commondb}"  into g227_15_commondb_path_g227_18

when:
params.run_checkAndBuild == "yes"

script:
gtf_dir  = params.gtf.substring(0, params.gtf.lastIndexOf('/')) 
genome_dir  = params.genome.substring(0, params.genome.lastIndexOf('/')) 
slashCount = params.commondb_url.count("/")
cutDir = slashCount - 3;

downGenomePath = ""
downGtfPath = ""
if ( !file("${params.genome}").exists() ) {
	downGenomePath=downFile(params.genome_url)
}
if ( !file("${params.gtf}").exists() ) {
	downGtfPath=downFile(params.gtf_url)
}

"""
if [ ! -e "${params.genome}" ] ; then
    echo "${params.genome} not found"
    mkdir -p ${genome_dir}
    cp -n $downGenomePath ${params.genome}
fi
if [ ! -e "${params.gtf}" ] ; then
    echo "${params.gtf} not found"
    mkdir -p ${gtf_dir}
    cp -n $downGtfPath ${params.gtf}
fi
if [ ! -e "${params.commondb}" ] ; then
    echo "${params.commondb} not found"
    mkdir -p ${params.commondb}
    wget -l inf -nc -nH --cut-dirs=$cutDir -R 'index.html*' -r --no-parent --directory-prefix=${params.commondb} ${params.commondb_url}
fi

"""




}


process Check_and_Build_Module_Kallisto_Index {

input:
 val genome from g227_15_genomePath_g227_19
 val gtf from g227_15_gtfPath_g227_19

output:
 val resultDir  into g227_19_genomeIndexPath_g240_21

when:
(params.use_Kallisto_Index == "yes") && (params.run_checkAndBuild == "yes")

script:
gtf_dir  = gtf.substring(0, gtf.lastIndexOf('/'))
indexbasedir  = gtf_dir.substring(0, gtf_dir.lastIndexOf('/'))
newDirName = "KallistoIndex"
resultDir = indexbasedir +"/"+ newDirName 
tmpResultDir = indexbasedir +"/_tmp_"+ newDirName

"""
if [ ! -e "${resultDir}/transcripts.idx" ] ; then
    echo "${resultDir}/transcripts.idx Kallisto index not found"
    rm -rf $tmpResultDir $resultDir && mkdir -p $tmpResultDir && cd $tmpResultDir
    filter_gtf_for_genes_in_genome.py --gtf ${params.gtf} --fasta ${params.genome} -o genome_filtered_genes.gtf
    gffread -F -w transcripts.fa -g ${params.genome} genome_filtered_genes.gtf
    gzip transcripts.fa
    kallisto index -i transcripts.idx transcripts.fa.gz
    cd .. && mv $tmpResultDir $resultDir 
fi
"""



}


process Check_and_Build_Module_Bowtie_Index {

input:
 val genome from g227_15_genomePath_g227_13
 val gtf from g227_15_gtfPath_g227_13

output:
 val resultDir  into g227_13_genomeIndexPath_g227_18

when:
(params.use_Bowtie_Index == "yes") && (params.run_checkAndBuild == "yes")

script:
bowtie_build_parameters = params.Check_and_Build_Module_Bowtie_Index.bowtie_build_parameters
basedir  = genome.substring(0, genome.lastIndexOf('/'))
gtf_dir  = gtf.substring(0, gtf.lastIndexOf('/')) 
basename = genome.substring(genome.lastIndexOf('/')+1,genome.lastIndexOf('.'))
filename = genome.substring(genome.lastIndexOf('/')+1,genome.length())
newDirName = "BowtieIndex"
resultDir = basedir.substring(0, basedir.lastIndexOf('/')) +"/"+ newDirName 
tmpResultDir = basedir.substring(0, basedir.lastIndexOf('/')) +"/_tmp_"+ newDirName 

"""
if [ ! -e "${resultDir}/${basename}.rev.2.ebwt" ] ; then
    echo "${resultDir}/${basename}.rev.2.ebwt Bowtie index not found"
    rm -rf $tmpResultDir $resultDir && mkdir -p $tmpResultDir && cd $tmpResultDir
    ln -s ../main/${filename} ${filename}
    bowtie-build ${bowtie_build_parameters} ${filename} ${basename}
    cd .. && mv $tmpResultDir $resultDir
fi
"""




}

//* params.gtf2bed_path =  ""  //* @input

process Check_and_Build_Module_Check_GTF2BED12 {

input:
 val gtf from g227_15_gtfPath_g227_4


when:
params.run_checkAndBuild == "yes"

script:
gtf_dir  = gtf.substring(0, gtf.lastIndexOf('/')) 
"""
if [ ! -e "${params.bed}" ] ; then
    echo "${params.bed} not found"
    perl ${params.gtf2bed_path} $gtf > ${params.bed}
fi
"""




}

//* params.gtf2bed_path =  ""  //* @input

process Check_and_Build_Module_Check_chrom_sizes_and_index {

input:
 val genome from g227_15_genomePath_g227_5


when:
params.run_checkAndBuild == "yes"

script:
genome_dir  = genome.substring(0, genome.lastIndexOf('/')) 
basename_and_path  = genome.substring(0, genome.lastIndexOf('.'))

"""
if [ ! -e "${params.genome_sizes}" ] ; then
    echo "${params.genome_sizes} not found"
    cat ${genome} | awk '\$0 ~ ">" {print c; c=0;printf substr(\$0,2,100) "\\t"; } \$0 !~ ">" {c+=length(\$0);} END { print c; }' > ${basename_and_path}.chrom.sizes
    ##clean first empty line
    sed -i '1{/^\$/d}' ${basename_and_path}.chrom.sizes
fi
"""




}


process Check_and_Build_Module_Check_Build_Rsem_Index {

input:
 val genome from g227_15_genomePath_g227_10
 val gtf from g227_15_gtfPath_g227_10

output:
 val cmdAr  into g227_10_command_g227_11
 val resultDirAr  into g227_10_path_g227_11

when:
(params.use_RSEM_Index == "yes") && (params.run_checkAndBuild == "yes")

script:
create_bowtie_rsem_index = params.Check_and_Build_Module_Check_Build_Rsem_Index.create_bowtie_rsem_index
create_bowtie2_rsem_index = params.Check_and_Build_Module_Check_Build_Rsem_Index.create_bowtie2_rsem_index
create_star_rsem_index = params.Check_and_Build_Module_Check_Build_Rsem_Index.create_star_rsem_index
transcript_to_gene_map = params.Check_and_Build_Module_Check_Build_Rsem_Index.transcript_to_gene_map
RSEM_build_parameters = params.Check_and_Build_Module_Check_Build_Rsem_Index.RSEM_build_parameters

genome_dir  = genome.substring(0, genome.lastIndexOf('/'))
gtf_dir  = gtf.substring(0, gtf.lastIndexOf('/')) 
basenameGenome = genome.substring(genome.lastIndexOf('/')+1,genome.lastIndexOf('.'))
newDirNameAr = []
cmdAr = []
resultDirAr = []
if (create_bowtie_rsem_index == "true"){ newDirNameAr.push('RSEM_ref_Bowtie') }
if (create_bowtie2_rsem_index == "true"){ newDirNameAr.push('RSEM_ref_Bowtie2') }
if (create_star_rsem_index == "true"){ newDirNameAr.push('RSEM_ref_STAR') }

transcript_to_gene_mapText = ""
if (transcript_to_gene_map?.trim()){
    transcript_to_gene_mapText = "--transcript-to-gene-map " + transcript_to_gene_map
}

for (i = 0; i < newDirNameAr.size(); i++) {
    resultDir = gtf_dir.substring(0, gtf_dir.lastIndexOf('/')) +"/"+ newDirNameAr[i]
    tmpResultDir = gtf_dir.substring(0, gtf_dir.lastIndexOf('/')) +"/_tmp_"+ newDirNameAr[i]
    resultDirAr.push(resultDir)
    cmd = ""
    indexType = ""
    if (newDirNameAr[i] == 'RSEM_ref_Bowtie'){
        indexType = "--bowtie "
        checkFile = "${basenameGenome}.rev.2.ebwt" 
    } else if (newDirNameAr[i] == 'RSEM_ref_Bowtie2'){
        indexType = "--bowtie2 "
        checkFile = "${basenameGenome}.rev.2.bt2" 
    } else if (newDirNameAr[i] == 'RSEM_ref_STAR'){
        indexType = "--star "
        checkFile = "genomeParameters.txt" 
    }
    cmd = "if [ ! -e \"${resultDir}/${checkFile}\" ] ; then rm -rf $tmpResultDir $resultDir && mkdir -p $tmpResultDir && cd $tmpResultDir && rsem-prepare-reference ${RSEM_build_parameters} --gtf ${gtf} ${transcript_to_gene_mapText} ${indexType} ${genome} ${basenameGenome} && cd .. && mv $tmpResultDir $resultDir; fi"
    cmdAr.push(cmd)
}


"""

"""

}


process Check_and_Build_Module_Build_Index_RSEM_run {

input:
 val resultDir from g227_10_path_g227_11.flatten()
 val command from g227_10_command_g227_11.flatten()

output:
 val resultDir  into g227_11_genomeIndexPath_g242_22

script:
"""    
$command
"""
}


process Check_and_Build_Module_Hisat2_Index {

input:
 val genome from g227_15_genomePath_g227_8
 val gtf from g227_15_gtfPath_g227_8

output:
 val resultDir  into g227_8_genomeIndexPath_g216_11

when:
(params.use_Hisat2_Index == "yes") && (params.run_checkAndBuild == "yes")

script:
hisat2_build_parameters = params.Check_and_Build_Module_Hisat2_Index.hisat2_build_parameters
genome_dir  = genome.substring(0, genome.lastIndexOf('/'))
gtf_dir  = gtf.substring(0, gtf.lastIndexOf('/')) 
basenameGenome = genome.substring(genome.lastIndexOf('/')+1,genome.lastIndexOf('.'))
basenameGTF = gtf.substring(gtf.lastIndexOf('/')+1,gtf.lastIndexOf('.'))
filename = genome.substring(genome.lastIndexOf('/')+1,genome.length())
newDirName = "Hisat2Index"
resultDir = gtf_dir.substring(0, gtf_dir.lastIndexOf('/')) +"/"+ newDirName 
tmpResultDir = gtf_dir.substring(0, gtf_dir.lastIndexOf('/')) +"/_tmp_"+ newDirName 

extract_splice_sites = "hisat2_extract_splice_sites.py ${gtf} > ${tmpResultDir}/${basenameGTF}.hisat2_splice_sites.txt"
extract_exons = "hisat2_extract_exons.py ${gtf}> ${tmpResultDir}/${basenameGTF}.hisat2_exons.txt"
ss = "--ss ${basenameGTF}.hisat2_splice_sites.txt"
exon = "--exon ${basenameGTF}.hisat2_exons.txt"

"""
if [ ! -e "${resultDir}/${basenameGenome}.8.ht2" ] ; then
    echo "${resultDir}/${basenameGenome}.8.ht2 Hisat2 index not found"
    rm -rf $tmpResultDir $resultDir && mkdir -p $tmpResultDir && cd $tmpResultDir 
    $extract_splice_sites
    $extract_exons
    hisat2-build ${hisat2_build_parameters} $ss $exon ${genome} ${basenameGenome}
    cd .. && mv $tmpResultDir $resultDir 
fi
"""

}

bowtie2_build_parameters = params.Check_and_Build_Module_Bowtie2_Index.bowtie2_build_parameters

process Check_and_Build_Module_Bowtie2_Index {

input:
 val genome from g227_15_genomePath_g227_6
 val gtf from g227_15_gtfPath_g227_6

output:
 val resultDir  into g227_6_genomeIndexPath_g227_18, g227_6_genomeIndexPath_g218_11

when:
(params.use_Bowtie2_Index == "yes") && (params.run_checkAndBuild == "yes")

script:

basedir  = genome.substring(0, genome.lastIndexOf('/'))
gtf_dir  = gtf.substring(0, gtf.lastIndexOf('/')) 
basename = genome.substring(genome.lastIndexOf('/')+1,genome.lastIndexOf('.'))
filename = genome.substring(genome.lastIndexOf('/')+1,genome.length())
newDirName = "Bowtie2Index"
resultDir = basedir.substring(0, basedir.lastIndexOf('/')) +"/"+ newDirName 
tmpResultDir = basedir.substring(0, basedir.lastIndexOf('/')) +"/_tmp_"+ newDirName 

"""
if [ ! -e "${resultDir}/${basename}.rev.1.bt2" ] ; then
    echo "${resultDir}/${basename}.rev.1.bt2 Bowtie2 index not found"
    rm -rf $tmpResultDir $resultDir && mkdir -p $tmpResultDir && cd $tmpResultDir
    ln -s ../main/${filename} ${filename}
    bowtie2-build ${bowtie2_build_parameters} ${filename} ${basename}
    cd .. && mv $tmpResultDir $resultDir 
fi
"""



}


process Check_and_Build_Module_STAR_Index_Check_Build {

input:
 val genome from g227_15_genomePath_g227_0
 val gtf from g227_15_gtfPath_g227_0

output:
 val resultDir  into g227_0_genomeIndexPath_g227_18, g227_0_genomeIndexPath_g217_16

when:
(params.use_STAR_Index == "yes") && (params.run_checkAndBuild == "yes")

script:
star_build_parameters = params.Check_and_Build_Module_STAR_Index_Check_Build.star_build_parameters
gtf_dir  = gtf.substring(0, gtf.lastIndexOf('/')) 
indexbasedir  = gtf_dir.substring(0, gtf_dir.lastIndexOf('/'))
genome_dir  = genome.substring(0, genome.lastIndexOf('/')) 
filename = genome.substring(genome.lastIndexOf('/')+1,genome.length())
newDirName = "STARIndex" 
resultDir = indexbasedir +"/"+ newDirName 
tmpResultDir = indexbasedir +"/_tmp_"+ newDirName
"""
if [ ! -e "${resultDir}/SA" ] ; then
    echo "STAR index not found"
    rm -rf $tmpResultDir ${resultDir} && mkdir -p $tmpResultDir && cd $tmpResultDir
    STAR --runMode genomeGenerate ${star_build_parameters} --genomeDir $tmpResultDir --genomeFastaFiles ${genome} --sjdbGTFfile ${gtf}
    cd .. && mv $tmpResultDir $resultDir && cd ${resultDir}
    ln -s ../../main/${filename} ${filename}
    ln -s ../../main/${filename}.fai ${filename}.fai
fi
"""





}

//* params.gtf =  ""  //* @input
//* params.genome =  ""  //* @input
//* params.commondb =  ""  //* @input
if (!(params.run_checkAndBuild == "yes" && params.run_Sequential_Mapping  == "yes")){
g227_15_commondb_path_g227_18.into{g227_18_commondb_path_g241_34}
} else {


process Check_and_Build_Module_Check_Sequential_Mapping_Indexes {

input:
 val commondb from g227_15_commondb_path_g227_18
 val bowtieIndex from g227_13_genomeIndexPath_g227_18
 val bowtie2Index from g227_6_genomeIndexPath_g227_18
 val starIndex from g227_0_genomeIndexPath_g227_18

output:
 val commondb  into g227_18_commondb_path_g241_34

when:
params.run_checkAndBuild == "yes" && params.run_Sequential_Mapping  == "yes"

script:
"""
"""
}
}


g227_18_commondb_path_g241_34= g227_18_commondb_path_g241_34.ifEmpty([""]) 

//* params.run_Sequential_Mapping =   "yes"   //* @dropdown @options:"yes","no" @show_settings:"Sequential_Mapping" @description:"Filters out or quantify given sequence sets."
//* params.bowtieInd_rRNA =  ""  //* @input
//* params.bowtieInd_ercc =  ""  //* @input
//* params.bowtieInd_miRNA =  ""  //* @input
//* params.bowtieInd_tRNA =  ""  //* @input
//* params.bowtieInd_piRNA =  ""  //* @input
//* params.bowtieInd_snRNA =  ""  //* @input
//* params.bowtieInd_rmsk =  ""  //* @input
//* params.bowtie_index =  ""  //* @input
//* params.bowtie2_index =  ""  //* @input
//* params.star_index =  ""  //* @input

//both bowtie and bowtie2 indexes located in same path
bowtieIndexes = [rRNA: params.bowtieInd_rRNA, 
                 ercc: params.bowtieInd_ercc,
                 miRNA: params.bowtieInd_miRNA,
                 tRNA: params.bowtieInd_tRNA,
                 piRNA: params.bowtieInd_piRNA,
                 snRNA: params.bowtieInd_snRNA,
                 rmsk: params.bowtieInd_rmsk]
                 
genomeIndexes = [bowtie: params.bowtie_index,
                 bowtie2: params.bowtie2_index,
                 STAR: params.star_index+"/genome"]


//_nucleicAcidType="dna" should be defined in the autofill section of pipeline header in case dna is used.
_select_sequence = params.Sequential_Mapping_Module_Sequential_Mapping._select_sequence
index_directory = params.Sequential_Mapping_Module_Sequential_Mapping.index_directory
name_of_the_index_file = params.Sequential_Mapping_Module_Sequential_Mapping.name_of_the_index_file
_aligner = params.Sequential_Mapping_Module_Sequential_Mapping._aligner
aligner_Parameters = params.Sequential_Mapping_Module_Sequential_Mapping.aligner_Parameters
description = params.Sequential_Mapping_Module_Sequential_Mapping.description
filter_Out = params.Sequential_Mapping_Module_Sequential_Mapping.filter_Out
sense_antisense = params.Sequential_Mapping_Module_Sequential_Mapping.sense_antisense

desc_all=[]
description.eachWithIndex() {param,i -> 
    if (param.isEmpty()){
        desc_all[i] = name_of_the_index_file[i]
    }  else {
        desc_all[i] = param.replaceAll("[ |.|;]", "_")
    }
}
custom_index=[]
index_directory.eachWithIndex() {param,i -> 
    if (_select_sequence[i] == "genome"){
        custom_index[i] = genomeIndexes[_aligner[i]]
    }else if (_select_sequence[i] == "custom"){
        custom_index[i] = param+"/"+name_of_the_index_file[i]
    }else {
        custom_index[i] = bowtieIndexes[_select_sequence[i]]
    }
}

mapList = []
paramList = []
alignerList = []
filterList = []
indexList = []
senseList = []

//concat default mapping and custom mapping
mapList = (desc_all) 
paramList = (aligner_Parameters)
alignerList = (_aligner)
filterList = (filter_Out)
indexList = (custom_index)
senseList = (sense_antisense)

mappingList = mapList.join(" ") // convert into space separated format in order to use in bash for loop
paramsList = paramList.join(",") // convert into comma separated format in order to use in as array in bash
alignersList = alignerList.join(",") 
filtersList = filterList.join(",") 
indexesList = indexList.join(",") 
senseList = senseList.join(",")

//* @style @condition:{remove_duplicates="yes",remove_duplicates_based_on_UMI_after_mapping},{remove_duplicates="no"},{_select_sequence="custom", index_directory,name_of_the_index_file,description,_aligner,aligner_Parameters,filter_Out,sense_antisense},{_select_sequence=("rRNA","ercc","miRNA","tRNA","piRNA","snRNA","rmsk","genome"),_aligner,aligner_Parameters,filter_Out,sense_antisense}  @array:{_select_sequence,_select_sequence, index_directory,name_of_the_index_file,_aligner,aligner_Parameters,filter_Out,sense_antisense,description} @multicolumn:{_select_sequence,_select_sequence,index_directory,name_of_the_index_file,_aligner,aligner_Parameters,filter_Out, sense_antisense, description},{remove_duplicates,remove_duplicates_based_on_UMI_after_mapping}

//* autofill
if ($HOSTNAME == "default"){
    $CPU  = 4
    $MEMORY = 20
}
//* platform
if ($HOSTNAME == "ghpcc06.umassrc.org"){
    $TIME = 2000
    $CPU  = 4
    $MEMORY = 20
    $QUEUE = "long"
}
//* platform
//* autofill
if (!(params.run_Sequential_Mapping == "yes")){
g213_20_reads_g241_34.into{g241_34_reads_g_127; g241_34_reads_g240_21; g241_34_reads_g242_22}
g241_34_bowfiles_g241_26 = Channel.empty()
g241_34_bowfiles_g_177 = Channel.empty()
g241_34_bam_file_g241_35 = Channel.empty()
g241_34_bam_file_g241_36 = Channel.empty()
g241_34_bam_index_g241_35 = Channel.empty()
g241_34_bam_index_g241_36 = Channel.empty()
g241_34_filter_g241_26 = Channel.empty()
g241_34_log_file_g241_30 = Channel.empty()
} else {


process Sequential_Mapping_Module_Sequential_Mapping {

publishDir params.outdir, overwrite: true, mode: 'copy',
	saveAs: {filename ->
	if (filename =~ /.*\/.*_sorted.bam$/) "sequential_mapping/$filename"
	else if (filename =~ /.*\/.*_sorted.bam.bai$/) "sequential_mapping/$filename"
	else if (filename =~ /.*\/.*_duplicates_stats.log$/) "sequential_mapping/$filename"
}

input:
 set val(name), file(reads) from g213_20_reads_g241_34
 val mate from g_229_mate_g241_34
 val commondb_path from g227_18_commondb_path_g241_34

output:
 set val(name), file("final_reads/*q")  into g241_34_reads_g_127, g241_34_reads_g240_21, g241_34_reads_g242_22
 set val(name), file("bowfiles/?*") optional true  into g241_34_bowfiles_g241_26, g241_34_bowfiles_g_177
 file "*/*_sorted.bam" optional true  into g241_34_bam_file_g241_35
 file "*/*_sorted.bam.bai" optional true  into g241_34_bam_index_g241_35
 val filtersList  into g241_34_filter_g241_26
 file "*/*_sorted.dedup.bam" optional true  into g241_34_bam_file_g241_36
 file "*/*_sorted.dedup.bam.bai" optional true  into g241_34_bam_index_g241_36
 file "*/*_duplicates_stats.log" optional true  into g241_34_log_file_g241_30

errorStrategy 'retry'
maxRetries 2

when:
params.run_Sequential_Mapping == "yes"

script:
nameAll = reads.toString()
nameArray = nameAll.split(' ')
def file2;

if (nameAll.contains('.gz')) {
    newName =  nameArray[0] - ~/(\.fastq.gz)?(\.fq.gz)?$/
    file1 =  nameArray[0] - '.gz' 
    if (mate == "pair") {file2 =  nameArray[1] - '.gz'}
    runGzip = "ls *.gz | xargs -i echo gzip -df {} | sh"
} else {
    newName =  nameArray[0] - ~/(\.fastq)?(\.fq)?$/
    file1 =  nameArray[0]
    if (mate == "pair") {file2 =  nameArray[1]}
    runGzip = ''
}

remove_duplicates = params.Sequential_Mapping_Module_Sequential_Mapping.remove_duplicates
remove_duplicates_based_on_UMI_after_mapping = params.Sequential_Mapping_Module_Sequential_Mapping.remove_duplicates_based_on_UMI_after_mapping
remove_previous_reads = params.Sequential_Mapping_Module_Sequential_Mapping.remove_previous_reads

"""
#!/bin/bash
mkdir reads final_reads bowfiles
workflowWorkDir=\$(cd ../../ && pwd)
if [ -n "${mappingList}" ]; then
    $runGzip
    #rename files to standart format
    if [ "${mate}" == "pair" ]; then
        mv $file1 ${name}.1.fastq 2>/dev/null
        mv $file2 ${name}.2.fastq 2>/dev/null
        mv ${name}.1.fastq ${name}.2.fastq reads/.
    else
        mv $file1 ${name}.fastq 2>/dev/null
        mv ${name}.fastq reads/.
    fi
    #sequential mapping
    k=0
    prev="reads"
    IFS=',' read -r -a paramsListAr <<< "${paramsList}" #create comma separated array 
    IFS=',' read -r -a filtersListAr <<< "${filtersList}"
    IFS=',' read -r -a indexesListAr <<< "${indexesList}"
    IFS=',' read -r -a alignersListAr <<< "${alignersList}"
    IFS=',' read -r -a senseListAr <<< "${senseList}"
    wrkDir=\$(pwd)
    for rna_set in ${mappingList}
    do
        ((k++))
        printf -v k2 "%02d" "\$k" #turn into two digit format
        mkdir -p \${rna_set}/unmapped
        cd \$rna_set
        ## create link of the target file to prevent "too many symlinks error"
        for r in \${wrkDir}/\${prev}/*; do
            targetRead=\$(readlink -e \$r)
            rname=\$(basename \$r)
            echo "INFO: ln -s \$targetRead \$rname"
            ln -s \$targetRead \$rname
        done
        genomeDir=`dirname "\${indexesListAr[\$k-1]}"`
        echo "INFO: genomeDir: \$genomeDir"
        if [ -e "\${indexesListAr[\$k-1]}.1.bt2" -o  -e "\${indexesListAr[\$k-1]}.fa"  -o  -e "\${indexesListAr[\$k-1]}.fasta"  -o  -e "\$genomeDir/SAindex" ]; then
            if [ -e "\${indexesListAr[\$k-1]}.fa" ] ; then
                fasta=\${indexesListAr[\$k-1]}.fa
            elif [ -e "\${indexesListAr[\$k-1]}.fasta" ] ; then
                fasta=\${indexesListAr[\$k-1]}.fasta
            fi
            echo "INFO: fasta: \$fasta"
            if [ -e "\${indexesListAr[\$k-1]}.1.bt2" -a "\${alignersListAr[\$k-1]}" == "bowtie2" ] ; then
                echo "INFO: \${indexesListAr[\$k-1]}.1.bt2 Bowtie2 index found."
            elif [ -e "\${indexesListAr[\$k-1]}.1.ebwt" -a "\${alignersListAr[\$k-1]}" == "bowtie" ] ; then
                echo "INFO: \${indexesListAr[\$k-1]}.1.ebwt Bowtie index found."
            elif [ -e "\$genomeDir/SAindex" -a "\${alignersListAr[\$k-1]}" == "STAR" ] ; then
                echo "INFO: \$genomeDir/SAindex STAR index found."
            elif [ -e "\${indexesListAr[\$k-1]}.fa" -o  -e "\${indexesListAr[\$k-1]}.fasta" ] ; then
                if [ "\${alignersListAr[\$k-1]}" == "bowtie2" ]; then
                    bowtie2-build \$fasta \${indexesListAr[\$k-1]}
                elif [ "\${alignersListAr[\$k-1]}" == "STAR" ]; then
                    if [ -e "\${indexesListAr[\$k-1]}.gtf" ]; then
                        STAR --runMode genomeGenerate --genomeDir \$genomeDir --genomeFastaFiles \$fasta --sjdbGTFfile \${indexesListAr[\$k-1]}.gtf --genomeSAindexNbases 5
                    else
                        echo "WARNING: \${indexesListAr[\$k-1]}.gtf not found. STAR index is not generated."
                    fi
                elif [ "\${alignersListAr[\$k-1]}" == "bowtie" ]; then
                    bowtie-build \$fasta \${indexesListAr[\$k-1]}
                fi
            fi
                
            if [ "${mate}" == "pair" ]; then
                if [ "\${alignersListAr[\$k-1]}" == "bowtie2" ]; then
                    bowtie2 \${paramsListAr[\$k-1]} -x \${indexesListAr[\$k-1]} --no-unal --un-conc unmapped/${name}.unmapped.fastq -1 ${name}.1.fastq -2 ${name}.2.fastq --al-conc ${name}.fq.mapped -S \${rna_set}_${name}_alignment.sam 2>&1 | tee \${k2}_${name}.bow_\${rna_set}
                elif [ "\${alignersListAr[\$k-1]}" == "STAR" ]; then
                    STAR \${paramsListAr[\$k-1]}  --genomeDir \$genomeDir --readFilesIn ${name}.1.fastq ${name}.2.fastq --outSAMtype SAM  --outFileNamePrefix ${name}.star --outReadsUnmapped Fastx
                    mv ${name}.starAligned.out.sam \${rna_set}_${name}_alignment.sam
                    mv ${name}.starUnmapped.out.mate1 unmapped/${name}.unmapped.1.fastq
                    mv ${name}.starUnmapped.out.mate2 unmapped/${name}.unmapped.2.fastq
                    mv ${name}.starLog.final.out \${k2}_${name}.star_\${rna_set}
                elif [ "\${alignersListAr[\$k-1]}" == "bowtie" ]; then
                    bowtie \${paramsListAr[\$k-1]}   \${indexesListAr[\$k-1]}  --un  unmapped/${name}.unmapped.fastq -1 ${name}.1.fastq -2 ${name}.2.fastq -S  \${rna_set}_${name}_alignment.sam 2>&1 | tee \${k2}_${name}.bow1_\${rna_set}  
                    mv unmapped/${name}.unmapped_1.fastq unmapped/${name}.unmapped.1.fastq
                    mv unmapped/${name}.unmapped_2.fastq unmapped/${name}.unmapped.2.fastq
                fi
            else
                if [ "\${alignersListAr[\$k-1]}" == "bowtie2" ]; then
                    bowtie2 \${paramsListAr[\$k-1]} -x \${indexesListAr[\$k-1]} --no-unal --un  unmapped/${name}.unmapped.fastq -U ${name}.fastq --al ${name}.fq.mapped -S \${rna_set}_${name}_alignment.sam 2>&1 | tee \${k2}_${name}.bow_\${rna_set}  
                elif [ "\${alignersListAr[\$k-1]}" == "STAR" ]; then
                    STAR \${paramsListAr[\$k-1]}  --genomeDir \$genomeDir --readFilesIn ${name}.fastq --outSAMtype SAM  --outFileNamePrefix ${name}.star --outReadsUnmapped Fastx
                    mv ${name}.starAligned.out.sam \${rna_set}_${name}_alignment.sam
                    mv ${name}.starUnmapped.out.mate1 unmapped/${name}.unmapped.fastq
                    mv ${name}.starLog.final.out \${k2}_${name}.star_\${rna_set}
                elif [ "\${alignersListAr[\$k-1]}" == "bowtie" ]; then
                    bowtie \${paramsListAr[\$k-1]}  \${indexesListAr[\$k-1]}  --un  unmapped/${name}.unmapped.fastq  ${name}.fastq  -S \${rna_set}_${name}_alignment.sam 2>&1 | tee \${k2}_${name}.bow1_\${rna_set}  
                    
                fi
            fi
            echo "INFO: samtools view -bT \${fasta} \${rna_set}_${name}_alignment.sam > \${rna_set}_${name}_alignment.bam"
            samtools view -bT \${fasta} \${rna_set}_${name}_alignment.sam > \${rna_set}_${name}_alignment.bam
            rm -f \${rna_set}_${name}_alignment.sam
            if [ "\${alignersListAr[\$k-1]}" == "bowtie" ]; then
                mv \${rna_set}_${name}_alignment.bam \${rna_set}_${name}_tmp0.bam
                echo "INFO: samtools view -F 0x04 -b \${rna_set}_${name}_tmp0.bam > \${rna_set}_${name}_alignment.bam"
                samtools view -F 0x04 -b \${rna_set}_${name}_tmp0.bam > \${rna_set}_${name}_alignment.bam  # Remove unmapped reads
                if [ "${mate}" == "pair" ]; then
                    echo "# unique mapped reads: \$(samtools view -f 0x40 -F 0x4 -q 255 \${rna_set}_${name}_alignment.bam | cut -f 1 | sort -T '.' | uniq | wc -l)" >> \${k2}_${name}.bow1_\${rna_set}
                else
                    echo "# unique mapped reads: \$(samtools view -F 0x40 -q 255 \${rna_set}_${name}_alignment.bam | cut -f 1 | sort -T '.' | uniq | wc -l)" >> \${k2}_${name}.bow1_\${rna_set}
                fi
            fi
            if [ "${mate}" == "pair" ]; then
                mv \${rna_set}_${name}_alignment.bam \${rna_set}_${name}_alignment.tmp1.bam
                echo "INFO: samtools sort -n -o \${rna_set}_${name}_alignment.tmp2 \${rna_set}_${name}_alignment.tmp1.bam"
                samtools sort -n -o \${rna_set}_${name}_alignment.tmp2.bam \${rna_set}_${name}_alignment.tmp1.bam 
                echo "INFO: samtools view -bf 0x02 \${rna_set}_${name}_alignment.tmp2.bam >\${rna_set}_${name}_alignment.bam"
                samtools view -bf 0x02 \${rna_set}_${name}_alignment.tmp2.bam >\${rna_set}_${name}_alignment.bam
                rm \${rna_set}_${name}_alignment.tmp1.bam \${rna_set}_${name}_alignment.tmp2.bam
            fi
            echo "INFO: samtools sort -o \${rna_set}@${name}_sorted.bam \${rna_set}_${name}_alignment.bam"
            samtools sort -o \${rna_set}@${name}_sorted.bam \${rna_set}_${name}_alignment.bam 
            echo "INFO: samtools index \${rna_set}@${name}_sorted.bam"
            samtools index \${rna_set}@${name}_sorted.bam
            
             # split sense and antisense bam files. 
            if [ "\${senseListAr[\$k-1]}" == "Yes" ]; then
                if [ "${mate}" == "pair" ]; then
                    echo "INFO: paired end sense antisense separation"
                	samtools view -f 65 -b \${rna_set}@${name}_sorted.bam >\${rna_set}@${name}_forward_sorted.bam
	                samtools index \${rna_set}@${name}_forward_sorted.bam
	                samtools view -F 16 -b \${rna_set}@${name}_forward_sorted.bam >\${rna_set}@${name}_sense_sorted.bam
	                samtools index \${rna_set}@${name}_sense_sorted.bam
	                samtools view -f 16 -b \${rna_set}@${name}_forward_sorted.bam >\${rna_set}@${name}_antisense_sorted.bam
	                samtools index \${rna_set}@${name}_antisense_sorted.bam
                else
	                echo "INFO: single end sense antisense separation"
	                samtools view -F 16 -b \${rna_set}@${name}_sorted.bam >\${rna_set}@${name}_sense_sorted.bam
	                samtools index \${rna_set}@${name}_sense_sorted.bam
	                samtools view -f 16 -b \${rna_set}@${name}_sorted.bam >\${rna_set}@${name}_antisense_sorted.bam
	                samtools index \${rna_set}@${name}_antisense_sorted.bam
                fi
            fi
            
            if [ "${remove_duplicates}" == "yes" ]; then
                ## check read header whether they have UMI tags which are separated with underscore.(eg. NS5HGY:2:11_GTATAACCTT)
                umiCheck=\$(samtools view \${rna_set}@${name}_sorted.bam |head -n 1 | awk 'BEGIN {FS="\\t"}; {print \$1}' | awk 'BEGIN {FS=":"}; \$NF ~ /_/ {print \$NF}')
                
                # based on remove_duplicates_based_on_UMI_after_mapping
                if [ "${remove_duplicates_based_on_UMI_after_mapping}" == "yes" -a ! -z "\$umiCheck" ]; then
                    echo "INFO: umi_mark_duplicates.py will be executed for removing duplicates from bam file"
                    echo "python umi_mark_duplicates.py -f \${rna_set}@${name}_sorted.bam -p 4"
                    python umi_mark_duplicates.py -f \${rna_set}@${name}_sorted.bam -p 4
                else
                    echo "INFO: Picard MarkDuplicates will be executed for removing duplicates from bam file"
                    if [ "${remove_duplicates_based_on_UMI_after_mapping}" == "yes"  ]; then
                        echo "WARNING: Read header have no UMI tags which are separated with underscore. Picard MarkDuplicates will be executed to remove duplicates from alignment file (bam) instead of remove_duplicates_based_on_UMI_after_mapping."
                    fi
                    echo "INFO: picard MarkDuplicates OUTPUT=\${rna_set}@${name}_sorted.deumi.sorted.bam METRICS_FILE=${name}_picard_PCR_duplicates.log  VALIDATION_STRINGENCY=LENIENT REMOVE_DUPLICATES=false INPUT=\${rna_set}@${name}_sorted.bam"
                    picard MarkDuplicates OUTPUT=\${rna_set}@${name}_sorted.deumi.sorted.bam METRICS_FILE=${name}_picard_PCR_duplicates.log  VALIDATION_STRINGENCY=LENIENT REMOVE_DUPLICATES=false INPUT=\${rna_set}@${name}_sorted.bam 
                fi
                #get duplicates stats (read the sam flags)
                samtools flagstat \${rna_set}@${name}_sorted.deumi.sorted.bam > \${k2}@\${rna_set}@${name}_duplicates_stats.log
                #remove alignments marked as duplicates
                samtools view -b -F 0x400 \${rna_set}@${name}_sorted.deumi.sorted.bam > \${rna_set}@${name}_sorted.deumi.sorted.bam.x_dup
                #sort deduplicated files by chrom pos
                echo "INFO: samtools sort -o \${rna_set}@${name}_sorted.dedup.bam \${rna_set}@${name}_sorted.deumi.sorted.bam.x_dup"
                samtools sort -o \${rna_set}@${name}_sorted.dedup.bam \${rna_set}@${name}_sorted.deumi.sorted.bam.x_dup 
                samtools index \${rna_set}@${name}_sorted.dedup.bam
                #get flagstat after dedup
                echo "##After Deduplication##" >> \${k2}@\${rna_set}@${name}_duplicates_stats.log
                samtools flagstat \${rna_set}@${name}_sorted.dedup.bam >> \${k2}@\${rna_set}@${name}_duplicates_stats.log
                # split sense and antisense bam files. 
	            if [ "\${senseListAr[\$k-1]}" == "Yes" ]; then
	                if [ "${mate}" == "pair" ]; then
	                    echo "INFO: paired end sense antisense separation"
	                	samtools view -f 65 -b \${rna_set}@${name}_sorted.dedup.bam >\${rna_set}@${name}_forward_sorted.dedup.bam
		                samtools index \${rna_set}@${name}_forward_sorted.dedup.bam
		                samtools view -F 16 -b \${rna_set}@${name}_forward_sorted.dedup.bam>\${rna_set}@${name}_sense_sorted.dedup.bam
		                samtools index \${rna_set}@${name}_sense_sorted.dedup.bam
		                samtools view -f 16 -b \${rna_set}@${name}_forward_sorted.dedup.bam >\${rna_set}@${name}_antisense_sorted.dedup.bam
		                samtools index \${rna_set}@${name}_antisense_sorted.dedup.bam
	                else
		                echo "INFO: single end sense antisense separation"
		                samtools view -F 16 -b \${rna_set}@${name}_sorted.dedup.bam >\${rna_set}@${name}_sense_sorted.dedup.bam
		                samtools index \${rna_set}@${name}_sense_sorted.dedup.bam
		                samtools view -f 16 -b \${rna_set}@${name}_sorted.dedup.bam >\${rna_set}@${name}_antisense_sorted.dedup.bam
		                samtools index \${rna_set}@${name}_antisense_sorted.dedup.bam
	                fi
	            fi
            fi
            
        
            for file in unmapped/*; do mv \$file \${file/.unmapped/}; done ##remove .unmapped from filename
            if [ "\${alignersListAr[\$k-1]}" == "bowtie2" ]; then
                grep -v Warning \${k2}_${name}.bow_\${rna_set} > ${name}.tmp
                mv ${name}.tmp \${k2}_${name}.bow_\${rna_set}
                cp \${k2}_${name}.bow_\${rna_set} ./../bowfiles/.
            elif [ "\${alignersListAr[\$k-1]}" == "bowtie" ]; then
                cp \${k2}_${name}.bow1_\${rna_set} ./../bowfiles/.
            elif [ "\${alignersListAr[\$k-1]}" == "STAR" ]; then
                cp \${k2}_${name}.star_\${rna_set} ./../bowfiles/.
            fi
            cd ..
            # if filter is on, remove previously created unmapped fastq. 
            if [ "\${filtersListAr[\$k-1]}" == "Yes" ]; then
                if [ "\${prev}" != "reads" ]; then
                    echo "INFO: remove prev: \${prev}/*"
                    rm -rf \${prev}/*
                elif  [ "${remove_previous_reads}" == "true" ]; then
                    echo "INFO: inputs reads will be removed if they are located in the workdir"
                    for f in \${prev}/*; do
                        targetFile=\$(readlink -e \$f)
                        echo "INFO: targetFile: \$targetFile"
                        if [[ \$targetFile == *"\${workflowWorkDir}"* ]]; then
                            rm -f \$targetFile
                            echo "INFO: \$targetFile located in workdir and deleted."
                        fi
                    done
                fi
            # if filter is off remove current unmapped fastq
            else
                echo "INFO: remove \${rna_set}/unmapped/*"
                rm -rf \${rna_set}/unmapped/*
            fi
        else
            echo "WARNING: \${indexesListAr[\$k-1]} Mapping skipped. File not found."
            cd unmapped 
            ln -s \${wrkDir}/\${rna_set}/*fastq .
            cd ..
            cd ..
        fi
        
        if [ "\${filtersListAr[\$k-1]}" == "Yes" ]; then
            prev=\${rna_set}/unmapped
        fi
    done
    cd final_reads && ln -s \${wrkDir}/\${prev}/* .
else 
    mv ${reads} final_reads/.
fi
"""

}
}


//* params.run_Split_Fastq =  "no"  //* @dropdown @options:"yes","no" @show_settings:"SplitFastq" @description:"Splits Fastq files before aligning with Star, Hisat2 or Tophat2 to speed up the process. However, it will require more disk space."
readsPerFile = params.SplitFastq.readsPerFile
//Since splitFastq operator requires flat file structure, first convert grouped structure to flat, execute splitFastq, and then return back to original grouped structure
//.map(flatPairsClosure).splitFastq(splitFastqParams).map(groupPairsClosure)

//Mapping grouped read structure to flat structure
flatPairsClosure = {row -> if(row[1] instanceof Collection) {
        if (row[1][1]){
            tuple(row[0], file(row[1][0]), file(row[1][1]))
        } else {
            tuple(row[0], file(row[1][0]))
        }
    } else {
        tuple(row[0], file(row[1]))
    }
}

//Mapping flat read structure to grouped read structure
groupPairsClosure = {row -> tuple(row[0], (row[2]) ? [file(row[1]), file(row[2])] : [file(row[1])])}

// if mate of split process different than rest of the pipeline, use "mate_split" as input parameter. Otherwise use default "mate" as input parameter
mateParamName = (params.mate_split) ? "mate_split" : "mate"
splitFastqParams = ""
if (params[mateParamName] != "pair"){
    splitFastqParams = [by: readsPerFile, file:true]
}else {
    splitFastqParams = [by: readsPerFile, pe:true, file:true]
}

//* autofill
//* platform
if ($HOSTNAME == "ghpcc06.umassrc.org"){
    $TIME = 240
    $CPU  = 1
    $MEMORY = 8
    $QUEUE = "short"
}
//* platform
//* autofill
if (!(params.run_Split_Fastq == "yes")){
g241_34_reads_g_127.into{g_127_reads_g216_11; g_127_reads_g217_16; g_127_reads_g218_11}
} else {


process SplitFastq {

input:
 set val(name), file(reads) from g241_34_reads_g_127.map(flatPairsClosure).splitFastq(splitFastqParams).map(groupPairsClosure)

output:
 set val(name), file("split/*q")  into g_127_reads_g216_11, g_127_reads_g217_16, g_127_reads_g218_11

errorStrategy 'retry'
maxRetries 3

when:
params.run_Split_Fastq == "yes"

script:
"""    
mkdir -p split
mv ${reads} split/.
"""
}
}


g227_6_genomeIndexPath_g218_11= g227_6_genomeIndexPath_g218_11.ifEmpty([""]) 

//* params.bowtie2_index =  ""  //* @input
//* params.gtf =  ""  //* @input


//* autofill
if ($HOSTNAME == "default"){
    $CPU  = 3
    $MEMORY = 24
}
//* platform
if ($HOSTNAME == "ghpcc06.umassrc.org"){
    $TIME = 2500
    $CPU  = 4
    $MEMORY = 24
    $QUEUE = "long"
}
//* platform
//* autofill

process Tophat2_Module_Map_Tophat2 {

input:
 val mate from g_229_mate_g218_11
 set val(name), file(reads) from g_127_reads_g218_11
 val tophat2_index from g227_6_genomeIndexPath_g218_11

output:
 set val(name), file("${newName}.bam")  into g218_11_mapped_reads_g218_13
 set val(name), file("${newName}_unmapped.bam")  into g218_11_unmapped_reads
 set val(name), file("${newName}_align_summary.txt")  into g218_11_summary_g218_3, g218_11_summary_g_177

errorStrategy 'retry'

when:
(params.run_Tophat && (params.run_Tophat == "yes")) || !params.run_Tophat

script:
tophat2_parameters = params.Tophat2_Module_Map_Tophat2.tophat2_parameters
nameAll = reads.toString()
nameArray = nameAll.split(' ')

if (nameAll.contains('.gz')) {
    newName =  nameArray[0] - ~/(\.fastq.gz)?(\.fq.gz)?$/
    file =  nameAll - '.gz' - '.gz'
    runGzip = "ls *.gz | xargs -i echo gzip -df {} | sh"
} else {
    newName =  nameArray[0] - ~/(\.fastq)?(\.fq)?$/
    file =  nameAll 
    runGzip = ''
}

"""
$runGzip
if [ "${mate}" == "pair" ]; then
    tophat2 ${tophat2_parameters}  --keep-tmp -G ${params.gtf} -o . ${params.bowtie2_index} $file
else
    tophat2 ${tophat2_parameters}  --keep-tmp -G ${params.gtf} -o . ${params.bowtie2_index} $file
fi

if [ -f unmapped.bam ]; then
    mv unmapped.bam ${newName}_unmapped.bam
else
    touch ${newName}_unmapped.bam
fi

mv accepted_hits.bam ${newName}.bam
mv align_summary.txt ${newName}_align_summary.txt
"""

}


//* autofill
//* platform
if ($HOSTNAME == "ghpcc06.umassrc.org"){
    $TIME = 2000
    $CPU  = 1
    $MEMORY = 8
    $QUEUE = "long"
}
//* platform
//* autofill

process Tophat2_Module_Merge_Bam {

publishDir params.outdir, overwrite: true, mode: 'copy',
	saveAs: {filename ->
	if (filename =~ /.*_sorted.*bai$/) "tophat2/$filename"
	else if (filename =~ /.*_sorted.*bam$/) "tophat2/$filename"
}

input:
 set val(oldname), file(bamfiles) from g218_11_mapped_reads_g218_13.groupTuple()

output:
 set val(oldname), file("${oldname}.bam")  into g218_13_merged_bams
 set val(oldname), file("*_sorted*bai")  into g218_13_bam_index
 set val(oldname), file("*_sorted*bam")  into g218_13_sorted_bam_g222_121, g218_13_sorted_bam_g222_122, g218_13_sorted_bam_g222_123, g218_13_sorted_bam_g222_124, g218_13_sorted_bam_g222_126

errorStrategy 'retry'
maxRetries 2

shell:
'''
num=$(echo "!{bamfiles.join(" ")}" | awk -F" " '{print NF-1}')
if [ "${num}" -gt 0 ]; then
    samtools merge !{oldname}.bam !{bamfiles.join(" ")} && samtools sort -o !{oldname}_sorted.bam !{oldname}.bam && samtools index !{oldname}_sorted.bam
else
    mv !{bamfiles.join(" ")} !{oldname}.bam 2>/dev/null || true
    samtools sort  -o !{oldname}_sorted.bam !{oldname}.bam && samtools index !{oldname}_sorted.bam
fi
'''
}

//* params.gtf =  ""  //* @input


process BAM_Analysis_Tophat2_featureCounts {

publishDir params.outdir, overwrite: true, mode: 'copy',
	saveAs: {filename ->
	if (filename =~ /.*$/) "featureCounts_after_Tophat2/$filename"
}

input:
 set val(name), file(bam) from g218_13_sorted_bam_g222_126
 val paired from g_229_mate_g222_126
 each run_params from g222_125_run_parameters_g222_126

output:
 file "*"  into g222_126_outputFileTSV_g222_117

script:
pairText = ""
if (paired == "pair"){
    pairText = "-p"
}

run_name = run_params["run_name"] 
run_parameters = run_params["run_parameters"] 

"""
featureCounts ${pairText} ${run_parameters} -a ${params.gtf} -o ${name}@${run_name}@fCounts.txt ${bam}
## remove first line
sed -i '1d' ${name}@${run_name}@fCounts.txt

"""
}

//* autofill
//* platform
if ($HOSTNAME == "ghpcc06.umassrc.org"){
    $TIME = 30
    $CPU  = 1
    $MEMORY = 10
    $QUEUE = "short"
}
//* platform
//* autofill

process BAM_Analysis_Tophat2_summary_featureCounts {

publishDir params.outdir, overwrite: true, mode: 'copy',
	saveAs: {filename ->
	if (filename =~ /.*_featureCounts.tsv$/) "featureCounts_after_Tophat2_summary/$filename"
	else if (filename =~ /.*_featureCounts.sum.tsv$/) "featureCounts_after_Tophat2_details/$filename"
}

input:
 file featureCountsOut from g222_126_outputFileTSV_g222_117.collect()

output:
 file "*_featureCounts.tsv"  into g222_117_outputFile
 file "*_featureCounts.sum.tsv"  into g222_117_outFileTSV

shell:
'''
#!/usr/bin/env perl

# Step 1: Merge count files
my %tf = ( expected_count => 6 );
my @run_name=();
chomp(my $contents = `ls *@fCounts.txt`);
my @files = split(/[\\n]+/, $contents);
foreach my $file (@files){
    $file=~/(.*)\\@(.*)\\@fCounts\\.txt/;
    my $runname = $2;
    push(@run_name, $runname) unless grep{$_ eq $runname} @run_name;
}


my @expectedCount_ar = ("expected_count");
for($l = 0; $l <= $#run_name; $l++) {
    my $runName = $run_name[$l];
    for($ll = 0; $ll <= $#expectedCount_ar; $ll++) {
        my $expectedCount = $expectedCount_ar[$ll];
    
        my @a=();
        my %b=();
        my %c=();
        my $i=0;
        chomp(my $contents = `ls *\\@${runName}\\@fCounts.txt`);
        my @files = split(/[\\n]+/, $contents);
        foreach my $file (@files){
        $i++;
        $file=~/(.*)\\@${runName}\\@fCounts\\.txt/;
        my $libname = $1; 
        $a[$i]=$libname;
        open IN, $file;
            $_=<IN>;
            while(<IN>){
                my @v=split; 
                $b{$v[0]}{$i}=$v[$tf{$expectedCount}];
                $c{$v[0]}=$v[5]; #length column
            }
            close IN;
        }
        my $outfile="$runName"."_featureCounts.tsv";
        open OUT, ">$outfile";
        if ($runName eq "transcript_id") {
            print OUT "transcript\tlength";
        } else {
            print OUT "gene\tlength";
        }
    
        for(my $j=1;$j<=$i;$j++) {
            print OUT "\t$a[$j]";
        }
        print OUT "\n";
    
        foreach my $key (keys %b) {
            print OUT "$key\t$c{$key}";
            for(my $j=1;$j<=$i;$j++){
                print OUT "\t$b{$key}{$j}";
            }
            print OUT "\n";
        }
        close OUT;
         
    }
}


	

# Step 2: Merge summary files
for($l = 0; $l <= $#run_name; $l++) {
    my $runName = $run_name[$l];
    my @a=();
    my %b=();
    my $i=0;
    chomp(my $contents = `ls *\\@${runName}\\@fCounts.txt.summary`);
    my @files = split(/[\\n]+/, $contents);
    foreach my $file (@files){
        $i++;
        $file=~/(.*)\\@${runName}\\@fCounts\\.txt\\.summary/;
        my $libname = $1; 
        $a[$i]=$libname;
        open IN, $file;
        $_=<IN>;
        while(<IN>){
            my @v=split; 
            $b{$v[0]}{$i}=$v[1];
        }
        close IN;
    }
    my $outfile="$runName"."_featureCounts.sum.tsv";
    open OUT, ">$outfile";
    print OUT "criteria";
    for(my $j=1;$j<=$i;$j++) {
        print OUT "\t$a[$j]";
    }
    print OUT "\n";
    
    foreach my $key (keys %b) {
        print OUT "$key";
        for(my $j=1;$j<=$i;$j++){
            print OUT "\t$b{$key}{$j}";
        }
        print OUT "\n";
    }
    close OUT;
}

'''
}

//* params.genome_sizes =  ""  //* @input

//* autofill
if ($HOSTNAME == "default"){
    $CPU  = 1
    $MEMORY = 48
}
//* platform
if ($HOSTNAME == "ghpcc06.umassrc.org"){
    $TIME = 240
    $CPU  = 1
    $MEMORY = 48
    $QUEUE = "short"
} 
//* platform
//* autofill

process BAM_Analysis_Tophat2_UCSC_BAM2BigWig_converter {

publishDir params.outdir, overwrite: true, mode: 'copy',
	saveAs: {filename ->
	if (filename =~ /.*.bw$/) "bigwig_tophat/$filename"
}

input:
 set val(name), file(bam) from g218_13_sorted_bam_g222_124

output:
 file "*.bw"  into g222_124_outputFileBw

when:
(params.run_BigWig_Conversion && (params.run_BigWig_Conversion == "yes")) || !params.run_BigWig_Conversion

script:
nameAll = bam.toString()
if (nameAll.contains('_sorted.bam')) {
    runSamtools = "samtools index ${nameAll}"
    nameFinal = nameAll
} else {
    runSamtools = "samtools sort -o ${name}_sorted.bam $bam && samtools index ${name}_sorted.bam "
    nameFinal = "${name}_sorted.bam"
}

"""
$runSamtools
bedtools genomecov -split -bg -ibam ${nameFinal} -g ${params.genome_sizes} > ${name}.bg 
wigToBigWig -clip -itemsPerSlot=1 ${name}.bg ${params.genome_sizes} ${name}.bw 
"""
}

igv_extention_factor = params.BAM_Analysis_Tophat2_IGV_BAM2TDF_converter.igv_extention_factor
igv_window_size = params.BAM_Analysis_Tophat2_IGV_BAM2TDF_converter.igv_window_size

//* params.genome =  ""  //* @input

//* autofill
if ($HOSTNAME == "default"){
    $CPU  = 1
    $MEMORY = 32
}
//* platform
if ($HOSTNAME == "ghpcc06.umassrc.org"){
    $TIME = 400
    $CPU  = 1
    $MEMORY = 32
    $QUEUE = "long"
} 
//* platform
//* autofill

process BAM_Analysis_Tophat2_IGV_BAM2TDF_converter {

publishDir params.outdir, overwrite: true, mode: 'copy',
	saveAs: {filename ->
	if (filename =~ /.*.tdf$/) "igvtools_tophat2/$filename"
}

input:
 val mate from g_229_mate_g222_123
 set val(name), file(bam) from g218_13_sorted_bam_g222_123

output:
 file "*.tdf"  into g222_123_outputFileOut

when:
(params.run_IGV_TDF_Conversion && (params.run_IGV_TDF_Conversion == "yes")) || !params.run_IGV_TDF_Conversion

script:
pairedText = (params.nucleicAcidType == "dna" && mate == "pair") ? " --pairs " : ""
nameAll = bam.toString()
if (nameAll.contains('_sorted.bam')) {
    runSamtools = "samtools index ${nameAll}"
    nameFinal = nameAll
} else {
    runSamtools = "samtools sort -o ${name}_sorted.bam $bam && samtools index ${name}_sorted.bam "
    nameFinal = "${name}_sorted.bam"
}
"""
$runSamtools
igvtools count -w ${igv_window_size} -e ${igv_extention_factor} ${pairedText} ${nameFinal} ${name}.tdf ${params.genome_sizes}
"""
}

//* params.bed =  ""  //* @input

process BAM_Analysis_Tophat2_RSeQC {

publishDir params.outdir, overwrite: true, mode: 'copy',
	saveAs: {filename ->
	if (filename =~ /result\/.*.out$/) "rseqc_tophat2/$filename"
}

input:
 set val(name), file(bam) from g218_13_sorted_bam_g222_122

output:
 file "result/*.out"  into g222_122_outputFileOut_g222_95, g222_122_outputFileOut_g_177

when:
(params.run_RSeQC && (params.run_RSeQC == "yes")) || !params.run_RSeQC

script:
"""
mkdir result
read_distribution.py  -i ${bam} -r ${params.bed}> result/RSeQC.${name}.out
"""
}


process BAM_Analysis_Tophat2_RSeQC_Summary {

publishDir params.outdir, overwrite: true, mode: 'copy',
	saveAs: {filename ->
	if (filename =~ /.*.tsv$/) "rseqc_summary_tophat2/$filename"
}

input:
 file rseqcOut from g222_122_outputFileOut_g222_95.collect()
 val mate from g_229_mate_g222_95

output:
 file "*.tsv"  into g222_95_outputFileTSV

shell:
'''
#!/usr/bin/env perl
use List::Util qw[min max];
use strict;
use File::Basename;
use Getopt::Long;
use Pod::Usage; 
use Data::Dumper;

my $indir = $ENV{'PWD'};
my $outd = $ENV{'PWD'};
my @files = ();
my @outtypes = ("RSeQC");
my @order=( "Total Reads", "Total Tags" , "Total Assigned Tags", "CDS_Exons", "5'UTR_Exons", "3'UTR_Exons", "Introns", "TSS_up_1kb", "TSS_up_5kb", "TSS_up_10kb", "TES_down_1kb", "TES_down_5kb", "TES_down_10kb");
my %lines=(
  "Total Reads" => 1,
  "Total Tags" => 1,
  "Total Assigned Tags" => 1,
  "CDS_Exons" => 2,
  "5'UTR_Exons" => 2,
  "3'UTR_Exons" => 2,
  "Introns" => 2,
  "TSS_up_1kb" => 2,
  "TSS_up_5kb" => 2,
  "TSS_up_10kb" => 2,
  "TES_down_1kb" => 2,
  "TES_down_5kb" => 2,
  "TES_down_10kb" => 2
);


foreach my $outtype (@outtypes)
{

my $ext=".out";
@files = <$indir/$outtype*$ext>;

my @rowheaders=();
my @libs=();
my %vals=();
my %normvals=();
my $type = "rsem";

foreach my $d (@files){
  my $libname=basename($d, $ext);
  $libname=~s/RSeQC.//g;
  $libname=~s/rsem.out.//g;
  $libname=~s/.genome//g;
  print $libname."\\n";
  push(@libs, $libname); 
  getVals($d, $libname, \\%vals, \\%normvals, \\%lines);
}
#print Dumper(%vals);
#print Dumper(%normvals);

my $sizemetrics = keys %vals;
write_results("$outd/$outtype.$type.counts.tsv", \\@libs,\\%vals, \\@order, "region") if ($sizemetrics>0);
write_results("$outd/$outtype.$type.tagskb.tsv", \\@libs,\\%normvals, \\@order, "region") if ($sizemetrics>0);

}

sub write_results
{
  my ($outfile, $libs, $vals, $order, $name )=@_;
  open(OUT, ">$outfile");
  print OUT "$name\\t".join("\\t", @{$libs})."\\n";

  my $lib=${$libs}[0];
  foreach my $key ( @order )
  {
    if (exists ${$vals}{$lib}{$key}) {
    print OUT $key;
    foreach my $lib (@{$libs})
    {
      print OUT "\\t".${$vals}{$lib}{$key};
    } 
    print OUT "\\n";
    }
  }
  close(OUT);
}

sub getVals{
  my ($filename, $libname, $vals, $normvals, $lines)=@_;
  if (-e $filename){
     open(IN, $filename);
     while(my $line=<IN>)
     {
       chomp($line);
       my @vals_arr=split(/\\s{2,}/,$line);
       if (exists ${$lines}{$vals_arr[0]}) {
         my $idx=${$lines}{$vals_arr[0]};
         ${$vals}{$libname}{$vals_arr[0]}=$vals_arr[$idx] if (exists $vals_arr[$idx]);
         if ($idx==2) {
             ${$normvals}{$libname}{$vals_arr[0]}=$vals_arr[3] if (exists $vals_arr[3]);
         }
       }
     } 
  }
  
}
'''

}

//* params.pdfbox_path =  ""  //* @input
//* autofill
if ($HOSTNAME == "default"){
    $CPU  = 1
    $MEMORY = 32
}
//* platform
if ($HOSTNAME == "ghpcc06.umassrc.org"){
    $TIME = 240
    $CPU  = 1
    $MEMORY = 32
    $QUEUE = "short"
}
//* platform
//* autofill

process BAM_Analysis_Tophat2_Picard {

input:
 set val(name), file(bam) from g218_13_sorted_bam_g222_121

output:
 file "*_metrics"  into g222_121_outputFileOut_g222_82
 file "results/*.pdf"  into g222_121_outputFilePdf_g222_82

when:
(params.run_Picard_CollectMultipleMetrics && (params.run_Picard_CollectMultipleMetrics == "yes")) || !params.run_Picard_CollectMultipleMetrics

script:
"""
picard CollectMultipleMetrics OUTPUT=${name}_multiple.out VALIDATION_STRINGENCY=LENIENT INPUT=${bam}
mkdir results && java -jar ${params.pdfbox_path} PDFMerger *.pdf results/${name}_multi_metrics.pdf
"""
}


process BAM_Analysis_Tophat2_Picard_Summary {

publishDir params.outdir, overwrite: true, mode: 'copy',
	saveAs: {filename ->
	if (filename =~ /.*.tsv$/) "picard_summary_tophat2/$filename"
	else if (filename =~ /results\/.*.pdf$/) "picard_summary_pdf_tophat2/$filename"
}

input:
 file picardOut from g222_121_outputFileOut_g222_82.collect()
 val mate from g_229_mate_g222_82
 file picardPdf from g222_121_outputFilePdf_g222_82.collect()

output:
 file "*.tsv"  into g222_82_outputFileTSV
 file "results/*.pdf"  into g222_82_outputFilePdf

shell:
'''
#!/usr/bin/env perl
use List::Util qw[min max];
use strict;
use File::Basename;
use Getopt::Long;
use Pod::Usage; 
use Data::Dumper;

system("mkdir results && mv *.pdf results/. ");

my $indir = $ENV{'PWD'};
my $outd = $ENV{'PWD'};
my @files = ();
my @outtypes = ("CollectRnaSeqMetrics", "alignment_summary_metrics", "base_distribution_by_cycle_metrics", "insert_size_metrics", "quality_by_cycle_metrics", "quality_distribution_metrics" );

foreach my $outtype (@outtypes)
{
my $ext="_multiple.out";
$ext.=".$outtype" if ($outtype ne "CollectRnaSeqMetrics");
@files = <$indir/*$ext>;

my @rowheaders=();
my @libs=();
my %metricvals=();
my %histvals=();

my $pdffile="";
my $libname="";
foreach my $d (@files){
  my $libname=basename($d, $ext);
  print $libname."\\n";
  push(@libs, $libname); 
  getMetricVals($d, $libname, \\%metricvals, \\%histvals, \\@rowheaders);
}

my $sizemetrics = keys %metricvals;
write_results("$outd/$outtype.stats.tsv", \\@libs,\\%metricvals, \\@rowheaders, "metric") if ($sizemetrics>0);
my $sizehist = keys %histvals;
write_results("$outd/$outtype.hist.tsv", \\@libs,\\%histvals, "none", "nt") if ($sizehist>0);

}

sub write_results
{
  my ($outfile, $libs, $vals, $rowheaders, $name )=@_;
  open(OUT, ">$outfile");
  print OUT "$name\\t".join("\\t", @{$libs})."\\n";
  my $size=0;
  $size=scalar(@{${$vals}{${$libs}[0]}}) if(exists ${$libs}[0] and exists ${$vals}{${$libs}[0]} );
  
  for (my $i=0; $i<$size;$i++)
  { 
    my $rowname=$i;
    $rowname = ${$rowheaders}[$i] if ($name=~/metric/);
    print OUT $rowname;
    foreach my $lib (@{$libs})
    {
      print OUT "\\t".${${$vals}{$lib}}[$i];
    } 
    print OUT "\\n";
  }
  close(OUT);
}

sub getMetricVals{
  my ($filename, $libname, $metricvals, $histvals,$rowheaders)=@_;
  if (-e $filename){
     my $nextisheader=0;
     my $nextisvals=0;
     my $nexthist=0;
     open(IN, $filename);
     while(my $line=<IN>)
     {
       chomp($line);
       @{$rowheaders}=split(/\\t/, $line) if ($nextisheader && !scalar(@{$rowheaders})); 
       if ($nextisvals) {
         @{${$metricvals}{$libname}}=split(/\\t/, $line);
         $nextisvals=0;
       }
       if($nexthist){
          my @vals=split(/[\\s\\t]+/,$line); 
          push(@{${$histvals}{$libname}}, $vals[1]) if (exists $vals[1]);
       }
       $nextisvals=1 if ($nextisheader); $nextisheader=0;
       $nextisheader=1 if ($line=~/METRICS CLASS/);
       $nexthist=1 if ($line=~/normalized_position/);
     } 
  }
  
}
'''

}

//* autofill
//* platform
if ($HOSTNAME == "ghpcc06.umassrc.org"){
    $TIME = 240
    $CPU  = 1
    $MEMORY = 8
    $QUEUE = "short"
}
//* platform
//* autofill

process Tophat2_Module_Merge_Tophat_Summary {

input:
 set val(name), file(alignSum) from g218_11_summary_g218_3.groupTuple()
 val mate from g_229_mate_g218_3

output:
 set val(name), file("${name}_tophat_sum.tsv")  into g218_3_report_g218_9
 val "tophat2_alignment_sum"  into g218_3_name_g218_9

errorStrategy 'retry'
maxRetries 3

shell:
'''
#!/usr/bin/env perl
use List::Util qw[min max];
use strict;
use File::Basename;
use Getopt::Long;
use Pod::Usage; 
use Data::Dumper;

my %tsv;
my @headers = ();
my $name = "!{name}";

alteredAligned();

my @keys = keys %tsv;
my $summary = "$name"."_tophat_sum.tsv";
my $header_string = join("\\t", @headers);
`echo "$header_string" > $summary`;
foreach my $key (@keys){
	my $values = join("\\t", @{ $tsv{$key} });
	`echo "$values" >> $summary`;
}


sub alteredAligned
{
	my @files = qw(!{alignSum});
	my $multimappedSum;
	my $alignedSum;
	my $inputCountSum;
	push(@headers, "Sample");
    push(@headers, "Total Reads");
	push(@headers, "Multimapped Reads Aligned (Tophat2)");
	push(@headers, "Unique Reads Aligned (Tophat2)");
	foreach my $file (@files){
		my $multimapped;
		my $aligned;
		my $inputCount;
		chomp($aligned = `cat $file | grep 'Aligned pairs:' | awk '{sum=\\$3} END {print sum}'`);
		if ($aligned eq "") { # then it is single-end
		        chomp($inputCount = `cat $file | grep 'Input' | awk '{sum=\\$3} END {print sum}'`);
				chomp($aligned = `cat $file | grep 'Mapped' | awk '{sum=\\$3} END {print sum}'`);
				chomp($multimapped = `cat $file | grep 'multiple alignments' | awk '{sum+=\\$3} END {print sum}'`);
			}else{ # continue to pair end
			    chomp($inputCount = `cat $file | grep 'Input' | awk '{sum=\\$3} END {print sum}'`);
				chomp($multimapped = `cat $file | grep -A 1 'Aligned pairs:' | awk 'NR % 3 == 2 {sum+=\\$3} END {print sum}'`);
			}
        $multimappedSum += int($multimapped);
        $alignedSum += (int($aligned) - int($multimapped));
        $inputCountSum += int($inputCount);
        if ($alignedSum < 0){
            $alignedSum = 0;
        }
	}
	$tsv{$name} = [$name, $inputCountSum];
	push(@{$tsv{$name}}, $multimappedSum);
	push(@{$tsv{$name}}, $alignedSum);
}
'''

}


process Tophat2_Module_Merge_TSV_Files {

input:
 file tsv from g218_3_report_g218_9.collect()
 val outputFileName from g218_3_name_g218_9.collect()

output:
 file "${name}.tsv"  into g218_9_outputFileTSV_g_198

errorStrategy 'retry'
maxRetries 3

script:
name = outputFileName[0]
"""    
awk 'FNR==1 && NR!=1 {  getline; } 1 {print} ' *.tsv > ${name}.tsv
"""
}

g227_0_genomeIndexPath_g217_16= g227_0_genomeIndexPath_g217_16.ifEmpty([""]) 

//* params.star_index =  ""  //* @input

//* autofill
if ($HOSTNAME == "default"){
    $CPU  = 3
    $MEMORY = 32
}
//* platform
if ($HOSTNAME == "ghpcc06.umassrc.org"){
    $TIME = 1000
    $CPU  = 3
    $MEMORY = 32
    $QUEUE = "long"
} 
//* platform
//* autofill

process STAR_Module_Map_STAR {

input:
 val mate from g_229_mate_g217_16
 set val(name), file(reads) from g_127_reads_g217_16
 val star_index from g227_0_genomeIndexPath_g217_16

output:
 set val(name), file("${newName}Log.final.out")  into g217_16_outputFileOut_g217_18
 set val(name), file("${newName}.flagstat.txt")  into g217_16_outputFileTxt
 set val(name), file("${newName}Log.out")  into g217_16_logOut_g217_18
 set val(name), file("${newName}.bam")  into g217_16_mapped_reads_g217_14
 set val(name), file("${newName}SJ.out.tab")  into g217_16_outputFileTab_g217_18
 set val(name), file("${newName}Log.progress.out")  into g217_16_progressOut_g217_18
 set val(name), file("${newName}Aligned.toTranscriptome.out.bam") optional true  into g217_16_transcriptome_bam_g217_15

errorStrategy 'retry'

when:
(params.run_STAR && (params.run_STAR == "yes")) || !params.run_STAR

script:
params_STAR = params.STAR_Module_Map_STAR.params_STAR
nameAll = reads.toString()
nameArray = nameAll.split(' ')

if (nameAll.contains('.gz')) {
    newName =  nameArray[0] - ~/(\.fastq.gz)?(\.fq.gz)?$/
    file =  nameAll - '.gz' - '.gz'
    runGzip = "ls *.gz | xargs -i echo gzip -df {} | sh"
} else {
    newName =  nameArray[0] - ~/(\.fastq)?(\.fq)?$/
    file =  nameAll 
    runGzip = ''
}

"""
$runGzip
STAR ${params_STAR}  --genomeDir ${params.star_index} --readFilesIn $file --outFileNamePrefix ${newName}
echo "Alignment completed."
if [ ! -e "${newName}Aligned.toTranscriptome.out.bam" -a -e "${newName}Aligned.toTranscriptome.out.sam" ] ; then
    samtools view -S -b ${newName}Aligned.toTranscriptome.out.sam > ${newName}Aligned.toTranscriptome.out.bam
elif [ ! -e "${newName}Aligned.out.bam" -a -e "${newName}Aligned.out.sam" ] ; then
    samtools view -S -b ${newName}Aligned.out.sam > ${newName}Aligned.out.bam
fi
rm -rf *.sam
if [ -e "${newName}Aligned.sortedByCoord.out.bam" ] ; then
    mv ${newName}Aligned.sortedByCoord.out.bam ${newName}.bam
elif [ -e "${newName}Aligned.out.bam" ] ; then
    mv ${newName}Aligned.out.bam ${newName}.bam
fi

samtools flagstat ${newName}.bam > ${newName}.flagstat.txt
"""


}


process STAR_Module_STAR_Summary {

publishDir params.outdir, overwrite: true, mode: 'copy',
	saveAs: {filename ->
	if (filename =~ /.*.(out|tab)$/) "star/$filename"
}

input:
 set val(name), file(alignSum) from g217_16_outputFileOut_g217_18.groupTuple()
 set val(name), file(LogOut) from g217_16_logOut_g217_18.groupTuple()
 set val(name), file(progressOut) from g217_16_progressOut_g217_18.groupTuple()
 set val(name), file(TabOut) from g217_16_outputFileTab_g217_18.groupTuple()

output:
 file "*.tsv"  into g217_18_outputFile_g217_11
 file "*.{out,tab}"  into g217_18_logOut_g_177
 val "star_alignment_sum"  into g217_18_name_g217_11

shell:
'''
#!/usr/bin/env perl
use List::Util qw[min max];
use strict;
use File::Basename;
use Getopt::Long;
use Pod::Usage; 
use Data::Dumper;

my %tsv;
my @headers = ();
my $name = "!{name}";

# merge output files 
`cat !{alignSum} >${name}_Merged_Log.final.out`;
`cat !{LogOut} >${name}_Merged_Log.out`;
`cat !{progressOut} >${name}_Merged_Log.progress.out`;
`cat !{TabOut} >${name}_Merged_SJ.out.tab`;

alteredAligned();

my @keys = keys %tsv;
my $summary = "$name"."_star_sum.tsv";
my $header_string = join("\\t", @headers);
`echo "$header_string" > $summary`;
foreach my $key (@keys){
	my $values = join("\\t", @{ $tsv{$key} });
	`echo "$values" >> $summary`;
}


sub alteredAligned
{
	my @files = qw(!{alignSum});
	my $multimappedSum;
	my $alignedSum;
	my $inputCountSum;
	push(@headers, "Sample");
    push(@headers, "Total Reads");
	push(@headers, "Multimapped Reads Aligned (STAR)");
	push(@headers, "Unique Reads Aligned (STAR)");
	foreach my $file (@files){
		my $multimapped;
		my $aligned;
		my $inputCount;
		chomp($inputCount = `cat $file | grep 'Number of input reads' | awk '{sum+=\\$6} END {print sum}'`);
		chomp($aligned = `cat $file | grep 'Uniquely mapped reads number' | awk '{sum+=\\$6} END {print sum}'`);
		chomp($multimapped = `cat $file | grep 'Number of reads mapped to multiple loci' | awk '{sum+=\\$9} END {print sum}'`);
		$multimappedSum += int($multimapped);
        $alignedSum += int($aligned);
        $inputCountSum += int($inputCount);
	}
	$tsv{$name} = [$name, $inputCountSum];
	push(@{$tsv{$name}}, $multimappedSum);
	push(@{$tsv{$name}}, $alignedSum);
}
'''

}


process STAR_Module_merge_tsv_files_with_same_header {

input:
 file tsv from g217_18_outputFile_g217_11.collect()
 val outputFileName from g217_18_name_g217_11.collect()

output:
 file "${name}.tsv"  into g217_11_outputFileTSV_g_198

errorStrategy 'retry'
maxRetries 3

script:
name = outputFileName[0]
"""    
awk 'FNR==1 && NR!=1 {  getline; } 1 {print} ' *.tsv > ${name}.tsv
"""
}


//* autofill
//* platform
if ($HOSTNAME == "ghpcc06.umassrc.org"){
    $TIME = 2000
    $CPU  = 1
    $MEMORY = 8
    $QUEUE = "long"
}
//* platform
//* autofill

process STAR_Module_Merge_Bam {

publishDir params.outdir, overwrite: true, mode: 'copy',
	saveAs: {filename ->
	if (filename =~ /.*_sorted.*bai$/) "star/$filename"
	else if (filename =~ /.*_sorted.*bam$/) "star/$filename"
}

input:
 set val(oldname), file(bamfiles) from g217_16_mapped_reads_g217_14.groupTuple()

output:
 set val(oldname), file("${oldname}.bam")  into g217_14_merged_bams
 set val(oldname), file("*_sorted*bai")  into g217_14_bam_index
 set val(oldname), file("*_sorted*bam")  into g217_14_sorted_bam_g221_121, g217_14_sorted_bam_g221_122, g217_14_sorted_bam_g221_123, g217_14_sorted_bam_g221_124, g217_14_sorted_bam_g221_126

errorStrategy 'retry'
maxRetries 2

shell:
'''
num=$(echo "!{bamfiles.join(" ")}" | awk -F" " '{print NF-1}')
if [ "${num}" -gt 0 ]; then
    samtools merge !{oldname}.bam !{bamfiles.join(" ")} && samtools sort -o !{oldname}_sorted.bam !{oldname}.bam && samtools index !{oldname}_sorted.bam
else
    mv !{bamfiles.join(" ")} !{oldname}.bam 2>/dev/null || true
    samtools sort  -o !{oldname}_sorted.bam !{oldname}.bam && samtools index !{oldname}_sorted.bam
fi
'''
}

//* params.gtf =  ""  //* @input


process BAM_Analysis_STAR_featureCounts {

publishDir params.outdir, overwrite: true, mode: 'copy',
	saveAs: {filename ->
	if (filename =~ /.*$/) "featureCounts_after_STAR/$filename"
}

input:
 set val(name), file(bam) from g217_14_sorted_bam_g221_126
 val paired from g_229_mate_g221_126
 each run_params from g221_125_run_parameters_g221_126

output:
 file "*"  into g221_126_outputFileTSV_g221_117

script:
pairText = ""
if (paired == "pair"){
    pairText = "-p"
}

run_name = run_params["run_name"] 
run_parameters = run_params["run_parameters"] 

"""
featureCounts ${pairText} ${run_parameters} -a ${params.gtf} -o ${name}@${run_name}@fCounts.txt ${bam}
## remove first line
sed -i '1d' ${name}@${run_name}@fCounts.txt

"""
}

//* autofill
//* platform
if ($HOSTNAME == "ghpcc06.umassrc.org"){
    $TIME = 30
    $CPU  = 1
    $MEMORY = 10
    $QUEUE = "short"
}
//* platform
//* autofill

process BAM_Analysis_STAR_summary_featureCounts {

publishDir params.outdir, overwrite: true, mode: 'copy',
	saveAs: {filename ->
	if (filename =~ /.*_featureCounts.tsv$/) "featureCounts_after_STAR_summary/$filename"
	else if (filename =~ /.*_featureCounts.sum.tsv$/) "featureCounts_after_STAR_details/$filename"
}

input:
 file featureCountsOut from g221_126_outputFileTSV_g221_117.collect()

output:
 file "*_featureCounts.tsv"  into g221_117_outputFile
 file "*_featureCounts.sum.tsv"  into g221_117_outFileTSV

shell:
'''
#!/usr/bin/env perl

# Step 1: Merge count files
my %tf = ( expected_count => 6 );
my @run_name=();
chomp(my $contents = `ls *@fCounts.txt`);
my @files = split(/[\\n]+/, $contents);
foreach my $file (@files){
    $file=~/(.*)\\@(.*)\\@fCounts\\.txt/;
    my $runname = $2;
    push(@run_name, $runname) unless grep{$_ eq $runname} @run_name;
}


my @expectedCount_ar = ("expected_count");
for($l = 0; $l <= $#run_name; $l++) {
    my $runName = $run_name[$l];
    for($ll = 0; $ll <= $#expectedCount_ar; $ll++) {
        my $expectedCount = $expectedCount_ar[$ll];
    
        my @a=();
        my %b=();
        my %c=();
        my $i=0;
        chomp(my $contents = `ls *\\@${runName}\\@fCounts.txt`);
        my @files = split(/[\\n]+/, $contents);
        foreach my $file (@files){
        $i++;
        $file=~/(.*)\\@${runName}\\@fCounts\\.txt/;
        my $libname = $1; 
        $a[$i]=$libname;
        open IN, $file;
            $_=<IN>;
            while(<IN>){
                my @v=split; 
                $b{$v[0]}{$i}=$v[$tf{$expectedCount}];
                $c{$v[0]}=$v[5]; #length column
            }
            close IN;
        }
        my $outfile="$runName"."_featureCounts.tsv";
        open OUT, ">$outfile";
        if ($runName eq "transcript_id") {
            print OUT "transcript\tlength";
        } else {
            print OUT "gene\tlength";
        }
    
        for(my $j=1;$j<=$i;$j++) {
            print OUT "\t$a[$j]";
        }
        print OUT "\n";
    
        foreach my $key (keys %b) {
            print OUT "$key\t$c{$key}";
            for(my $j=1;$j<=$i;$j++){
                print OUT "\t$b{$key}{$j}";
            }
            print OUT "\n";
        }
        close OUT;
         
    }
}


	

# Step 2: Merge summary files
for($l = 0; $l <= $#run_name; $l++) {
    my $runName = $run_name[$l];
    my @a=();
    my %b=();
    my $i=0;
    chomp(my $contents = `ls *\\@${runName}\\@fCounts.txt.summary`);
    my @files = split(/[\\n]+/, $contents);
    foreach my $file (@files){
        $i++;
        $file=~/(.*)\\@${runName}\\@fCounts\\.txt\\.summary/;
        my $libname = $1; 
        $a[$i]=$libname;
        open IN, $file;
        $_=<IN>;
        while(<IN>){
            my @v=split; 
            $b{$v[0]}{$i}=$v[1];
        }
        close IN;
    }
    my $outfile="$runName"."_featureCounts.sum.tsv";
    open OUT, ">$outfile";
    print OUT "criteria";
    for(my $j=1;$j<=$i;$j++) {
        print OUT "\t$a[$j]";
    }
    print OUT "\n";
    
    foreach my $key (keys %b) {
        print OUT "$key";
        for(my $j=1;$j<=$i;$j++){
            print OUT "\t$b{$key}{$j}";
        }
        print OUT "\n";
    }
    close OUT;
}

'''
}

//* params.genome_sizes =  ""  //* @input

//* autofill
if ($HOSTNAME == "default"){
    $CPU  = 1
    $MEMORY = 48
}
//* platform
if ($HOSTNAME == "ghpcc06.umassrc.org"){
    $TIME = 240
    $CPU  = 1
    $MEMORY = 48
    $QUEUE = "short"
} 
//* platform
//* autofill

process BAM_Analysis_STAR_UCSC_BAM2BigWig_converter {

publishDir params.outdir, overwrite: true, mode: 'copy',
	saveAs: {filename ->
	if (filename =~ /.*.bw$/) "bigwig_star/$filename"
}

input:
 set val(name), file(bam) from g217_14_sorted_bam_g221_124

output:
 file "*.bw"  into g221_124_outputFileBw

when:
(params.run_BigWig_Conversion && (params.run_BigWig_Conversion == "yes")) || !params.run_BigWig_Conversion

script:
nameAll = bam.toString()
if (nameAll.contains('_sorted.bam')) {
    runSamtools = "samtools index ${nameAll}"
    nameFinal = nameAll
} else {
    runSamtools = "samtools sort -o ${name}_sorted.bam $bam && samtools index ${name}_sorted.bam "
    nameFinal = "${name}_sorted.bam"
}

"""
$runSamtools
bedtools genomecov -split -bg -ibam ${nameFinal} -g ${params.genome_sizes} > ${name}.bg 
wigToBigWig -clip -itemsPerSlot=1 ${name}.bg ${params.genome_sizes} ${name}.bw 
"""
}

igv_extention_factor = params.BAM_Analysis_STAR_IGV_BAM2TDF_converter.igv_extention_factor
igv_window_size = params.BAM_Analysis_STAR_IGV_BAM2TDF_converter.igv_window_size

//* params.genome =  ""  //* @input

//* autofill
if ($HOSTNAME == "default"){
    $CPU  = 1
    $MEMORY = 32
}
//* platform
if ($HOSTNAME == "ghpcc06.umassrc.org"){
    $TIME = 400
    $CPU  = 1
    $MEMORY = 32
    $QUEUE = "long"
} 
//* platform
//* autofill

process BAM_Analysis_STAR_IGV_BAM2TDF_converter {

publishDir params.outdir, overwrite: true, mode: 'copy',
	saveAs: {filename ->
	if (filename =~ /.*.tdf$/) "igvtools_star/$filename"
}

input:
 val mate from g_229_mate_g221_123
 set val(name), file(bam) from g217_14_sorted_bam_g221_123

output:
 file "*.tdf"  into g221_123_outputFileOut

when:
(params.run_IGV_TDF_Conversion && (params.run_IGV_TDF_Conversion == "yes")) || !params.run_IGV_TDF_Conversion

script:
pairedText = (params.nucleicAcidType == "dna" && mate == "pair") ? " --pairs " : ""
nameAll = bam.toString()
if (nameAll.contains('_sorted.bam')) {
    runSamtools = "samtools index ${nameAll}"
    nameFinal = nameAll
} else {
    runSamtools = "samtools sort -o ${name}_sorted.bam $bam && samtools index ${name}_sorted.bam "
    nameFinal = "${name}_sorted.bam"
}
"""
$runSamtools
igvtools count -w ${igv_window_size} -e ${igv_extention_factor} ${pairedText} ${nameFinal} ${name}.tdf ${params.genome_sizes}
"""
}

//* params.bed =  ""  //* @input

process BAM_Analysis_STAR_RSeQC {

publishDir params.outdir, overwrite: true, mode: 'copy',
	saveAs: {filename ->
	if (filename =~ /result\/.*.out$/) "rseqc_star/$filename"
}

input:
 set val(name), file(bam) from g217_14_sorted_bam_g221_122

output:
 file "result/*.out"  into g221_122_outputFileOut_g221_95, g221_122_outputFileOut_g_177

when:
(params.run_RSeQC && (params.run_RSeQC == "yes")) || !params.run_RSeQC

script:
"""
mkdir result
read_distribution.py  -i ${bam} -r ${params.bed}> result/RSeQC.${name}.out
"""
}


process BAM_Analysis_STAR_RSeQC_Summary {

publishDir params.outdir, overwrite: true, mode: 'copy',
	saveAs: {filename ->
	if (filename =~ /.*.tsv$/) "rseqc_summary_star/$filename"
}

input:
 file rseqcOut from g221_122_outputFileOut_g221_95.collect()
 val mate from g_229_mate_g221_95

output:
 file "*.tsv"  into g221_95_outputFileTSV

shell:
'''
#!/usr/bin/env perl
use List::Util qw[min max];
use strict;
use File::Basename;
use Getopt::Long;
use Pod::Usage; 
use Data::Dumper;

my $indir = $ENV{'PWD'};
my $outd = $ENV{'PWD'};
my @files = ();
my @outtypes = ("RSeQC");
my @order=( "Total Reads", "Total Tags" , "Total Assigned Tags", "CDS_Exons", "5'UTR_Exons", "3'UTR_Exons", "Introns", "TSS_up_1kb", "TSS_up_5kb", "TSS_up_10kb", "TES_down_1kb", "TES_down_5kb", "TES_down_10kb");
my %lines=(
  "Total Reads" => 1,
  "Total Tags" => 1,
  "Total Assigned Tags" => 1,
  "CDS_Exons" => 2,
  "5'UTR_Exons" => 2,
  "3'UTR_Exons" => 2,
  "Introns" => 2,
  "TSS_up_1kb" => 2,
  "TSS_up_5kb" => 2,
  "TSS_up_10kb" => 2,
  "TES_down_1kb" => 2,
  "TES_down_5kb" => 2,
  "TES_down_10kb" => 2
);


foreach my $outtype (@outtypes)
{

my $ext=".out";
@files = <$indir/$outtype*$ext>;

my @rowheaders=();
my @libs=();
my %vals=();
my %normvals=();
my $type = "rsem";

foreach my $d (@files){
  my $libname=basename($d, $ext);
  $libname=~s/RSeQC.//g;
  $libname=~s/rsem.out.//g;
  $libname=~s/.genome//g;
  print $libname."\\n";
  push(@libs, $libname); 
  getVals($d, $libname, \\%vals, \\%normvals, \\%lines);
}
#print Dumper(%vals);
#print Dumper(%normvals);

my $sizemetrics = keys %vals;
write_results("$outd/$outtype.$type.counts.tsv", \\@libs,\\%vals, \\@order, "region") if ($sizemetrics>0);
write_results("$outd/$outtype.$type.tagskb.tsv", \\@libs,\\%normvals, \\@order, "region") if ($sizemetrics>0);

}

sub write_results
{
  my ($outfile, $libs, $vals, $order, $name )=@_;
  open(OUT, ">$outfile");
  print OUT "$name\\t".join("\\t", @{$libs})."\\n";

  my $lib=${$libs}[0];
  foreach my $key ( @order )
  {
    if (exists ${$vals}{$lib}{$key}) {
    print OUT $key;
    foreach my $lib (@{$libs})
    {
      print OUT "\\t".${$vals}{$lib}{$key};
    } 
    print OUT "\\n";
    }
  }
  close(OUT);
}

sub getVals{
  my ($filename, $libname, $vals, $normvals, $lines)=@_;
  if (-e $filename){
     open(IN, $filename);
     while(my $line=<IN>)
     {
       chomp($line);
       my @vals_arr=split(/\\s{2,}/,$line);
       if (exists ${$lines}{$vals_arr[0]}) {
         my $idx=${$lines}{$vals_arr[0]};
         ${$vals}{$libname}{$vals_arr[0]}=$vals_arr[$idx] if (exists $vals_arr[$idx]);
         if ($idx==2) {
             ${$normvals}{$libname}{$vals_arr[0]}=$vals_arr[3] if (exists $vals_arr[3]);
         }
       }
     } 
  }
  
}
'''

}

//* params.pdfbox_path =  ""  //* @input
//* autofill
if ($HOSTNAME == "default"){
    $CPU  = 1
    $MEMORY = 32
}
//* platform
if ($HOSTNAME == "ghpcc06.umassrc.org"){
    $TIME = 240
    $CPU  = 1
    $MEMORY = 32
    $QUEUE = "short"
}
//* platform
//* autofill

process BAM_Analysis_STAR_Picard {

input:
 set val(name), file(bam) from g217_14_sorted_bam_g221_121

output:
 file "*_metrics"  into g221_121_outputFileOut_g221_82
 file "results/*.pdf"  into g221_121_outputFilePdf_g221_82

when:
(params.run_Picard_CollectMultipleMetrics && (params.run_Picard_CollectMultipleMetrics == "yes")) || !params.run_Picard_CollectMultipleMetrics

script:
"""
picard CollectMultipleMetrics OUTPUT=${name}_multiple.out VALIDATION_STRINGENCY=LENIENT INPUT=${bam}
mkdir results && java -jar ${params.pdfbox_path} PDFMerger *.pdf results/${name}_multi_metrics.pdf
"""
}


process BAM_Analysis_STAR_Picard_Summary {

publishDir params.outdir, overwrite: true, mode: 'copy',
	saveAs: {filename ->
	if (filename =~ /.*.tsv$/) "picard_summary_star/$filename"
	else if (filename =~ /results\/.*.pdf$/) "picard_summary_pdf_star/$filename"
}

input:
 file picardOut from g221_121_outputFileOut_g221_82.collect()
 val mate from g_229_mate_g221_82
 file picardPdf from g221_121_outputFilePdf_g221_82.collect()

output:
 file "*.tsv"  into g221_82_outputFileTSV
 file "results/*.pdf"  into g221_82_outputFilePdf

shell:
'''
#!/usr/bin/env perl
use List::Util qw[min max];
use strict;
use File::Basename;
use Getopt::Long;
use Pod::Usage; 
use Data::Dumper;

system("mkdir results && mv *.pdf results/. ");

my $indir = $ENV{'PWD'};
my $outd = $ENV{'PWD'};
my @files = ();
my @outtypes = ("CollectRnaSeqMetrics", "alignment_summary_metrics", "base_distribution_by_cycle_metrics", "insert_size_metrics", "quality_by_cycle_metrics", "quality_distribution_metrics" );

foreach my $outtype (@outtypes)
{
my $ext="_multiple.out";
$ext.=".$outtype" if ($outtype ne "CollectRnaSeqMetrics");
@files = <$indir/*$ext>;

my @rowheaders=();
my @libs=();
my %metricvals=();
my %histvals=();

my $pdffile="";
my $libname="";
foreach my $d (@files){
  my $libname=basename($d, $ext);
  print $libname."\\n";
  push(@libs, $libname); 
  getMetricVals($d, $libname, \\%metricvals, \\%histvals, \\@rowheaders);
}

my $sizemetrics = keys %metricvals;
write_results("$outd/$outtype.stats.tsv", \\@libs,\\%metricvals, \\@rowheaders, "metric") if ($sizemetrics>0);
my $sizehist = keys %histvals;
write_results("$outd/$outtype.hist.tsv", \\@libs,\\%histvals, "none", "nt") if ($sizehist>0);

}

sub write_results
{
  my ($outfile, $libs, $vals, $rowheaders, $name )=@_;
  open(OUT, ">$outfile");
  print OUT "$name\\t".join("\\t", @{$libs})."\\n";
  my $size=0;
  $size=scalar(@{${$vals}{${$libs}[0]}}) if(exists ${$libs}[0] and exists ${$vals}{${$libs}[0]} );
  
  for (my $i=0; $i<$size;$i++)
  { 
    my $rowname=$i;
    $rowname = ${$rowheaders}[$i] if ($name=~/metric/);
    print OUT $rowname;
    foreach my $lib (@{$libs})
    {
      print OUT "\\t".${${$vals}{$lib}}[$i];
    } 
    print OUT "\\n";
  }
  close(OUT);
}

sub getMetricVals{
  my ($filename, $libname, $metricvals, $histvals,$rowheaders)=@_;
  if (-e $filename){
     my $nextisheader=0;
     my $nextisvals=0;
     my $nexthist=0;
     open(IN, $filename);
     while(my $line=<IN>)
     {
       chomp($line);
       @{$rowheaders}=split(/\\t/, $line) if ($nextisheader && !scalar(@{$rowheaders})); 
       if ($nextisvals) {
         @{${$metricvals}{$libname}}=split(/\\t/, $line);
         $nextisvals=0;
       }
       if($nexthist){
          my @vals=split(/[\\s\\t]+/,$line); 
          push(@{${$histvals}{$libname}}, $vals[1]) if (exists $vals[1]);
       }
       $nextisvals=1 if ($nextisheader); $nextisheader=0;
       $nextisheader=1 if ($line=~/METRICS CLASS/);
       $nexthist=1 if ($line=~/normalized_position/);
     } 
  }
  
}
'''

}


//* autofill
//* platform
if ($HOSTNAME == "ghpcc06.umassrc.org"){
    $TIME = 2000
    $CPU  = 1
    $MEMORY = 8
    $QUEUE = "long"
}
//* platform
//* autofill

process STAR_Module_merge_transcriptome_bam {

input:
 set val(oldname), file(bamfiles) from g217_16_transcriptome_bam_g217_15.groupTuple()

output:
 set val(oldname), file("${oldname}.bam")  into g217_15_merged_bams
 set val(oldname), file("*_sorted*bai")  into g217_15_bam_index
 set val(oldname), file("*_sorted*bam")  into g217_15_sorted_bam

errorStrategy 'retry'
maxRetries 2

shell:
'''
num=$(echo "!{bamfiles.join(" ")}" | awk -F" " '{print NF-1}')
if [ "${num}" -gt 0 ]; then
    samtools merge !{oldname}.bam !{bamfiles.join(" ")} && samtools sort -o !{oldname}_sorted.bam !{oldname}.bam && samtools index !{oldname}_sorted.bam
else
    mv !{bamfiles.join(" ")} !{oldname}.bam 2>/dev/null || true
    samtools sort  -o !{oldname}_sorted.bam !{oldname}.bam && samtools index !{oldname}_sorted.bam
fi
'''
}

g227_8_genomeIndexPath_g216_11= g227_8_genomeIndexPath_g216_11.ifEmpty([""]) 

//* params.hisat2_index =  ""  //* @input


//* autofill
if ($HOSTNAME == "default"){
    $CPU  = 3
    $MEMORY = 32
}
//* platform
if ($HOSTNAME == "ghpcc06.umassrc.org"){
    $TIME = 240
    $CPU  = 4
    $MEMORY = 32
    $QUEUE = "short"
} 
//* platform
//* autofill

process HISAT2_Module_Map_HISAT2 {

input:
 val mate from g_229_mate_g216_11
 set val(name), file(reads) from g_127_reads_g216_11
 val hisat2index from g227_8_genomeIndexPath_g216_11

output:
 set val(name), file("${newName}.bam")  into g216_11_mapped_reads_g216_13
 set val(name), file("${newName}.align_summary.txt")  into g216_11_outputFileTxt_g216_2
 set val(name), file("${newName}.flagstat.txt")  into g216_11_outputFileOut

when:
(params.run_HISAT2 && (params.run_HISAT2 == "yes")) || !params.run_HISAT2

script:
HISAT2_parameters = params.HISAT2_Module_Map_HISAT2.HISAT2_parameters
nameAll = reads.toString()
nameArray = nameAll.split(' ')
file2 = ""

if (nameAll.contains('.gz')) {
    newName =  nameArray[0] - ~/(\.fastq.gz)?(\.fq.gz)?$/
    file1 =  nameArray[0] - '.gz' 
    if (mate == "pair") {file2 =  nameArray[1] - '.gz'}
    runGzip = "ls *.gz | xargs -i echo gzip -df {} | sh"
} else {
    newName =  nameArray[0] - ~/(\.fastq)?(\.fq)?$/
    file1 =  nameArray[0]
    if (mate == "pair") {file2 =  nameArray[1]}
    runGzip = ''
}

"""
$runGzip
if [ "${mate}" == "pair" ]; then
    hisat2 ${HISAT2_parameters} -x ${params.hisat2_index} -1 ${file1} -2 ${file2} -S ${newName}.sam &> ${newName}.align_summary.txt
else
    hisat2 ${HISAT2_parameters} -x ${params.hisat2_index} -U ${file1} -S ${newName}.sam &> ${newName}.align_summary.txt
fi
samtools view -bS ${newName}.sam > ${newName}.bam
samtools flagstat ${newName}.bam > ${newName}.flagstat.txt
"""

}


//* autofill
//* platform
if ($HOSTNAME == "ghpcc06.umassrc.org"){
    $TIME = 2000
    $CPU  = 1
    $MEMORY = 8
    $QUEUE = "long"
}
//* platform
//* autofill

process HISAT2_Module_Merge_Bam {

publishDir params.outdir, overwrite: true, mode: 'copy',
	saveAs: {filename ->
	if (filename =~ /.*_sorted.*bam$/) "hisat2/$filename"
}

input:
 set val(oldname), file(bamfiles) from g216_11_mapped_reads_g216_13.groupTuple()

output:
 set val(oldname), file("${oldname}.bam")  into g216_13_merged_bams
 set val(oldname), file("*_sorted*bai")  into g216_13_bam_index
 set val(oldname), file("*_sorted*bam")  into g216_13_sorted_bam_g220_121, g216_13_sorted_bam_g220_122, g216_13_sorted_bam_g220_123, g216_13_sorted_bam_g220_124, g216_13_sorted_bam_g220_126

errorStrategy 'retry'
maxRetries 2

shell:
'''
num=$(echo "!{bamfiles.join(" ")}" | awk -F" " '{print NF-1}')
if [ "${num}" -gt 0 ]; then
    samtools merge !{oldname}.bam !{bamfiles.join(" ")} && samtools sort -o !{oldname}_sorted.bam !{oldname}.bam && samtools index !{oldname}_sorted.bam
else
    mv !{bamfiles.join(" ")} !{oldname}.bam 2>/dev/null || true
    samtools sort  -o !{oldname}_sorted.bam !{oldname}.bam && samtools index !{oldname}_sorted.bam
fi
'''
}

//* params.gtf =  ""  //* @input


process BAM_Analysis_Hisat2_featureCounts {

publishDir params.outdir, overwrite: true, mode: 'copy',
	saveAs: {filename ->
	if (filename =~ /.*$/) "featureCounts_after_hisat2/$filename"
}

input:
 set val(name), file(bam) from g216_13_sorted_bam_g220_126
 val paired from g_229_mate_g220_126
 each run_params from g220_125_run_parameters_g220_126

output:
 file "*"  into g220_126_outputFileTSV_g220_117

script:
pairText = ""
if (paired == "pair"){
    pairText = "-p"
}

run_name = run_params["run_name"] 
run_parameters = run_params["run_parameters"] 

"""
featureCounts ${pairText} ${run_parameters} -a ${params.gtf} -o ${name}@${run_name}@fCounts.txt ${bam}
## remove first line
sed -i '1d' ${name}@${run_name}@fCounts.txt

"""
}

//* autofill
//* platform
if ($HOSTNAME == "ghpcc06.umassrc.org"){
    $TIME = 30
    $CPU  = 1
    $MEMORY = 10
    $QUEUE = "short"
}
//* platform
//* autofill

process BAM_Analysis_Hisat2_summary_featureCounts {

publishDir params.outdir, overwrite: true, mode: 'copy',
	saveAs: {filename ->
	if (filename =~ /.*_featureCounts.tsv$/) "featureCounts_after_hisat2_summary/$filename"
	else if (filename =~ /.*_featureCounts.sum.tsv$/) "featureCounts_after_hisat2_details/$filename"
}

input:
 file featureCountsOut from g220_126_outputFileTSV_g220_117.collect()

output:
 file "*_featureCounts.tsv"  into g220_117_outputFile
 file "*_featureCounts.sum.tsv"  into g220_117_outFileTSV

shell:
'''
#!/usr/bin/env perl

# Step 1: Merge count files
my %tf = ( expected_count => 6 );
my @run_name=();
chomp(my $contents = `ls *@fCounts.txt`);
my @files = split(/[\\n]+/, $contents);
foreach my $file (@files){
    $file=~/(.*)\\@(.*)\\@fCounts\\.txt/;
    my $runname = $2;
    push(@run_name, $runname) unless grep{$_ eq $runname} @run_name;
}


my @expectedCount_ar = ("expected_count");
for($l = 0; $l <= $#run_name; $l++) {
    my $runName = $run_name[$l];
    for($ll = 0; $ll <= $#expectedCount_ar; $ll++) {
        my $expectedCount = $expectedCount_ar[$ll];
    
        my @a=();
        my %b=();
        my %c=();
        my $i=0;
        chomp(my $contents = `ls *\\@${runName}\\@fCounts.txt`);
        my @files = split(/[\\n]+/, $contents);
        foreach my $file (@files){
        $i++;
        $file=~/(.*)\\@${runName}\\@fCounts\\.txt/;
        my $libname = $1; 
        $a[$i]=$libname;
        open IN, $file;
            $_=<IN>;
            while(<IN>){
                my @v=split; 
                $b{$v[0]}{$i}=$v[$tf{$expectedCount}];
                $c{$v[0]}=$v[5]; #length column
            }
            close IN;
        }
        my $outfile="$runName"."_featureCounts.tsv";
        open OUT, ">$outfile";
        if ($runName eq "transcript_id") {
            print OUT "transcript\tlength";
        } else {
            print OUT "gene\tlength";
        }
    
        for(my $j=1;$j<=$i;$j++) {
            print OUT "\t$a[$j]";
        }
        print OUT "\n";
    
        foreach my $key (keys %b) {
            print OUT "$key\t$c{$key}";
            for(my $j=1;$j<=$i;$j++){
                print OUT "\t$b{$key}{$j}";
            }
            print OUT "\n";
        }
        close OUT;
         
    }
}


	

# Step 2: Merge summary files
for($l = 0; $l <= $#run_name; $l++) {
    my $runName = $run_name[$l];
    my @a=();
    my %b=();
    my $i=0;
    chomp(my $contents = `ls *\\@${runName}\\@fCounts.txt.summary`);
    my @files = split(/[\\n]+/, $contents);
    foreach my $file (@files){
        $i++;
        $file=~/(.*)\\@${runName}\\@fCounts\\.txt\\.summary/;
        my $libname = $1; 
        $a[$i]=$libname;
        open IN, $file;
        $_=<IN>;
        while(<IN>){
            my @v=split; 
            $b{$v[0]}{$i}=$v[1];
        }
        close IN;
    }
    my $outfile="$runName"."_featureCounts.sum.tsv";
    open OUT, ">$outfile";
    print OUT "criteria";
    for(my $j=1;$j<=$i;$j++) {
        print OUT "\t$a[$j]";
    }
    print OUT "\n";
    
    foreach my $key (keys %b) {
        print OUT "$key";
        for(my $j=1;$j<=$i;$j++){
            print OUT "\t$b{$key}{$j}";
        }
        print OUT "\n";
    }
    close OUT;
}

'''
}

//* params.genome_sizes =  ""  //* @input

//* autofill
if ($HOSTNAME == "default"){
    $CPU  = 1
    $MEMORY = 48
}
//* platform
if ($HOSTNAME == "ghpcc06.umassrc.org"){
    $TIME = 240
    $CPU  = 1
    $MEMORY = 48
    $QUEUE = "short"
} 
//* platform
//* autofill

process BAM_Analysis_Hisat2_UCSC_BAM2BigWig_converter {

publishDir params.outdir, overwrite: true, mode: 'copy',
	saveAs: {filename ->
	if (filename =~ /.*.bw$/) "bigwig_hisat2/$filename"
}

input:
 set val(name), file(bam) from g216_13_sorted_bam_g220_124

output:
 file "*.bw"  into g220_124_outputFileBw

when:
(params.run_BigWig_Conversion && (params.run_BigWig_Conversion == "yes")) || !params.run_BigWig_Conversion

script:
nameAll = bam.toString()
if (nameAll.contains('_sorted.bam')) {
    runSamtools = "samtools index ${nameAll}"
    nameFinal = nameAll
} else {
    runSamtools = "samtools sort -o ${name}_sorted.bam $bam && samtools index ${name}_sorted.bam "
    nameFinal = "${name}_sorted.bam"
}

"""
$runSamtools
bedtools genomecov -split -bg -ibam ${nameFinal} -g ${params.genome_sizes} > ${name}.bg 
wigToBigWig -clip -itemsPerSlot=1 ${name}.bg ${params.genome_sizes} ${name}.bw 
"""
}

igv_extention_factor = params.BAM_Analysis_Hisat2_IGV_BAM2TDF_converter.igv_extention_factor
igv_window_size = params.BAM_Analysis_Hisat2_IGV_BAM2TDF_converter.igv_window_size

//* params.genome =  ""  //* @input

//* autofill
if ($HOSTNAME == "default"){
    $CPU  = 1
    $MEMORY = 32
}
//* platform
if ($HOSTNAME == "ghpcc06.umassrc.org"){
    $TIME = 400
    $CPU  = 1
    $MEMORY = 32
    $QUEUE = "long"
} 
//* platform
//* autofill

process BAM_Analysis_Hisat2_IGV_BAM2TDF_converter {

publishDir params.outdir, overwrite: true, mode: 'copy',
	saveAs: {filename ->
	if (filename =~ /.*.tdf$/) "igvtools_hisat2/$filename"
}

input:
 val mate from g_229_mate_g220_123
 set val(name), file(bam) from g216_13_sorted_bam_g220_123

output:
 file "*.tdf"  into g220_123_outputFileOut

when:
(params.run_IGV_TDF_Conversion && (params.run_IGV_TDF_Conversion == "yes")) || !params.run_IGV_TDF_Conversion

script:
pairedText = (params.nucleicAcidType == "dna" && mate == "pair") ? " --pairs " : ""
nameAll = bam.toString()
if (nameAll.contains('_sorted.bam')) {
    runSamtools = "samtools index ${nameAll}"
    nameFinal = nameAll
} else {
    runSamtools = "samtools sort -o ${name}_sorted.bam $bam && samtools index ${name}_sorted.bam "
    nameFinal = "${name}_sorted.bam"
}
"""
$runSamtools
igvtools count -w ${igv_window_size} -e ${igv_extention_factor} ${pairedText} ${nameFinal} ${name}.tdf ${params.genome_sizes}
"""
}

//* params.bed =  ""  //* @input

process BAM_Analysis_Hisat2_RSeQC {

publishDir params.outdir, overwrite: true, mode: 'copy',
	saveAs: {filename ->
	if (filename =~ /result\/.*.out$/) "rseqc_hisat2/$filename"
}

input:
 set val(name), file(bam) from g216_13_sorted_bam_g220_122

output:
 file "result/*.out"  into g220_122_outputFileOut_g220_95, g220_122_outputFileOut_g_177

when:
(params.run_RSeQC && (params.run_RSeQC == "yes")) || !params.run_RSeQC

script:
"""
mkdir result
read_distribution.py  -i ${bam} -r ${params.bed}> result/RSeQC.${name}.out
"""
}


process BAM_Analysis_Hisat2_RSeQC_Summary {

publishDir params.outdir, overwrite: true, mode: 'copy',
	saveAs: {filename ->
	if (filename =~ /.*.tsv$/) "rseqc_summary_hisat2/$filename"
}

input:
 file rseqcOut from g220_122_outputFileOut_g220_95.collect()
 val mate from g_229_mate_g220_95

output:
 file "*.tsv"  into g220_95_outputFileTSV

shell:
'''
#!/usr/bin/env perl
use List::Util qw[min max];
use strict;
use File::Basename;
use Getopt::Long;
use Pod::Usage; 
use Data::Dumper;

my $indir = $ENV{'PWD'};
my $outd = $ENV{'PWD'};
my @files = ();
my @outtypes = ("RSeQC");
my @order=( "Total Reads", "Total Tags" , "Total Assigned Tags", "CDS_Exons", "5'UTR_Exons", "3'UTR_Exons", "Introns", "TSS_up_1kb", "TSS_up_5kb", "TSS_up_10kb", "TES_down_1kb", "TES_down_5kb", "TES_down_10kb");
my %lines=(
  "Total Reads" => 1,
  "Total Tags" => 1,
  "Total Assigned Tags" => 1,
  "CDS_Exons" => 2,
  "5'UTR_Exons" => 2,
  "3'UTR_Exons" => 2,
  "Introns" => 2,
  "TSS_up_1kb" => 2,
  "TSS_up_5kb" => 2,
  "TSS_up_10kb" => 2,
  "TES_down_1kb" => 2,
  "TES_down_5kb" => 2,
  "TES_down_10kb" => 2
);


foreach my $outtype (@outtypes)
{

my $ext=".out";
@files = <$indir/$outtype*$ext>;

my @rowheaders=();
my @libs=();
my %vals=();
my %normvals=();
my $type = "rsem";

foreach my $d (@files){
  my $libname=basename($d, $ext);
  $libname=~s/RSeQC.//g;
  $libname=~s/rsem.out.//g;
  $libname=~s/.genome//g;
  print $libname."\\n";
  push(@libs, $libname); 
  getVals($d, $libname, \\%vals, \\%normvals, \\%lines);
}
#print Dumper(%vals);
#print Dumper(%normvals);

my $sizemetrics = keys %vals;
write_results("$outd/$outtype.$type.counts.tsv", \\@libs,\\%vals, \\@order, "region") if ($sizemetrics>0);
write_results("$outd/$outtype.$type.tagskb.tsv", \\@libs,\\%normvals, \\@order, "region") if ($sizemetrics>0);

}

sub write_results
{
  my ($outfile, $libs, $vals, $order, $name )=@_;
  open(OUT, ">$outfile");
  print OUT "$name\\t".join("\\t", @{$libs})."\\n";

  my $lib=${$libs}[0];
  foreach my $key ( @order )
  {
    if (exists ${$vals}{$lib}{$key}) {
    print OUT $key;
    foreach my $lib (@{$libs})
    {
      print OUT "\\t".${$vals}{$lib}{$key};
    } 
    print OUT "\\n";
    }
  }
  close(OUT);
}

sub getVals{
  my ($filename, $libname, $vals, $normvals, $lines)=@_;
  if (-e $filename){
     open(IN, $filename);
     while(my $line=<IN>)
     {
       chomp($line);
       my @vals_arr=split(/\\s{2,}/,$line);
       if (exists ${$lines}{$vals_arr[0]}) {
         my $idx=${$lines}{$vals_arr[0]};
         ${$vals}{$libname}{$vals_arr[0]}=$vals_arr[$idx] if (exists $vals_arr[$idx]);
         if ($idx==2) {
             ${$normvals}{$libname}{$vals_arr[0]}=$vals_arr[3] if (exists $vals_arr[3]);
         }
       }
     } 
  }
  
}
'''

}

//* params.pdfbox_path =  ""  //* @input
//* autofill
if ($HOSTNAME == "default"){
    $CPU  = 1
    $MEMORY = 32
}
//* platform
if ($HOSTNAME == "ghpcc06.umassrc.org"){
    $TIME = 240
    $CPU  = 1
    $MEMORY = 32
    $QUEUE = "short"
}
//* platform
//* autofill

process BAM_Analysis_Hisat2_Picard {

input:
 set val(name), file(bam) from g216_13_sorted_bam_g220_121

output:
 file "*_metrics"  into g220_121_outputFileOut_g220_82
 file "results/*.pdf"  into g220_121_outputFilePdf_g220_82

when:
(params.run_Picard_CollectMultipleMetrics && (params.run_Picard_CollectMultipleMetrics == "yes")) || !params.run_Picard_CollectMultipleMetrics

script:
"""
picard CollectMultipleMetrics OUTPUT=${name}_multiple.out VALIDATION_STRINGENCY=LENIENT INPUT=${bam}
mkdir results && java -jar ${params.pdfbox_path} PDFMerger *.pdf results/${name}_multi_metrics.pdf
"""
}


process BAM_Analysis_Hisat2_Picard_Summary {

publishDir params.outdir, overwrite: true, mode: 'copy',
	saveAs: {filename ->
	if (filename =~ /.*.tsv$/) "picard_summary_hisat2/$filename"
	else if (filename =~ /results\/.*.pdf$/) "picard_summary_pdf_hisat2/$filename"
}

input:
 file picardOut from g220_121_outputFileOut_g220_82.collect()
 val mate from g_229_mate_g220_82
 file picardPdf from g220_121_outputFilePdf_g220_82.collect()

output:
 file "*.tsv"  into g220_82_outputFileTSV
 file "results/*.pdf"  into g220_82_outputFilePdf

shell:
'''
#!/usr/bin/env perl
use List::Util qw[min max];
use strict;
use File::Basename;
use Getopt::Long;
use Pod::Usage; 
use Data::Dumper;

system("mkdir results && mv *.pdf results/. ");

my $indir = $ENV{'PWD'};
my $outd = $ENV{'PWD'};
my @files = ();
my @outtypes = ("CollectRnaSeqMetrics", "alignment_summary_metrics", "base_distribution_by_cycle_metrics", "insert_size_metrics", "quality_by_cycle_metrics", "quality_distribution_metrics" );

foreach my $outtype (@outtypes)
{
my $ext="_multiple.out";
$ext.=".$outtype" if ($outtype ne "CollectRnaSeqMetrics");
@files = <$indir/*$ext>;

my @rowheaders=();
my @libs=();
my %metricvals=();
my %histvals=();

my $pdffile="";
my $libname="";
foreach my $d (@files){
  my $libname=basename($d, $ext);
  print $libname."\\n";
  push(@libs, $libname); 
  getMetricVals($d, $libname, \\%metricvals, \\%histvals, \\@rowheaders);
}

my $sizemetrics = keys %metricvals;
write_results("$outd/$outtype.stats.tsv", \\@libs,\\%metricvals, \\@rowheaders, "metric") if ($sizemetrics>0);
my $sizehist = keys %histvals;
write_results("$outd/$outtype.hist.tsv", \\@libs,\\%histvals, "none", "nt") if ($sizehist>0);

}

sub write_results
{
  my ($outfile, $libs, $vals, $rowheaders, $name )=@_;
  open(OUT, ">$outfile");
  print OUT "$name\\t".join("\\t", @{$libs})."\\n";
  my $size=0;
  $size=scalar(@{${$vals}{${$libs}[0]}}) if(exists ${$libs}[0] and exists ${$vals}{${$libs}[0]} );
  
  for (my $i=0; $i<$size;$i++)
  { 
    my $rowname=$i;
    $rowname = ${$rowheaders}[$i] if ($name=~/metric/);
    print OUT $rowname;
    foreach my $lib (@{$libs})
    {
      print OUT "\\t".${${$vals}{$lib}}[$i];
    } 
    print OUT "\\n";
  }
  close(OUT);
}

sub getMetricVals{
  my ($filename, $libname, $metricvals, $histvals,$rowheaders)=@_;
  if (-e $filename){
     my $nextisheader=0;
     my $nextisvals=0;
     my $nexthist=0;
     open(IN, $filename);
     while(my $line=<IN>)
     {
       chomp($line);
       @{$rowheaders}=split(/\\t/, $line) if ($nextisheader && !scalar(@{$rowheaders})); 
       if ($nextisvals) {
         @{${$metricvals}{$libname}}=split(/\\t/, $line);
         $nextisvals=0;
       }
       if($nexthist){
          my @vals=split(/[\\s\\t]+/,$line); 
          push(@{${$histvals}{$libname}}, $vals[1]) if (exists $vals[1]);
       }
       $nextisvals=1 if ($nextisheader); $nextisheader=0;
       $nextisheader=1 if ($line=~/METRICS CLASS/);
       $nexthist=1 if ($line=~/normalized_position/);
     } 
  }
  
}
'''

}


process HISAT2_Module_HISAT2_Summary {

input:
 set val(name), file(alignSum) from g216_11_outputFileTxt_g216_2.groupTuple()

output:
 file "*.tsv"  into g216_2_outputFile_g216_10
 val "hisat2_alignment_sum"  into g216_2_name_g216_10

shell:
'''
#!/usr/bin/env perl
use List::Util qw[min max];
use strict;
use File::Basename;
use Getopt::Long;
use Pod::Usage; 
use Data::Dumper;

my %tsv;
my @headers = ();
my $name = "!{name}";


alteredAligned();

my @keys = keys %tsv;
my $summary = "$name"."_hisat_sum.tsv";
my $header_string = join("\\t", @headers);
`echo "$header_string" > $summary`;
foreach my $key (@keys){
	my $values = join("\\t", @{ $tsv{$key} });
	`echo "$values" >> $summary`;
}


sub alteredAligned
{
	my @files = qw(!{alignSum});
	my $multimappedSum;
	my $alignedSum;
	my $inputCountSum;
	push(@headers, "Sample");
    push(@headers, "Total Reads");
	push(@headers, "Multimapped Reads Aligned (HISAT2)");
	push(@headers, "Unique Reads Aligned (HISAT2)");
	foreach my $file (@files){
		my $multimapped;
		my $aligned;
		my $inputCount;
		chomp($inputCount = `cat $file | grep 'reads; of these:' | awk '{sum+=\\$1} END {print sum}'`);
		chomp($aligned = `cat $file | grep 'aligned.*exactly 1 time' | awk '{sum+=\\$1} END {print sum}'`);
		chomp($multimapped = `cat $file | grep 'aligned.*>1 times' | awk '{sum+=\\$1} END {print sum}'`);
		$multimappedSum += int($multimapped);
        $alignedSum += int($aligned);
        $inputCountSum += int($inputCount);
	}
	$tsv{$name} = [$name, $inputCountSum];
	push(@{$tsv{$name}}, $multimappedSum);
	push(@{$tsv{$name}}, $alignedSum);
}
'''

}


process HISAT2_Module_Merge_TSV_Files {

publishDir params.outdir, overwrite: true, mode: 'copy',
	saveAs: {filename ->
	if (filename =~ /${name}.tsv$/) "hisat2_summary/$filename"
}

input:
 file tsv from g216_2_outputFile_g216_10.collect()
 val outputFileName from g216_2_name_g216_10.collect()

output:
 file "${name}.tsv"  into g216_10_outputFileTSV_g_198

errorStrategy 'retry'
maxRetries 3

script:
name = outputFileName[0]
"""    
awk 'FNR==1 && NR!=1 {  getline; } 1 {print} ' *.tsv > ${name}.tsv
"""
}

mappingListQuoteSep = mapList.collect{ '"' + it + '"'}.join(",") 
rawIndexList = indexList.collect{ '"' + it + '"'}.join(",") 
process Sequential_Mapping_Module_Sequential_Mapping_Bam_dedup_count {

publishDir params.outdir, overwrite: true, mode: 'copy',
	saveAs: {filename ->
	if (filename =~ /.*.counts.tsv$/) "Sequential_Mapping_Bam_dedup_count/$filename"
}

input:
 file bam from g241_34_bam_file_g241_36.collect()
 file index from g241_34_bam_index_g241_36.collect()

output:
 file "*.counts.tsv"  into g241_36_outputFileTSV

shell:
'''
#!/usr/bin/env perl
use List::Util qw[min max];
use File::Basename;
use Getopt::Long;
use Pod::Usage;
use Data::Dumper;

my @header;
my @header_antisense;
my @header_sense;
my %all_files;
my %sense_files;
my %antisense_files;

my @mappingList = (!{mappingListQuoteSep});
my @rawIndexList = (!{rawIndexList});
my %indexHash;
my $dedup = "";
@indexHash{@mappingList} = @rawIndexList;

chomp(my $contents = `ls *.bam`);
my @files = split(/[\\n]+/, $contents);
foreach my $file (@files){
        $file=~/(.*)@(.*)_sorted(.*)\\.bam/;
        my $mapper = $1; 
        my $name = $2; ##header
        #print $3;
        if ($3 eq ".dedup"){
            $dedup = ".dedup";
        }
        if ($name=~/_antisense$/){
        	push(@header_antisense, $name) unless grep{$_ eq $name} @header_antisense; #mapped element header
        	$antisense_files{$mapper} .= $file." ";
        }
        elsif ($name=~/_sense$/){
        	push(@header_sense, $name) unless grep{$_ eq $name} @header_sense; #mapped element header
        	$sense_files{$mapper} .= $file." ";
        }
        else{
			push(@header, $name) unless grep{$_ eq $name} @header; #mapped element header
	        $all_files{$mapper} .= $file." ";
        }
}

runCov(\\%all_files, \\@header, \\@indexHash, "", $dedup);
runCov(\\%sense_files, \\@header_sense, \\@indexHash, "sense", $dedup);
runCov(\\%antisense_files, \\@header_antisense, \\@indexHash, "antisense", $dedup);

sub runCov {
	my ( \$files, \$header, \$indexHash, \$sense_antisense, \$dedup) = @_;
	open OUT, ">header".\$sense_antisense.".tsv";
	print OUT join ("\\t", "id","len",@{\$header}),"\\n";
	close OUT;
	my $par = "";
	if ($sense_antisense=~/^sense\$/){
      $par = "-s";
    }elsif($sense_antisense=~/^antisense\$/){
      $par = "-S";
    }
	
	foreach my $key (sort keys %{\$files}) {  
	   my $bamFiles = ${\$files}{$key};
		unless (-e ${indexHash}{$key}.".bed") {
            print "2: bed not found run makeBed\n";
                if (-e ${indexHash}{$key}.".fa") {
                    makeBed(${indexHash}{$key}.".fa", $key, ${indexHash}{$key}.".bed");
                } elsif(-e ${indexHash}{$key}.".fasta"){
                    makeBed(${indexHash}{$key}.".fasta", $key, ${indexHash}{$key}.".bed");
                }
        }
	    
		my $com =  "bedtools multicov $par -bams $bamFiles -bed ".${indexHash}{$key}.".bed >$key${dedup}${sense_antisense}.counts.tmp\n";
        print $com;
        `$com`;
        my $iniResColumn = int(countColumn(${indexHash}{$key}.".bed")) + 1;
	    `awk -F \\"\\\\t\\" \\'{a=\\"\\";for (i=$iniResColumn;i<=NF;i++){a=a\\"\\\\t\\"\\$i;} print \\$4\\"\\\\t\\"(\\$3-\\$2)\\"\\"a}\\' $key${dedup}${sense_antisense}.counts.tmp> $key${dedup}${sense_antisense}.counts.tsv`;
	    `sort -k3,3nr $key${dedup}${sense_antisense}.counts.tsv>$key${dedup}${sense_antisense}.sorted.tsv`;
        `cat header${sense_antisense}.tsv $key${dedup}${sense_antisense}.sorted.tsv> $key${dedup}${sense_antisense}.counts.tsv`;
	}
}

sub countColumn {
    my ( \$file) = @_;
    open(IN, \$file);
    my $line=<IN>;
    chomp($line);
    my @cols = split('\\t', $line);
    my $n = @cols;
    close OUT;
    return $n;
}

sub makeBed {
    my ( \$fasta, \$type, \$bed) = @_;
    print "makeBed $fasta\\n";
    print "makeBed $bed\\n";
    open OUT, ">$bed";
    open(IN, \$fasta);
    my $name="";
    my $seq="";
    my $i=0;
    while(my $line=<IN>){
        chomp($line);
        if($line=~/^>(.*)/){
            $i++ if (length($seq)>0);
            print OUT "$name\\t1\\t".length($seq)."\\t$name\\t0\\t+\\n" if (length($seq)>0); 
            $name="$1";
            $seq="";
        } elsif($line=~/[ACGTNacgtn]+/){
            $seq.=$line;
        }
    }
    print OUT "$name\\t1\\t".length($seq)."\\t$name\\t0\\t+\\n" if (length($seq)>0); 
    close OUT;
}

'''
}

mappingListQuoteSep = mapList.collect{ '"' + it + '"'}.join(",") 
rawIndexList = indexList.collect{ '"' + it + '"'}.join(",") 
process Sequential_Mapping_Module_Sequential_Mapping_Bam_count {

publishDir params.outdir, overwrite: true, mode: 'copy',
	saveAs: {filename ->
	if (filename =~ /.*.counts.tsv$/) "Sequential_Mapping_Bam_count/$filename"
}

input:
 file bam from g241_34_bam_file_g241_35.collect()
 file index from g241_34_bam_index_g241_35.collect()

output:
 file "*.counts.tsv"  into g241_35_outputFileTSV

shell:
'''
#!/usr/bin/env perl
use List::Util qw[min max];
use File::Basename;
use Getopt::Long;
use Pod::Usage;
use Data::Dumper;

my @header;
my @header_antisense;
my @header_sense;
my %all_files;
my %sense_files;
my %antisense_files;

my @mappingList = (!{mappingListQuoteSep});
my @rawIndexList = (!{rawIndexList});
my %indexHash;
my $dedup = "";
@indexHash{@mappingList} = @rawIndexList;

chomp(my $contents = `ls *.bam`);
my @files = split(/[\\n]+/, $contents);
foreach my $file (@files){
        $file=~/(.*)@(.*)_sorted(.*)\\.bam/;
        my $mapper = $1; 
        my $name = $2; ##header
        #print $3;
        if ($3 eq ".dedup"){
            $dedup = ".dedup";
        }
        if ($name=~/_antisense$/){
        	push(@header_antisense, $name) unless grep{$_ eq $name} @header_antisense; #mapped element header
        	$antisense_files{$mapper} .= $file." ";
        }
        elsif ($name=~/_sense$/){
        	push(@header_sense, $name) unless grep{$_ eq $name} @header_sense; #mapped element header
        	$sense_files{$mapper} .= $file." ";
        }
        else{
			push(@header, $name) unless grep{$_ eq $name} @header; #mapped element header
	        $all_files{$mapper} .= $file." ";
        }
}

runCov(\\%all_files, \\@header, \\@indexHash, "", $dedup);
runCov(\\%sense_files, \\@header_sense, \\@indexHash, "sense", $dedup);
runCov(\\%antisense_files, \\@header_antisense, \\@indexHash, "antisense", $dedup);

sub runCov {
	my ( \$files, \$header, \$indexHash, \$sense_antisense, \$dedup) = @_;
	open OUT, ">header".\$sense_antisense.".tsv";
	print OUT join ("\\t", "id","len",@{\$header}),"\\n";
	close OUT;
	my $par = "";
	if ($sense_antisense=~/^sense\$/){
      $par = "-s";
    }elsif($sense_antisense=~/^antisense\$/){
      $par = "-S";
    }
	
	foreach my $key (sort keys %{\$files}) {  
	   my $bamFiles = ${\$files}{$key};
		unless (-e ${indexHash}{$key}.".bed") {
            print "2: bed not found run makeBed\n";
                if (-e ${indexHash}{$key}.".fa") {
                    makeBed(${indexHash}{$key}.".fa", $key, ${indexHash}{$key}.".bed");
                } elsif(-e ${indexHash}{$key}.".fasta"){
                    makeBed(${indexHash}{$key}.".fasta", $key, ${indexHash}{$key}.".bed");
                }
        }
	    
		my $com =  "bedtools multicov $par -bams $bamFiles -bed ".${indexHash}{$key}.".bed >$key${dedup}${sense_antisense}.counts.tmp\n";
        print $com;
        `$com`;
        my $iniResColumn = int(countColumn(${indexHash}{$key}.".bed")) + 1;
	    `awk -F \\"\\\\t\\" \\'{a=\\"\\";for (i=$iniResColumn;i<=NF;i++){a=a\\"\\\\t\\"\\$i;} print \\$4\\"\\\\t\\"(\\$3-\\$2)\\"\\"a}\\' $key${dedup}${sense_antisense}.counts.tmp> $key${dedup}${sense_antisense}.counts.tsv`;
	    `sort -k3,3nr $key${dedup}${sense_antisense}.counts.tsv>$key${dedup}${sense_antisense}.sorted.tsv`;
        `cat header${sense_antisense}.tsv $key${dedup}${sense_antisense}.sorted.tsv> $key${dedup}${sense_antisense}.counts.tsv`;
	}
}

sub countColumn {
    my ( \$file) = @_;
    open(IN, \$file);
    my $line=<IN>;
    chomp($line);
    my @cols = split('\\t', $line);
    my $n = @cols;
    close OUT;
    return $n;
}

sub makeBed {
    my ( \$fasta, \$type, \$bed) = @_;
    print "makeBed $fasta\\n";
    print "makeBed $bed\\n";
    open OUT, ">$bed";
    open(IN, \$fasta);
    my $name="";
    my $seq="";
    my $i=0;
    while(my $line=<IN>){
        chomp($line);
        if($line=~/^>(.*)/){
            $i++ if (length($seq)>0);
            print OUT "$name\\t1\\t".length($seq)."\\t$name\\t0\\t+\\n" if (length($seq)>0); 
            $name="$1";
            $seq="";
        } elsif($line=~/[ACGTNacgtn]+/){
            $seq.=$line;
        }
    }
    print OUT "$name\\t1\\t".length($seq)."\\t$name\\t0\\t+\\n" if (length($seq)>0); 
    close OUT;
}

'''
}


process Sequential_Mapping_Module_Deduplication_Summary {

publishDir params.outdir, overwrite: true, mode: 'copy',
	saveAs: {filename ->
	if (filename =~ /deduplication_summary.tsv$/) "sequential_mapping_summary/$filename"
}

input:
 file flagstat from g241_34_log_file_g241_30.collect()
 val mate from g_229_mate_g241_30

output:
 file "deduplication_summary.tsv"  into g241_30_outputFileTSV

errorStrategy 'retry'
maxRetries 2

shell:
'''
#!/usr/bin/env perl
use List::Util qw[min max];
use strict;
use File::Basename;
use Getopt::Long;
use Pod::Usage;
use Data::Dumper;

my @header;
my %all_files;
my %tsv;
my %headerHash;
my %headerText;

my $i=0;
chomp(my $contents = `ls *_duplicates_stats.log`);
my @files = split(/[\\n]+/, $contents);
foreach my $file (@files){
    $i++;
    $file=~/(.*)@(.*)@(.*)_duplicates_stats\\.log/;
    my $mapOrder = int($1); 
    my $mapper = $2; #mapped element 
    my $name = $3; ##sample name
    push(@header, $mapper) unless grep{$_ eq $mapper} @header; 
        
    # my $duplicates;
    my $aligned;
    my $dedup; #aligned reads after dedup
    my $percent=0;
    if ("!{mate}" eq "pair" ){
        #first flagstat belongs to first bam file
        chomp($aligned = `cat $file | grep 'properly paired (' | sed -n 1p | awk '{sum+=\\$1+\\$3} END {print sum}'`);
        #second flagstat belongs to dedup bam file
        chomp($dedup = `cat $file | grep 'properly paired (' | sed -n 2p | awk '{sum+=\\$1+\\$3} END {print sum}'`);
    } else {
        chomp($aligned = `cat $file | grep 'mapped (' | sed -n 1p | awk '{sum+=\\$1+\\$3} END {print sum}'`);
        chomp($dedup = `cat $file | grep 'mapped (' | sed -n 2p | awk '{sum+=\\$1+\\$3} END {print sum}'`);
    }
    # chomp($duplicates = `cat $file | grep 'duplicates' | awk '{sum+=\\$1+\\$3} END {print sum}'`);
    # $dedup = int($aligned) - int($duplicates);
    if ("!{mate}" eq "pair" ){
       $dedup = int($dedup/2);
       $aligned = int($aligned/2);
    } 
    $percent = "0.00";
    if (int($aligned)  > 0 ){
       $percent = sprintf("%.2f", ($aligned-$dedup)/$aligned*100); 
    } 
    $tsv{$name}{$mapper}=[$aligned,$dedup,"$percent%"];
    $headerHash{$mapOrder}=$mapper;
    $headerText{$mapOrder}=["$mapper (Before Dedup)", "$mapper (After Dedup)", "$mapper (Duplication Ratio %)"];
}

my @mapOrderArray = ( keys %headerHash );
my @sortedOrderArray = sort { $a <=> $b } @mapOrderArray;

my $summary = "deduplication_summary.tsv";
open(OUT, ">$summary");
print OUT "Sample\\t";
my @headArr = ();
for my $mapOrder (@sortedOrderArray) {
    push (@headArr, @{$headerText{$mapOrder}});
}
my $headArrAll = join("\\t", @headArr);
print OUT "$headArrAll\\n";

foreach my $name (keys %tsv){
    my @rowArr = ();
    for my $mapOrder (@sortedOrderArray) {
        push (@rowArr, @{$tsv{$name}{$headerHash{$mapOrder}}});
    }
    my $rowArrAll = join("\\t", @rowArr);
    print OUT "$name\\t$rowArrAll\\n";
}
close(OUT);
'''
}

g227_11_genomeIndexPath_g242_22= g227_11_genomeIndexPath_g242_22.ifEmpty([""]) 

//* params.rsem_ref_using_star_index =  ""  //* @input
//* params.rsem_ref_using_bowtie2_index =  ""  //* @input
//* params.rsem_ref_using_bowtie_index =  ""  //* @input
//* @style @condition:{no_bam_output="false", output_genome_bam}, {no_bam_output="true"} @multicolumn:{no_bam_output, output_genome_bam}

//* autofill
if ($HOSTNAME == "default"){
    $CPU  = 3
    $MEMORY = 32
}
//* platform
if ($HOSTNAME == "ghpcc06.umassrc.org"){
    $TIME = 3000
    $CPU  = 4
    $MEMORY = 20
    $QUEUE = "long"
}
//* platform
//* autofill


process RSEM_module_RSEM {

publishDir params.outdir, overwrite: true, mode: 'copy',
	saveAs: {filename ->
	if (filename =~ /pipe.rsem.${name}$/) "rsem/$filename"
}

input:
 val mate from g_229_mate_g242_22
 set val(name), file(reads) from g241_34_reads_g242_22
 val rsemIndex from g227_11_genomeIndexPath_g242_22.collect()

output:
 file "pipe.rsem.${name}"  into g242_22_rsemOut_g242_17, g242_22_rsemOut_g242_21, g242_22_rsemOut_g_177
 set val(name), file("pipe.rsem.*/*.genome.bam") optional true  into g242_22_bam_file
 set val(name), file("pipe.rsem.*/*.bam") optional true  into g242_22_mapped_reads
 file "pipe.rsem.*/*.genome.bam" optional true  into g242_22_genome_bam_g242_23

errorStrategy 'retry'
maxRetries 1

when:
(params.run_RSEM && (params.run_RSEM == "yes")) || !params.run_RSEM

script:
RSEM_reference_type = params.RSEM_module_RSEM.RSEM_reference_type
RSEM_parameters = params.RSEM_module_RSEM.RSEM_parameters
no_bam_output = params.RSEM_module_RSEM.no_bam_output
output_genome_bam = params.RSEM_module_RSEM.output_genome_bam
sense_antisense = params.RSEM_module_RSEM.sense_antisense

nameAll = reads.toString()
nameArray = nameAll.split(' ')
def file2;

if (nameAll.contains('.gz')) {
    newName =  nameArray[0] - ~/(\.fastq.gz)?(\.fq.gz)?$/
    file1 =  nameArray[0] - '.gz' 
    if (mate == "pair") {file2 =  nameArray[1] - '.gz'}
    runGzip = "ls *.gz | xargs -i echo gzip -df {} | sh"
} else {
    newName =  nameArray[0] - ~/(\.fastq)?(\.fq)?$/
    file1 =  nameArray[0]
    if (mate == "pair") {file2 =  nameArray[1]}
    runGzip = ''
}
genome_BamText = (output_genome_bam.toString() != "false") ? "--output-genome-bam" : ""
if (no_bam_output.toString() != "false"){
    noBamText = "--no-bam-output"
    genome_BamText = ""
} else {
    noBamText = ""
}


rsemRef = ""
refType = ""
if (RSEM_reference_type == "star"){
    rsemRef = params.rsem_ref_using_star_index
    refType = "--star"
} else if (RSEM_reference_type == "bowtie2"){
    rsemRef = params.rsem_ref_using_bowtie2_index
    refType = "--bowtie2"
} else if (RSEM_reference_type == "bowtie"){
    rsemRef = params.rsem_ref_using_bowtie_index
    refType = ""
}
"""
$runGzip
mkdir -p pipe.rsem.${name}

if [ "${mate}" == "pair" ]; then
    echo "rsem-calculate-expression ${refType} ${RSEM_parameters} ${genome_BamText} ${noBamText} --paired-end ${file1} ${file2} ${rsemRef} pipe.rsem.${name}/rsem.out.${name}"
    rsem-calculate-expression ${refType} ${RSEM_parameters} ${genome_BamText} ${noBamText} --paired-end ${file1} ${file2} ${rsemRef} pipe.rsem.${name}/rsem.out.${name}
	if [ "${sense_antisense}" == "Yes" ]; then
		 rsem-calculate-expression ${refType} ${RSEM_parameters} ${genome_BamText} ${noBamText} --forward-prob 1 --paired-end ${file1} ${file2} ${rsemRef} pipe.rsem.${name}/rsem.out.forward.${name}
		 rsem-calculate-expression ${refType} ${RSEM_parameters} ${genome_BamText} ${noBamText} --forward-prob 0 --paired-end ${file1} ${file2} ${rsemRef} pipe.rsem.${name}/rsem.out.reverse.${name}
	fi
else
    echo "rsem-calculate-expression ${refType} ${RSEM_parameters} ${genome_BamText} ${noBamText} --calc-ci ${file1} ${rsemRef} pipe.rsem.${name}/rsem.out.${name}"
    rsem-calculate-expression ${refType} ${RSEM_parameters} ${genome_BamText} ${noBamText} --calc-ci ${file1} ${rsemRef} pipe.rsem.${name}/rsem.out.${name}
	if [ "${sense_antisense}" == "Yes" ]; then
	    rsem-calculate-expression ${refType} ${RSEM_parameters} ${genome_BamText} ${noBamText} --forward-prob 1 --calc-ci ${file1} ${rsemRef} pipe.rsem.${name}/rsem.out.forward.${name}
    	rsem-calculate-expression ${refType} ${RSEM_parameters} ${genome_BamText} ${noBamText} --forward-prob 0 --calc-ci ${file1} ${rsemRef} pipe.rsem.${name}/rsem.out.reverse.${name}
	fi
fi
"""

}


process RSEM_module_file_to_set_conversion_for_bam {

input:
 file bam from g242_22_genome_bam_g242_23.flatten()

output:
 set val(name),file(bam)  into g242_23_bam_file_g219_121, g242_23_bam_file_g219_122, g242_23_bam_file_g219_123, g242_23_bam_file_g219_124, g242_23_bam_file_g219_126

script:
name = bam.baseName
"""
echo "done"	
"""
}

//* params.gtf =  ""  //* @input


process BAM_Analysis_RSEM_featureCounts {

publishDir params.outdir, overwrite: true, mode: 'copy',
	saveAs: {filename ->
	if (filename =~ /.*$/) "featureCounts_after_rsem/$filename"
}

input:
 set val(name), file(bam) from g242_23_bam_file_g219_126
 val paired from g_229_mate_g219_126
 each run_params from g219_125_run_parameters_g219_126

output:
 file "*"  into g219_126_outputFileTSV_g219_117

script:
pairText = ""
if (paired == "pair"){
    pairText = "-p"
}

run_name = run_params["run_name"] 
run_parameters = run_params["run_parameters"] 

"""
featureCounts ${pairText} ${run_parameters} -a ${params.gtf} -o ${name}@${run_name}@fCounts.txt ${bam}
## remove first line
sed -i '1d' ${name}@${run_name}@fCounts.txt

"""
}

//* autofill
//* platform
if ($HOSTNAME == "ghpcc06.umassrc.org"){
    $TIME = 30
    $CPU  = 1
    $MEMORY = 10
    $QUEUE = "short"
}
//* platform
//* autofill

process BAM_Analysis_RSEM_summary_featureCounts {

publishDir params.outdir, overwrite: true, mode: 'copy',
	saveAs: {filename ->
	if (filename =~ /.*_featureCounts.tsv$/) "featureCounts_after_rsem_summary/$filename"
	else if (filename =~ /.*_featureCounts.sum.tsv$/) "featureCounts_after_rsem_details/$filename"
}

input:
 file featureCountsOut from g219_126_outputFileTSV_g219_117.collect()

output:
 file "*_featureCounts.tsv"  into g219_117_outputFile
 file "*_featureCounts.sum.tsv"  into g219_117_outFileTSV

shell:
'''
#!/usr/bin/env perl

# Step 1: Merge count files
my %tf = ( expected_count => 6 );
my @run_name=();
chomp(my $contents = `ls *@fCounts.txt`);
my @files = split(/[\\n]+/, $contents);
foreach my $file (@files){
    $file=~/(.*)\\@(.*)\\@fCounts\\.txt/;
    my $runname = $2;
    push(@run_name, $runname) unless grep{$_ eq $runname} @run_name;
}


my @expectedCount_ar = ("expected_count");
for($l = 0; $l <= $#run_name; $l++) {
    my $runName = $run_name[$l];
    for($ll = 0; $ll <= $#expectedCount_ar; $ll++) {
        my $expectedCount = $expectedCount_ar[$ll];
    
        my @a=();
        my %b=();
        my %c=();
        my $i=0;
        chomp(my $contents = `ls *\\@${runName}\\@fCounts.txt`);
        my @files = split(/[\\n]+/, $contents);
        foreach my $file (@files){
        $i++;
        $file=~/(.*)\\@${runName}\\@fCounts\\.txt/;
        my $libname = $1; 
        $a[$i]=$libname;
        open IN, $file;
            $_=<IN>;
            while(<IN>){
                my @v=split; 
                $b{$v[0]}{$i}=$v[$tf{$expectedCount}];
                $c{$v[0]}=$v[5]; #length column
            }
            close IN;
        }
        my $outfile="$runName"."_featureCounts.tsv";
        open OUT, ">$outfile";
        if ($runName eq "transcript_id") {
            print OUT "transcript\tlength";
        } else {
            print OUT "gene\tlength";
        }
    
        for(my $j=1;$j<=$i;$j++) {
            print OUT "\t$a[$j]";
        }
        print OUT "\n";
    
        foreach my $key (keys %b) {
            print OUT "$key\t$c{$key}";
            for(my $j=1;$j<=$i;$j++){
                print OUT "\t$b{$key}{$j}";
            }
            print OUT "\n";
        }
        close OUT;
         
    }
}


	

# Step 2: Merge summary files
for($l = 0; $l <= $#run_name; $l++) {
    my $runName = $run_name[$l];
    my @a=();
    my %b=();
    my $i=0;
    chomp(my $contents = `ls *\\@${runName}\\@fCounts.txt.summary`);
    my @files = split(/[\\n]+/, $contents);
    foreach my $file (@files){
        $i++;
        $file=~/(.*)\\@${runName}\\@fCounts\\.txt\\.summary/;
        my $libname = $1; 
        $a[$i]=$libname;
        open IN, $file;
        $_=<IN>;
        while(<IN>){
            my @v=split; 
            $b{$v[0]}{$i}=$v[1];
        }
        close IN;
    }
    my $outfile="$runName"."_featureCounts.sum.tsv";
    open OUT, ">$outfile";
    print OUT "criteria";
    for(my $j=1;$j<=$i;$j++) {
        print OUT "\t$a[$j]";
    }
    print OUT "\n";
    
    foreach my $key (keys %b) {
        print OUT "$key";
        for(my $j=1;$j<=$i;$j++){
            print OUT "\t$b{$key}{$j}";
        }
        print OUT "\n";
    }
    close OUT;
}

'''
}

//* params.genome_sizes =  ""  //* @input

//* autofill
if ($HOSTNAME == "default"){
    $CPU  = 1
    $MEMORY = 48
}
//* platform
if ($HOSTNAME == "ghpcc06.umassrc.org"){
    $TIME = 240
    $CPU  = 1
    $MEMORY = 48
    $QUEUE = "short"
} 
//* platform
//* autofill

process BAM_Analysis_RSEM_UCSC_BAM2BigWig_converter {

publishDir params.outdir, overwrite: true, mode: 'copy',
	saveAs: {filename ->
	if (filename =~ /.*.bw$/) "bigwig_rsem/$filename"
}

input:
 set val(name), file(bam) from g242_23_bam_file_g219_124

output:
 file "*.bw"  into g219_124_outputFileBw

when:
(params.run_BigWig_Conversion && (params.run_BigWig_Conversion == "yes")) || !params.run_BigWig_Conversion

script:
nameAll = bam.toString()
if (nameAll.contains('_sorted.bam')) {
    runSamtools = "samtools index ${nameAll}"
    nameFinal = nameAll
} else {
    runSamtools = "samtools sort -o ${name}_sorted.bam $bam && samtools index ${name}_sorted.bam "
    nameFinal = "${name}_sorted.bam"
}

"""
$runSamtools
bedtools genomecov -split -bg -ibam ${nameFinal} -g ${params.genome_sizes} > ${name}.bg 
wigToBigWig -clip -itemsPerSlot=1 ${name}.bg ${params.genome_sizes} ${name}.bw 
"""
}

igv_extention_factor = params.BAM_Analysis_RSEM_IGV_BAM2TDF_converter.igv_extention_factor
igv_window_size = params.BAM_Analysis_RSEM_IGV_BAM2TDF_converter.igv_window_size

//* params.genome =  ""  //* @input

//* autofill
if ($HOSTNAME == "default"){
    $CPU  = 1
    $MEMORY = 32
}
//* platform
if ($HOSTNAME == "ghpcc06.umassrc.org"){
    $TIME = 400
    $CPU  = 1
    $MEMORY = 32
    $QUEUE = "long"
} 
//* platform
//* autofill

process BAM_Analysis_RSEM_IGV_BAM2TDF_converter {

publishDir params.outdir, overwrite: true, mode: 'copy',
	saveAs: {filename ->
	if (filename =~ /.*.tdf$/) "igvtools_rsem/$filename"
}

input:
 val mate from g_229_mate_g219_123
 set val(name), file(bam) from g242_23_bam_file_g219_123

output:
 file "*.tdf"  into g219_123_outputFileOut

when:
(params.run_IGV_TDF_Conversion && (params.run_IGV_TDF_Conversion == "yes")) || !params.run_IGV_TDF_Conversion

script:
pairedText = (params.nucleicAcidType == "dna" && mate == "pair") ? " --pairs " : ""
nameAll = bam.toString()
if (nameAll.contains('_sorted.bam')) {
    runSamtools = "samtools index ${nameAll}"
    nameFinal = nameAll
} else {
    runSamtools = "samtools sort -o ${name}_sorted.bam $bam && samtools index ${name}_sorted.bam "
    nameFinal = "${name}_sorted.bam"
}
"""
$runSamtools
igvtools count -w ${igv_window_size} -e ${igv_extention_factor} ${pairedText} ${nameFinal} ${name}.tdf ${params.genome_sizes}
"""
}

//* params.bed =  ""  //* @input

process BAM_Analysis_RSEM_RSeQC {

publishDir params.outdir, overwrite: true, mode: 'copy',
	saveAs: {filename ->
	if (filename =~ /result\/.*.out$/) "rseqc_rsem/$filename"
}

input:
 set val(name), file(bam) from g242_23_bam_file_g219_122

output:
 file "result/*.out"  into g219_122_outputFileOut_g219_95, g219_122_outputFileOut_g_177

when:
(params.run_RSeQC && (params.run_RSeQC == "yes")) || !params.run_RSeQC

script:
"""
mkdir result
read_distribution.py  -i ${bam} -r ${params.bed}> result/RSeQC.${name}.out
"""
}


process BAM_Analysis_RSEM_RSeQC_Summary {

publishDir params.outdir, overwrite: true, mode: 'copy',
	saveAs: {filename ->
	if (filename =~ /.*.tsv$/) "rseqc_summary_rsem/$filename"
}

input:
 file rseqcOut from g219_122_outputFileOut_g219_95.collect()
 val mate from g_229_mate_g219_95

output:
 file "*.tsv"  into g219_95_outputFileTSV

shell:
'''
#!/usr/bin/env perl
use List::Util qw[min max];
use strict;
use File::Basename;
use Getopt::Long;
use Pod::Usage; 
use Data::Dumper;

my $indir = $ENV{'PWD'};
my $outd = $ENV{'PWD'};
my @files = ();
my @outtypes = ("RSeQC");
my @order=( "Total Reads", "Total Tags" , "Total Assigned Tags", "CDS_Exons", "5'UTR_Exons", "3'UTR_Exons", "Introns", "TSS_up_1kb", "TSS_up_5kb", "TSS_up_10kb", "TES_down_1kb", "TES_down_5kb", "TES_down_10kb");
my %lines=(
  "Total Reads" => 1,
  "Total Tags" => 1,
  "Total Assigned Tags" => 1,
  "CDS_Exons" => 2,
  "5'UTR_Exons" => 2,
  "3'UTR_Exons" => 2,
  "Introns" => 2,
  "TSS_up_1kb" => 2,
  "TSS_up_5kb" => 2,
  "TSS_up_10kb" => 2,
  "TES_down_1kb" => 2,
  "TES_down_5kb" => 2,
  "TES_down_10kb" => 2
);


foreach my $outtype (@outtypes)
{

my $ext=".out";
@files = <$indir/$outtype*$ext>;

my @rowheaders=();
my @libs=();
my %vals=();
my %normvals=();
my $type = "rsem";

foreach my $d (@files){
  my $libname=basename($d, $ext);
  $libname=~s/RSeQC.//g;
  $libname=~s/rsem.out.//g;
  $libname=~s/.genome//g;
  print $libname."\\n";
  push(@libs, $libname); 
  getVals($d, $libname, \\%vals, \\%normvals, \\%lines);
}
#print Dumper(%vals);
#print Dumper(%normvals);

my $sizemetrics = keys %vals;
write_results("$outd/$outtype.$type.counts.tsv", \\@libs,\\%vals, \\@order, "region") if ($sizemetrics>0);
write_results("$outd/$outtype.$type.tagskb.tsv", \\@libs,\\%normvals, \\@order, "region") if ($sizemetrics>0);

}

sub write_results
{
  my ($outfile, $libs, $vals, $order, $name )=@_;
  open(OUT, ">$outfile");
  print OUT "$name\\t".join("\\t", @{$libs})."\\n";

  my $lib=${$libs}[0];
  foreach my $key ( @order )
  {
    if (exists ${$vals}{$lib}{$key}) {
    print OUT $key;
    foreach my $lib (@{$libs})
    {
      print OUT "\\t".${$vals}{$lib}{$key};
    } 
    print OUT "\\n";
    }
  }
  close(OUT);
}

sub getVals{
  my ($filename, $libname, $vals, $normvals, $lines)=@_;
  if (-e $filename){
     open(IN, $filename);
     while(my $line=<IN>)
     {
       chomp($line);
       my @vals_arr=split(/\\s{2,}/,$line);
       if (exists ${$lines}{$vals_arr[0]}) {
         my $idx=${$lines}{$vals_arr[0]};
         ${$vals}{$libname}{$vals_arr[0]}=$vals_arr[$idx] if (exists $vals_arr[$idx]);
         if ($idx==2) {
             ${$normvals}{$libname}{$vals_arr[0]}=$vals_arr[3] if (exists $vals_arr[3]);
         }
       }
     } 
  }
  
}
'''

}

//* params.pdfbox_path =  ""  //* @input
//* autofill
if ($HOSTNAME == "default"){
    $CPU  = 1
    $MEMORY = 32
}
//* platform
if ($HOSTNAME == "ghpcc06.umassrc.org"){
    $TIME = 240
    $CPU  = 1
    $MEMORY = 32
    $QUEUE = "short"
}
//* platform
//* autofill

process BAM_Analysis_RSEM_Picard {

input:
 set val(name), file(bam) from g242_23_bam_file_g219_121

output:
 file "*_metrics"  into g219_121_outputFileOut_g219_82
 file "results/*.pdf"  into g219_121_outputFilePdf_g219_82

when:
(params.run_Picard_CollectMultipleMetrics && (params.run_Picard_CollectMultipleMetrics == "yes")) || !params.run_Picard_CollectMultipleMetrics

script:
"""
picard CollectMultipleMetrics OUTPUT=${name}_multiple.out VALIDATION_STRINGENCY=LENIENT INPUT=${bam}
mkdir results && java -jar ${params.pdfbox_path} PDFMerger *.pdf results/${name}_multi_metrics.pdf
"""
}


process BAM_Analysis_RSEM_Picard_Summary {

publishDir params.outdir, overwrite: true, mode: 'copy',
	saveAs: {filename ->
	if (filename =~ /.*.tsv$/) "picard_summary_rsem/$filename"
	else if (filename =~ /results\/.*.pdf$/) "picard_summary_pdf_rsem/$filename"
}

input:
 file picardOut from g219_121_outputFileOut_g219_82.collect()
 val mate from g_229_mate_g219_82
 file picardPdf from g219_121_outputFilePdf_g219_82.collect()

output:
 file "*.tsv"  into g219_82_outputFileTSV
 file "results/*.pdf"  into g219_82_outputFilePdf

shell:
'''
#!/usr/bin/env perl
use List::Util qw[min max];
use strict;
use File::Basename;
use Getopt::Long;
use Pod::Usage; 
use Data::Dumper;

system("mkdir results && mv *.pdf results/. ");

my $indir = $ENV{'PWD'};
my $outd = $ENV{'PWD'};
my @files = ();
my @outtypes = ("CollectRnaSeqMetrics", "alignment_summary_metrics", "base_distribution_by_cycle_metrics", "insert_size_metrics", "quality_by_cycle_metrics", "quality_distribution_metrics" );

foreach my $outtype (@outtypes)
{
my $ext="_multiple.out";
$ext.=".$outtype" if ($outtype ne "CollectRnaSeqMetrics");
@files = <$indir/*$ext>;

my @rowheaders=();
my @libs=();
my %metricvals=();
my %histvals=();

my $pdffile="";
my $libname="";
foreach my $d (@files){
  my $libname=basename($d, $ext);
  print $libname."\\n";
  push(@libs, $libname); 
  getMetricVals($d, $libname, \\%metricvals, \\%histvals, \\@rowheaders);
}

my $sizemetrics = keys %metricvals;
write_results("$outd/$outtype.stats.tsv", \\@libs,\\%metricvals, \\@rowheaders, "metric") if ($sizemetrics>0);
my $sizehist = keys %histvals;
write_results("$outd/$outtype.hist.tsv", \\@libs,\\%histvals, "none", "nt") if ($sizehist>0);

}

sub write_results
{
  my ($outfile, $libs, $vals, $rowheaders, $name )=@_;
  open(OUT, ">$outfile");
  print OUT "$name\\t".join("\\t", @{$libs})."\\n";
  my $size=0;
  $size=scalar(@{${$vals}{${$libs}[0]}}) if(exists ${$libs}[0] and exists ${$vals}{${$libs}[0]} );
  
  for (my $i=0; $i<$size;$i++)
  { 
    my $rowname=$i;
    $rowname = ${$rowheaders}[$i] if ($name=~/metric/);
    print OUT $rowname;
    foreach my $lib (@{$libs})
    {
      print OUT "\\t".${${$vals}{$lib}}[$i];
    } 
    print OUT "\\n";
  }
  close(OUT);
}

sub getMetricVals{
  my ($filename, $libname, $metricvals, $histvals,$rowheaders)=@_;
  if (-e $filename){
     my $nextisheader=0;
     my $nextisvals=0;
     my $nexthist=0;
     open(IN, $filename);
     while(my $line=<IN>)
     {
       chomp($line);
       @{$rowheaders}=split(/\\t/, $line) if ($nextisheader && !scalar(@{$rowheaders})); 
       if ($nextisvals) {
         @{${$metricvals}{$libname}}=split(/\\t/, $line);
         $nextisvals=0;
       }
       if($nexthist){
          my @vals=split(/[\\s\\t]+/,$line); 
          push(@{${$histvals}{$libname}}, $vals[1]) if (exists $vals[1]);
       }
       $nextisvals=1 if ($nextisheader); $nextisheader=0;
       $nextisheader=1 if ($line=~/METRICS CLASS/);
       $nexthist=1 if ($line=~/normalized_position/);
     } 
  }
  
}
'''

}

//* autofill
//* platform
if ($HOSTNAME == "ghpcc06.umassrc.org"){
    $TIME = 30
    $CPU  = 1
    $MEMORY = 10
    $QUEUE = "short"
}
//* platform
//* autofill

process RSEM_module_RSEM_Alignment_Summary {

input:
 file rsemDir from g242_22_rsemOut_g242_17.collect()

output:
 file "rsem_alignment_sum.tsv"  into g242_17_outputFileTSV_g_198

shell:
'''

#!/usr/bin/env perl
use List::Util qw[min max];
use strict;
use File::Basename;
use Getopt::Long;
use Pod::Usage;
use Data::Dumper;
my $indir = $ENV{'PWD'};

opendir D, $indir or die "Could not open $indir";
my @alndirs = sort { $a cmp $b } grep /^pipe/, readdir(D);
closedir D;

my @a=();
my %b=();
my %c=();
my $i=0;
my @headers = ();
my %tsv;
foreach my $d (@alndirs){
    my $dir = "${indir}/$d";
    my $libname=$d;
    $libname=~s/pipe\\.rsem\\.//;
    my $multimapped;
    my $aligned;
    my $total;
    
    chomp($total = `awk 'NR == 1 {print \\$4}' ${dir}/rsem.out.$libname.stat/rsem.out.$libname.cnt`);
    chomp($aligned = `awk 'NR == 1 {print \\$2}' ${dir}/rsem.out.$libname.stat/rsem.out.$libname.cnt`);
    chomp($multimapped = `awk 'NR == 2 {print \\$2}' ${dir}/rsem.out.$libname.stat/rsem.out.$libname.cnt`);
    $tsv{$libname}=[$libname, $total];
    push(@{$tsv{$libname}}, $multimapped);
    push(@{$tsv{$libname}}, (int($aligned) - int($multimapped))."");
}


push(@headers, "Sample");
push(@headers, "Total Reads");
push(@headers, "Multimapped Reads Aligned (RSEM)");
push(@headers, "Unique Aligned Reads (RSEM)");


my @keys = keys %tsv;
my $summary = "rsem_alignment_sum.tsv";
my $header_string = join("\\t", @headers);
`echo "$header_string" > $summary`;
foreach my $key (@keys){
    my $values = join("\\t", @{ $tsv{$key} });
        `echo "$values" >> $summary`;
}
'''
}

g227_19_genomeIndexPath_g240_21= g227_19_genomeIndexPath_g240_21.ifEmpty([""]) 

//* params.kallisto_index =  ""  //* @input
//* params.genome_sizes =  ""  //* @input
//* params.gtf =  ""  //* @input
//* @style @multicolumn:{fragment_length,standard_deviation} @condition:{single_or_paired_end_reads="single", fragment_length,standard_deviation}, {single_or_paired_end_reads="pair"}


//* autofill
if ($HOSTNAME == "default"){
    $CPU  = 4
    $MEMORY = 20 
}
//* platform
if ($HOSTNAME == "ghpcc06.umassrc.org"){
    $TIME = 240
    $CPU  = 4
    $MEMORY = 15
    $QUEUE = "long"
}
//* platform
//* autofill


process Kallisto_module_kallisto_quant {

publishDir params.outdir, overwrite: true, mode: 'copy',
	saveAs: {filename ->
	if (filename =~ /kallisto_${name}$/) "kallisto/$filename"
}

input:
 val mate from g_229_mate_g240_21
 set val(name), file(reads) from g241_34_reads_g240_21
 val kallisto_index_path from g227_19_genomeIndexPath_g240_21

output:
 file "kallisto_${name}"  into g240_21_outputDir_g240_22, g240_21_outputDir_g240_28, g240_21_outputDir_g_177
 set val(name), file("kallisto_${name}/*.bam") optional true  into g240_21_bam_file_g224_121, g240_21_bam_file_g224_122, g240_21_bam_file_g224_123, g240_21_bam_file_g224_124, g240_21_bam_file_g224_126

when:
(params.run_Kallisto && (params.run_Kallisto == "yes")) || !params.run_Kallisto

script:
nameAll = reads.toString()
nameArray = nameAll.split(' ')
def file2;

if (nameAll.contains('.gz')) {
    newName =  nameArray[0] - ~/(\.fastq.gz)?(\.fq.gz)?$/
    file1 =  nameArray[0] - '.gz' 
    if (mate == "pair") {file2 =  nameArray[1] - '.gz'}
    runGzip = "ls *.gz | xargs -i echo gzip -df {} | sh"
} else {
    newName =  nameArray[0] - ~/(\.fastq)?(\.fq)?$/
    file1 =  nameArray[0]
    if (mate == "pair") {file2 =  nameArray[1]}
    runGzip = ''
}
single_or_paired_end_reads = params.Kallisto_module_kallisto_quant.single_or_paired_end_reads
fragment_length = params.Kallisto_module_kallisto_quant.fragment_length
standard_deviation = params.Kallisto_module_kallisto_quant.standard_deviation

kallisto_parameters = params.Kallisto_module_kallisto_quant.kallisto_parameters
genomebam = params.Kallisto_module_kallisto_quant.genomebam

genomebamText = (genomebam.toString() != "false") ? "--genomebam --gtf genes.gtf --chromosomes ${params.genome_sizes}" : ""
fragment_lengthText = (fragment_length.toString() != "" && single_or_paired_end_reads.toString() == "single") ? "-l ${fragment_length}" : ""
standard_deviationText = (standard_deviation.toString() != "" && single_or_paired_end_reads.toString() == "single") ? "-s ${standard_deviation}" : ""
"""
$runGzip
if [[ \$(awk '{print \$3}' ${params.gtf} | grep -c transcript) -le 1 ]]; then
    echo "transcript entries are not found in gtf file. gffread will add transcript entries."
    gffread -E --keep-genes ${params.gtf} -T -o- >genes.gtf 2>gffread.log
else
    ln -s ${params.gtf} genes.gtf
fi


mkdir -p kallisto_${name}
if [ "${mate}" == "pair" ]; then
    kallisto quant ${kallisto_parameters} -i ${params.kallisto_index} ${genomebamText} -o kallisto_${name} ${file1} ${file2} > kallisto_${name}/kallisto.log 2>&1
else
    kallisto quant --single ${kallisto_parameters} ${fragment_lengthText} ${standard_deviationText}  -i ${params.kallisto_index} ${genomebamText} -o kallisto_${name} ${file1} > kallisto_${name}/kallisto.log 2>&1
fi

if [ -f kallisto_${name}/pseudoalignments.bam ]; then
   mv kallisto_${name}/pseudoalignments.bam  kallisto_${name}/${name}.bam
fi
if [ -f kallisto_${name}/pseudoalignments.bam.bai ]; then
   mv kallisto_${name}/pseudoalignments.bam.bai  kallisto_${name}/${name}.bam.bai
fi


"""

}

//* params.gtf =  ""  //* @input


process BAM_Analysis_Module_Kallisto_featureCounts {

publishDir params.outdir, overwrite: true, mode: 'copy',
	saveAs: {filename ->
	if (filename =~ /.*$/) "featureCounts_after_kallisto/$filename"
}

input:
 set val(name), file(bam) from g240_21_bam_file_g224_126
 val paired from g_229_mate_g224_126
 each run_params from g224_125_run_parameters_g224_126

output:
 file "*"  into g224_126_outputFileTSV_g224_117

script:
pairText = ""
if (paired == "pair"){
    pairText = "-p"
}

run_name = run_params["run_name"] 
run_parameters = run_params["run_parameters"] 

"""
featureCounts ${pairText} ${run_parameters} -a ${params.gtf} -o ${name}@${run_name}@fCounts.txt ${bam}
## remove first line
sed -i '1d' ${name}@${run_name}@fCounts.txt

"""
}

//* autofill
//* platform
if ($HOSTNAME == "ghpcc06.umassrc.org"){
    $TIME = 30
    $CPU  = 1
    $MEMORY = 10
    $QUEUE = "short"
}
//* platform
//* autofill

process BAM_Analysis_Module_Kallisto_summary_featureCounts {

publishDir params.outdir, overwrite: true, mode: 'copy',
	saveAs: {filename ->
	if (filename =~ /.*_featureCounts.tsv$/) "featureCounts_after_Kallisto_summary/$filename"
	else if (filename =~ /.*_featureCounts.sum.tsv$/) "featureCounts_after_Kallisto_details/$filename"
}

input:
 file featureCountsOut from g224_126_outputFileTSV_g224_117.collect()

output:
 file "*_featureCounts.tsv"  into g224_117_outputFile
 file "*_featureCounts.sum.tsv"  into g224_117_outFileTSV

shell:
'''
#!/usr/bin/env perl

# Step 1: Merge count files
my %tf = ( expected_count => 6 );
my @run_name=();
chomp(my $contents = `ls *@fCounts.txt`);
my @files = split(/[\\n]+/, $contents);
foreach my $file (@files){
    $file=~/(.*)\\@(.*)\\@fCounts\\.txt/;
    my $runname = $2;
    push(@run_name, $runname) unless grep{$_ eq $runname} @run_name;
}


my @expectedCount_ar = ("expected_count");
for($l = 0; $l <= $#run_name; $l++) {
    my $runName = $run_name[$l];
    for($ll = 0; $ll <= $#expectedCount_ar; $ll++) {
        my $expectedCount = $expectedCount_ar[$ll];
    
        my @a=();
        my %b=();
        my %c=();
        my $i=0;
        chomp(my $contents = `ls *\\@${runName}\\@fCounts.txt`);
        my @files = split(/[\\n]+/, $contents);
        foreach my $file (@files){
        $i++;
        $file=~/(.*)\\@${runName}\\@fCounts\\.txt/;
        my $libname = $1; 
        $a[$i]=$libname;
        open IN, $file;
            $_=<IN>;
            while(<IN>){
                my @v=split; 
                $b{$v[0]}{$i}=$v[$tf{$expectedCount}];
                $c{$v[0]}=$v[5]; #length column
            }
            close IN;
        }
        my $outfile="$runName"."_featureCounts.tsv";
        open OUT, ">$outfile";
        if ($runName eq "transcript_id") {
            print OUT "transcript\tlength";
        } else {
            print OUT "gene\tlength";
        }
    
        for(my $j=1;$j<=$i;$j++) {
            print OUT "\t$a[$j]";
        }
        print OUT "\n";
    
        foreach my $key (keys %b) {
            print OUT "$key\t$c{$key}";
            for(my $j=1;$j<=$i;$j++){
                print OUT "\t$b{$key}{$j}";
            }
            print OUT "\n";
        }
        close OUT;
         
    }
}


	

# Step 2: Merge summary files
for($l = 0; $l <= $#run_name; $l++) {
    my $runName = $run_name[$l];
    my @a=();
    my %b=();
    my $i=0;
    chomp(my $contents = `ls *\\@${runName}\\@fCounts.txt.summary`);
    my @files = split(/[\\n]+/, $contents);
    foreach my $file (@files){
        $i++;
        $file=~/(.*)\\@${runName}\\@fCounts\\.txt\\.summary/;
        my $libname = $1; 
        $a[$i]=$libname;
        open IN, $file;
        $_=<IN>;
        while(<IN>){
            my @v=split; 
            $b{$v[0]}{$i}=$v[1];
        }
        close IN;
    }
    my $outfile="$runName"."_featureCounts.sum.tsv";
    open OUT, ">$outfile";
    print OUT "criteria";
    for(my $j=1;$j<=$i;$j++) {
        print OUT "\t$a[$j]";
    }
    print OUT "\n";
    
    foreach my $key (keys %b) {
        print OUT "$key";
        for(my $j=1;$j<=$i;$j++){
            print OUT "\t$b{$key}{$j}";
        }
        print OUT "\n";
    }
    close OUT;
}

'''
}

//* params.genome_sizes =  ""  //* @input

//* autofill
if ($HOSTNAME == "default"){
    $CPU  = 1
    $MEMORY = 48
}
//* platform
if ($HOSTNAME == "ghpcc06.umassrc.org"){
    $TIME = 240
    $CPU  = 1
    $MEMORY = 48
    $QUEUE = "short"
} 
//* platform
//* autofill

process BAM_Analysis_Module_Kallisto_UCSC_BAM2BigWig_converter {

publishDir params.outdir, overwrite: true, mode: 'copy',
	saveAs: {filename ->
	if (filename =~ /.*.bw$/) "bigwig_kallisto/$filename"
}

input:
 set val(name), file(bam) from g240_21_bam_file_g224_124

output:
 file "*.bw"  into g224_124_outputFileBw

when:
(params.run_BigWig_Conversion && (params.run_BigWig_Conversion == "yes")) || !params.run_BigWig_Conversion

script:
nameAll = bam.toString()
if (nameAll.contains('_sorted.bam')) {
    runSamtools = "samtools index ${nameAll}"
    nameFinal = nameAll
} else {
    runSamtools = "samtools sort -o ${name}_sorted.bam $bam && samtools index ${name}_sorted.bam "
    nameFinal = "${name}_sorted.bam"
}

"""
$runSamtools
bedtools genomecov -split -bg -ibam ${nameFinal} -g ${params.genome_sizes} > ${name}.bg 
wigToBigWig -clip -itemsPerSlot=1 ${name}.bg ${params.genome_sizes} ${name}.bw 
"""
}

igv_extention_factor = params.BAM_Analysis_Module_Kallisto_IGV_BAM2TDF_converter.igv_extention_factor
igv_window_size = params.BAM_Analysis_Module_Kallisto_IGV_BAM2TDF_converter.igv_window_size

//* params.genome =  ""  //* @input

//* autofill
if ($HOSTNAME == "default"){
    $CPU  = 1
    $MEMORY = 32
}
//* platform
if ($HOSTNAME == "ghpcc06.umassrc.org"){
    $TIME = 400
    $CPU  = 1
    $MEMORY = 32
    $QUEUE = "long"
} 
//* platform
//* autofill

process BAM_Analysis_Module_Kallisto_IGV_BAM2TDF_converter {

publishDir params.outdir, overwrite: true, mode: 'copy',
	saveAs: {filename ->
	if (filename =~ /.*.tdf$/) "igvtools_kallisto/$filename"
}

input:
 val mate from g_229_mate_g224_123
 set val(name), file(bam) from g240_21_bam_file_g224_123

output:
 file "*.tdf"  into g224_123_outputFileOut

when:
(params.run_IGV_TDF_Conversion && (params.run_IGV_TDF_Conversion == "yes")) || !params.run_IGV_TDF_Conversion

script:
pairedText = (params.nucleicAcidType == "dna" && mate == "pair") ? " --pairs " : ""
nameAll = bam.toString()
if (nameAll.contains('_sorted.bam')) {
    runSamtools = "samtools index ${nameAll}"
    nameFinal = nameAll
} else {
    runSamtools = "samtools sort -o ${name}_sorted.bam $bam && samtools index ${name}_sorted.bam "
    nameFinal = "${name}_sorted.bam"
}
"""
$runSamtools
igvtools count -w ${igv_window_size} -e ${igv_extention_factor} ${pairedText} ${nameFinal} ${name}.tdf ${params.genome_sizes}
"""
}

//* params.bed =  ""  //* @input

process BAM_Analysis_Module_Kallisto_RSeQC {

publishDir params.outdir, overwrite: true, mode: 'copy',
	saveAs: {filename ->
	if (filename =~ /result\/.*.out$/) "rseqc_kallisto/$filename"
}

input:
 set val(name), file(bam) from g240_21_bam_file_g224_122

output:
 file "result/*.out"  into g224_122_outputFileOut_g224_95, g224_122_outputFileOut_g_177

when:
(params.run_RSeQC && (params.run_RSeQC == "yes")) || !params.run_RSeQC

script:
"""
mkdir result
read_distribution.py  -i ${bam} -r ${params.bed}> result/RSeQC.${name}.out
"""
}


process BAM_Analysis_Module_Kallisto_RSeQC_Summary {

publishDir params.outdir, overwrite: true, mode: 'copy',
	saveAs: {filename ->
	if (filename =~ /.*.tsv$/) "rseqc_summary_kallisto/$filename"
}

input:
 file rseqcOut from g224_122_outputFileOut_g224_95.collect()
 val mate from g_229_mate_g224_95

output:
 file "*.tsv"  into g224_95_outputFileTSV

shell:
'''
#!/usr/bin/env perl
use List::Util qw[min max];
use strict;
use File::Basename;
use Getopt::Long;
use Pod::Usage; 
use Data::Dumper;

my $indir = $ENV{'PWD'};
my $outd = $ENV{'PWD'};
my @files = ();
my @outtypes = ("RSeQC");
my @order=( "Total Reads", "Total Tags" , "Total Assigned Tags", "CDS_Exons", "5'UTR_Exons", "3'UTR_Exons", "Introns", "TSS_up_1kb", "TSS_up_5kb", "TSS_up_10kb", "TES_down_1kb", "TES_down_5kb", "TES_down_10kb");
my %lines=(
  "Total Reads" => 1,
  "Total Tags" => 1,
  "Total Assigned Tags" => 1,
  "CDS_Exons" => 2,
  "5'UTR_Exons" => 2,
  "3'UTR_Exons" => 2,
  "Introns" => 2,
  "TSS_up_1kb" => 2,
  "TSS_up_5kb" => 2,
  "TSS_up_10kb" => 2,
  "TES_down_1kb" => 2,
  "TES_down_5kb" => 2,
  "TES_down_10kb" => 2
);


foreach my $outtype (@outtypes)
{

my $ext=".out";
@files = <$indir/$outtype*$ext>;

my @rowheaders=();
my @libs=();
my %vals=();
my %normvals=();
my $type = "rsem";

foreach my $d (@files){
  my $libname=basename($d, $ext);
  $libname=~s/RSeQC.//g;
  $libname=~s/rsem.out.//g;
  $libname=~s/.genome//g;
  print $libname."\\n";
  push(@libs, $libname); 
  getVals($d, $libname, \\%vals, \\%normvals, \\%lines);
}
#print Dumper(%vals);
#print Dumper(%normvals);

my $sizemetrics = keys %vals;
write_results("$outd/$outtype.$type.counts.tsv", \\@libs,\\%vals, \\@order, "region") if ($sizemetrics>0);
write_results("$outd/$outtype.$type.tagskb.tsv", \\@libs,\\%normvals, \\@order, "region") if ($sizemetrics>0);

}

sub write_results
{
  my ($outfile, $libs, $vals, $order, $name )=@_;
  open(OUT, ">$outfile");
  print OUT "$name\\t".join("\\t", @{$libs})."\\n";

  my $lib=${$libs}[0];
  foreach my $key ( @order )
  {
    if (exists ${$vals}{$lib}{$key}) {
    print OUT $key;
    foreach my $lib (@{$libs})
    {
      print OUT "\\t".${$vals}{$lib}{$key};
    } 
    print OUT "\\n";
    }
  }
  close(OUT);
}

sub getVals{
  my ($filename, $libname, $vals, $normvals, $lines)=@_;
  if (-e $filename){
     open(IN, $filename);
     while(my $line=<IN>)
     {
       chomp($line);
       my @vals_arr=split(/\\s{2,}/,$line);
       if (exists ${$lines}{$vals_arr[0]}) {
         my $idx=${$lines}{$vals_arr[0]};
         ${$vals}{$libname}{$vals_arr[0]}=$vals_arr[$idx] if (exists $vals_arr[$idx]);
         if ($idx==2) {
             ${$normvals}{$libname}{$vals_arr[0]}=$vals_arr[3] if (exists $vals_arr[3]);
         }
       }
     } 
  }
  
}
'''

}

//* params.pdfbox_path =  ""  //* @input
//* autofill
if ($HOSTNAME == "default"){
    $CPU  = 1
    $MEMORY = 32
}
//* platform
if ($HOSTNAME == "ghpcc06.umassrc.org"){
    $TIME = 240
    $CPU  = 1
    $MEMORY = 32
    $QUEUE = "short"
}
//* platform
//* autofill

process BAM_Analysis_Module_Kallisto_Picard {

input:
 set val(name), file(bam) from g240_21_bam_file_g224_121

output:
 file "*_metrics"  into g224_121_outputFileOut_g224_82
 file "results/*.pdf"  into g224_121_outputFilePdf_g224_82

when:
(params.run_Picard_CollectMultipleMetrics && (params.run_Picard_CollectMultipleMetrics == "yes")) || !params.run_Picard_CollectMultipleMetrics

script:
"""
picard CollectMultipleMetrics OUTPUT=${name}_multiple.out VALIDATION_STRINGENCY=LENIENT INPUT=${bam}
mkdir results && java -jar ${params.pdfbox_path} PDFMerger *.pdf results/${name}_multi_metrics.pdf
"""
}


process BAM_Analysis_Module_Kallisto_Picard_Summary {

publishDir params.outdir, overwrite: true, mode: 'copy',
	saveAs: {filename ->
	if (filename =~ /.*.tsv$/) "picard_summary_kallisto/$filename"
	else if (filename =~ /results\/.*.pdf$/) "picard_summary_pdf_kallisto/$filename"
}

input:
 file picardOut from g224_121_outputFileOut_g224_82.collect()
 val mate from g_229_mate_g224_82
 file picardPdf from g224_121_outputFilePdf_g224_82.collect()

output:
 file "*.tsv"  into g224_82_outputFileTSV
 file "results/*.pdf"  into g224_82_outputFilePdf

shell:
'''
#!/usr/bin/env perl
use List::Util qw[min max];
use strict;
use File::Basename;
use Getopt::Long;
use Pod::Usage; 
use Data::Dumper;

system("mkdir results && mv *.pdf results/. ");

my $indir = $ENV{'PWD'};
my $outd = $ENV{'PWD'};
my @files = ();
my @outtypes = ("CollectRnaSeqMetrics", "alignment_summary_metrics", "base_distribution_by_cycle_metrics", "insert_size_metrics", "quality_by_cycle_metrics", "quality_distribution_metrics" );

foreach my $outtype (@outtypes)
{
my $ext="_multiple.out";
$ext.=".$outtype" if ($outtype ne "CollectRnaSeqMetrics");
@files = <$indir/*$ext>;

my @rowheaders=();
my @libs=();
my %metricvals=();
my %histvals=();

my $pdffile="";
my $libname="";
foreach my $d (@files){
  my $libname=basename($d, $ext);
  print $libname."\\n";
  push(@libs, $libname); 
  getMetricVals($d, $libname, \\%metricvals, \\%histvals, \\@rowheaders);
}

my $sizemetrics = keys %metricvals;
write_results("$outd/$outtype.stats.tsv", \\@libs,\\%metricvals, \\@rowheaders, "metric") if ($sizemetrics>0);
my $sizehist = keys %histvals;
write_results("$outd/$outtype.hist.tsv", \\@libs,\\%histvals, "none", "nt") if ($sizehist>0);

}

sub write_results
{
  my ($outfile, $libs, $vals, $rowheaders, $name )=@_;
  open(OUT, ">$outfile");
  print OUT "$name\\t".join("\\t", @{$libs})."\\n";
  my $size=0;
  $size=scalar(@{${$vals}{${$libs}[0]}}) if(exists ${$libs}[0] and exists ${$vals}{${$libs}[0]} );
  
  for (my $i=0; $i<$size;$i++)
  { 
    my $rowname=$i;
    $rowname = ${$rowheaders}[$i] if ($name=~/metric/);
    print OUT $rowname;
    foreach my $lib (@{$libs})
    {
      print OUT "\\t".${${$vals}{$lib}}[$i];
    } 
    print OUT "\\n";
  }
  close(OUT);
}

sub getMetricVals{
  my ($filename, $libname, $metricvals, $histvals,$rowheaders)=@_;
  if (-e $filename){
     my $nextisheader=0;
     my $nextisvals=0;
     my $nexthist=0;
     open(IN, $filename);
     while(my $line=<IN>)
     {
       chomp($line);
       @{$rowheaders}=split(/\\t/, $line) if ($nextisheader && !scalar(@{$rowheaders})); 
       if ($nextisvals) {
         @{${$metricvals}{$libname}}=split(/\\t/, $line);
         $nextisvals=0;
       }
       if($nexthist){
          my @vals=split(/[\\s\\t]+/,$line); 
          push(@{${$histvals}{$libname}}, $vals[1]) if (exists $vals[1]);
       }
       $nextisvals=1 if ($nextisheader); $nextisheader=0;
       $nextisheader=1 if ($line=~/METRICS CLASS/);
       $nexthist=1 if ($line=~/normalized_position/);
     } 
  }
  
}
'''

}

//* autofill
//* platform
if ($HOSTNAME == "ghpcc06.umassrc.org"){
    $TIME = 240
    $CPU  = 1
    $MEMORY = 10
    $QUEUE = "short"
}
//* platform
//* autofill

process MultiQC {

publishDir params.outdir, overwrite: true, mode: 'copy',
	saveAs: {filename ->
	if (filename =~ /multiqc_report.html$/) "multiqc/$filename"
}

input:
 file "tophat/*" from g218_11_summary_g_177.flatten().toList()
 file "rsem/*" from g242_22_rsemOut_g_177.flatten().toList()
 file "star/*" from g217_18_logOut_g_177.flatten().toList()
 file "fastqc/*" from g213_3_FastQCout_g_177.flatten().toList()
 file "sequential_mapping/*" from g241_34_bowfiles_g_177.flatten().toList()
 file "rseqc_tophat/*" from g222_122_outputFileOut_g_177.flatten().toList()
 file "rseqc_rsem/*" from g219_122_outputFileOut_g_177.flatten().toList()
 file "rseqc_star/*" from g221_122_outputFileOut_g_177.flatten().toList()
 file "rseqc_hisat/*" from g220_122_outputFileOut_g_177.flatten().toList()
 file "kallisto/*" from g240_21_outputDir_g_177.flatten().toList()
 file "rseqc_kallisto/*" from g224_122_outputFileOut_g_177.flatten().toList()

output:
 file "multiqc_report.html" optional true  into g_177_outputHTML

errorStrategy 'retry'
maxRetries 2

script:
"""
multiqc -e general_stats -d -dd 2 .
"""
}

//* params.gtf =  ""  //* @input

//* autofill
//* platform
if ($HOSTNAME == "ghpcc06.umassrc.org"){
    $TIME = 30
    $CPU  = 1
    $MEMORY = 10
    $QUEUE = "short"
}
//* platform
//* autofill

process Kallisto_module_Kallisto_transcript_to_gene_count {

input:
 file outDir from g240_21_outputDir_g240_28

output:
 file newoutDir  into g240_28_outputDir_g240_27

shell:
newoutDir = "genes_" + outDir
'''
#!/usr/bin/env perl
use strict;
use Getopt::Long;
use IO::File;
use Data::Dumper;

my $gtf_file = "!{params.gtf}";
if (checkFile("!{outDir}/abundance.tsv")){
	rename ("!{outDir}/abundance.tsv", "!{outDir}/abundance_isoforms.tsv");
}
my $kallisto_transcript_matrix_in = "!{outDir}/abundance_isoforms.tsv";
my $kallisto_transcript_matrix_out = "!{outDir}/abundance_genes.tsv";
open(IN, "<$gtf_file") or die "Can't open $gtf_file.\\n";
my %all_genes; # save gene_id of transcript_id
while(<IN>){
  next if(/^##/); #ignore header
  chomp;
  my %attribs = ();
  my ($chr, $source, $type, $start, $end, $score,
    $strand, $phase, $attributes) = split("\\t");
  my @add_attributes = split(";", $attributes);
  # store ids and additional information in second hash
  foreach my $attr ( @add_attributes ) {
     next unless $attr =~ /^\\s*(.+)\\s(.+)$/;
     my $c_type  = $1;
     my $c_value = $2;
     $c_value =~ s/\\"//g;
     if($c_type  && $c_value){
       if(!exists($attribs{$c_type})){
         $attribs{$c_type} = [];
       }
       push(@{ $attribs{$c_type} }, $c_value);
     }
  }
  #work with the information from the two hashes...
  if(exists($attribs{'transcript_id'}->[0]) && exists($attribs{'gene_id'}->[0])){
    if(!exists($all_genes{$attribs{'transcript_id'}->[0]})){
        $all_genes{$attribs{'transcript_id'}->[0]} = $attribs{'gene_id'}->[0];
    }
  } 
}


# print Dumper \\%all_genes;

#Parse the kallisto input file, determine gene IDs for each transcript, and calculate sum TPM values
my %gene_exp;
my %gene_length;
my %samples;
my $ki_fh = IO::File->new($kallisto_transcript_matrix_in, 'r');
my $header = '';
my $h = 0;
while (my $ki_line = $ki_fh->getline) {
  $h++;
  chomp($ki_line);
  my @ki_entry = split("\\t", $ki_line);
  my $s = 0;
  if ($h == 1){
    $header = $ki_line;
    my $first_col = shift @ki_entry;
    my $second_col = shift @ki_entry;
    foreach my $sample (@ki_entry){
      $s++;
      $samples{$s}{name} = $sample;
    }
    next;
  }
  my $trans_id = shift @ki_entry;
  my $length = shift @ki_entry;
  my $gene_id;
  if ($all_genes{$trans_id}){
    $gene_id = $all_genes{$trans_id};
  }elsif($trans_id =~ /ERCC/){
    $gene_id = $trans_id;
  }else{
    print "\\n\\nCould not identify gene id from trans id: $trans_id\\n\\n";
  }

  $s = 0;
  foreach my $value (@ki_entry){
    $s++;
    $gene_exp{$gene_id}{$s} += $value;
  }
  if ($gene_length{$gene_id}){
    $gene_length{$gene_id} = $length if ($length > $gene_length{$gene_id});
  }else{
    $gene_length{$gene_id} = $length;
  }

}
$ki_fh->close;

my $ko_fh = IO::File->new($kallisto_transcript_matrix_out, 'w');
unless ($ko_fh) { die('Failed to open file: '. $kallisto_transcript_matrix_out); }

print $ko_fh "$header\\n";
foreach my $gene_id (sort {$a cmp $b} keys %gene_exp){
  print $ko_fh "$gene_id\\t$gene_length{$gene_id}\\t";
  my @vals;
  foreach my $s (sort {$a <=> $b} keys %samples){
     push(@vals, $gene_exp{$gene_id}{$s});
  }
  my $val_string = join("\\t", @vals);
  print $ko_fh "$val_string\\n";
}


$ko_fh->close;
if (checkFile("!{outDir}")){
	rename ("!{outDir}", "!{newoutDir}");
}

sub checkFile {
    my ($file) = @_;
    print "$file\\n";
    return 1 if ( -e $file );
    return 0;
}

'''
}

//* params.gtf =  ""  //* @input

//* autofill
//* platform
if ($HOSTNAME == "ghpcc06.umassrc.org"){
    $TIME = 30
    $CPU  = 1
    $MEMORY = 10
    $QUEUE = "short"
}
//* platform
//* autofill

process Kallisto_module_Kallisto_Summary {

publishDir params.outdir, overwrite: true, mode: 'copy',
	saveAs: {filename ->
	if (filename =~ /.*.tsv$/) "kallisto_count/$filename"
}

input:
 file kallistoOut from g240_28_outputDir_g240_27.collect()

output:
 file "*.tsv"  into g240_27_outputFile

shell:
'''
#!/usr/bin/env perl
use Data::Dumper;
use strict;

### Parse gtf file
my $gtf_file = "!{params.gtf}";
open(IN, "<$gtf_file") or die "Can't open $gtf_file.\\n";
my %all_genes; # save gene_id of transcript_id
my %all_trans; # map transcript_id of genes
while(<IN>){
  next if(/^##/); #ignore header
  chomp;
  my %attribs = ();
  my ($chr, $source, $type, $start, $end, $score,
    $strand, $phase, $attributes) = split("\\t");
  my @add_attributes = split(";", $attributes);
  # store ids and additional information in second hash
  foreach my $attr ( @add_attributes ) {
     next unless $attr =~ /^\\s*(.+)\\s(.+)$/;
     my $c_type  = $1;
     my $c_value = $2;
     $c_value =~ s/\\"//g;
     if($c_type  && $c_value){
       if(!exists($attribs{$c_type})){
         $attribs{$c_type} = [];
       }
       push(@{ $attribs{$c_type} }, $c_value);
     }
  }
  #work with the information from the two hashes...
  if(exists($attribs{'transcript_id'}->[0]) && exists($attribs{'gene_id'}->[0])){
    if(!exists($all_genes{$attribs{'transcript_id'}->[0]})){
        $all_genes{$attribs{'transcript_id'}->[0]} = $attribs{'gene_id'}->[0];
    }
    if(!exists($all_trans{$attribs{'gene_id'}->[0]})){
        $all_trans{$attribs{'gene_id'}->[0]} = $attribs{'transcript_id'}->[0];
    } else {
    	if (index($all_trans{$attribs{'gene_id'}->[0]}, $attribs{'transcript_id'}->[0]) == -1) {
			$all_trans{$attribs{'gene_id'}->[0]} = $all_trans{$attribs{'gene_id'}->[0]} . "," .$attribs{'transcript_id'}->[0];
		}
    	
    }
  } 
}


print Dumper \\%all_trans;



#### Create summary table

my %tf = (
        expected_count => 3,
        tpm => 4
    );

my $indir = $ENV{'PWD'};
my $outdir = $ENV{'PWD'};

my @gene_iso_ar = ("genes","isoforms");
my @tpm_fpkm_expectedCount_ar = ("expected_count", "tpm");
for(my $l = 0; $l <= $#gene_iso_ar; $l++) {
    my $gene_iso = $gene_iso_ar[$l];
    for(my $ll = 0; $ll <= $#tpm_fpkm_expectedCount_ar; $ll++) {
        my $tpm_fpkm_expectedCount = $tpm_fpkm_expectedCount_ar[$ll];

        opendir D, $indir or die "Could not open $indir\\n";
        my @alndirs = sort { $a cmp $b } grep /^genes_kallisto_/, readdir(D);
        closedir D;
    
        my @a=();
        my %b=();
        my %c=();
        my $i=0;
        foreach my $d (@alndirs){ 
            my $dir = "${indir}/$d";
            print $d."\\n";
            my $libname=$d;
            $libname=~s/genes_kallisto_//;
            $i++;
            $a[$i]=$libname;
            open IN,"${dir}/abundance_${gene_iso}.tsv";
            $_=<IN>;
            while(<IN>)
            {
                my @v=split; 
                # $v[0] -> transcript_id
                # $all_genes{$v[0]} -> $gene_id
                if ($gene_iso eq "isoforms"){
                	$c{$v[0]}=$all_genes{$v[0]};
                } elsif ($gene_iso eq "genes"){
                	$c{$v[0]}=$all_trans{$v[0]};
                } 
                $b{$v[0]}{$i}=$v[$tf{$tpm_fpkm_expectedCount}];
                 
            }
            close IN;
        }
        my $outfile="${indir}/"."$gene_iso"."_expression_"."$tpm_fpkm_expectedCount".".tsv";
        open OUT, ">$outfile";
        if ($gene_iso ne "isoforms") {
            print OUT "gene\\ttranscript";
        } else {
            print OUT "transcript\\tgene";
        }
    
        for(my $j=1;$j<=$i;$j++) {
            print OUT "\\t$a[$j]";
        }
        print OUT "\\n";
    
        foreach my $key (keys %b) {
            print OUT "$key\\t$c{$key}";
            for(my $j=1;$j<=$i;$j++){
                print OUT "\\t$b{$key}{$j}";
            }
            print OUT "\\n";
        }
        close OUT;
    }
}

'''
}

//* autofill
//* platform
if ($HOSTNAME == "ghpcc06.umassrc.org"){
    $TIME = 30
    $CPU  = 1
    $MEMORY = 10
    $QUEUE = "short"
}
//* platform
//* autofill

process Kallisto_module_Kallisto_Alignment_Summary {

input:
 file kallistoDir from g240_21_outputDir_g240_22.collect()

output:
 file "kallisto_alignment_sum.tsv"  into g240_22_outFileTSV_g_198

shell:
'''
#!/usr/bin/env perl
use List::Util qw[min max];
use strict;
use File::Basename;
use Getopt::Long;
use Pod::Usage;
use Data::Dumper;
my $indir = $ENV{'PWD'};

opendir D, $indir or die "Could not open $indir";
my @alndirs = sort { $a cmp $b } grep /^kallisto_/, readdir(D);
closedir D;

my @a=();
my %b=();
my %c=();
my $i=0;
my @headers = ();
my %tsv;
foreach my $d (@alndirs){
    my $dir = "${indir}/$d";
    my $libname=$d;
    $libname=~s/kallisto_//;
    my $multimapped;
    my $aligned;
    my $total;
    
    # eg. [quant] processed 24,788 reads, 19,238 reads pseudoaligned
	chomp($total   = `cat ${dir}/kallisto.log | grep 'pseudoaligned' | sed 's/,//g' | awk '{sum+=\\$3} END {print sum}'`);
	chomp($aligned = `cat ${dir}/kallisto.log | grep 'pseudoaligned' | sed 's/,//g' | awk '{sum+=\\$5} END {print sum}'`);
    $tsv{$libname}=[$libname, $total];
    push(@{$tsv{$libname}}, $aligned);
}

push(@headers, "Sample");
push(@headers, "Total Reads");
push(@headers, "Pseudoaligned Reads (Kallisto)");


my @keys = keys %tsv;
my $summary = "kallisto_alignment_sum.tsv";
my $header_string = join("\\t", @headers);
`echo "$header_string" > $summary`;
foreach my $key (@keys){
    my $values = join("\\t", @{ $tsv{$key} });
        `echo "$values" >> $summary`;
}
'''
}

//* autofill
//* platform
if ($HOSTNAME == "ghpcc06.umassrc.org"){
    $TIME = 240
    $CPU  = 1
    $MEMORY = 8
    $QUEUE = "short"
}
//* platform
//* autofill

process Sequential_Mapping_Module_Sequential_Mapping_Summary {

input:
 set val(name), file(bowfile) from g241_34_bowfiles_g241_26
 val mate from g_229_mate_g241_26
 val filtersList from g241_34_filter_g241_26

output:
 file '*.tsv'  into g241_26_outputFileTSV_g241_13
 val "sequential_mapping_sum"  into g241_26_name_g241_13

errorStrategy 'retry'
maxRetries 2

shell:
'''
#!/usr/bin/env perl
open(my \$fh, '>', "!{name}.tsv");
print $fh "Sample\\tGroup\\tTotal Reads\\tReads After Sequential Mapping\\tUniquely Mapped\\tMultimapped\\tMapped\\n";
my @bowArray = split(' ', "!{bowfile}");
my $group= "\\t";
my @filterArray = (!{filtersList});
foreach my $bowitem(@bowArray) {
    # get mapping id
    my @bowAr = $bowitem.split("_");
    $bowCount = $bowAr[0] + -1;
    # if bowfiles ends with underscore (eg. bow_rRNA), parse rRNA as a group.
    my ($RDS_In, $RDS_After, $RDS_Uniq, $RDS_Multi, $ALGN_T, $a, $b, $aPer, $bPer)=(0, 0, 0, 0, 0, 0, 0, 0, 0);
    if ($bowitem =~ m/bow_([^\\.]+)$/){
        $group = "$1\\t";
        open(IN, $bowitem);
        my $i = 0;
        while(my $line=<IN>){
            chomp($line);
            $line=~s/^ +//;
            my @arr=split(/ /, $line);
            $RDS_In=$arr[0] if ($i=~/^1$/);
            # Reads After Filtering column depends on filtering type
            if ($i == 2){
                if ($filterArray[$bowCount] eq "Yes"){
                    $RDS_After=$arr[0];
                } else {
                    $RDS_After=$RDS_In;
                }
            }
            if ($i == 3){
                $a=$arr[0];
                $aPer=$arr[1];
                $aPer=~ s/([()])//g;
                $RDS_Uniq=$arr[0];
            }
            if ($i == 4){
                $b=$arr[0];
                $bPer=$arr[1];
                $bPer=~ s/([()])//g;
                $RDS_Multi=$arr[0];
            }
            $ALGN_T=($a+$b);
            $i++;
        }
        close(IN);
    } elsif ($bowitem =~ m/star_([^\\.]+)$/){
        $group = "$1\\t";
        open(IN2, $bowitem);
        my $multimapped;
		my $aligned;
		my $inputCount;
		chomp($inputCount = `cat $bowitem | grep 'Number of input reads' | awk '{sum+=\\$6} END {print sum}'`);
		chomp($uniqAligned = `cat $bowitem | grep 'Uniquely mapped reads number' | awk '{sum+=\\$6} END {print sum}'`);
		chomp($multimapped = `cat $bowitem | grep 'Number of reads mapped to multiple loci' | awk '{sum+=\\$9} END {print sum}'`);
		## Here we exclude "Number of reads mapped to too many loci" from multimapped reads since in bam file it called as unmapped.
		## Besides, these "too many loci" reads exported as unmapped reads from STAR.
		$RDS_In = int($inputCount);
		$RDS_Multi = int($multimapped);
        $RDS_Uniq = int($uniqAligned);
        $ALGN_T = $RDS_Uniq+$RDS_Multi;
		if ($filterArray[$bowCount] eq "Yes"){
            $RDS_After=$RDS_In-$ALGN_T;
        } else {
            $RDS_After=$RDS_In;
        }
    } elsif ($bowitem =~ m/bow1_([^\\.]+)$/){
        $group = "$1\\t";
        open(IN2, $bowitem);
        my $multimapped;
		my $aligned;
		my $inputCount;
		my $uniqAligned;
		chomp($inputCount = `cat $bowitem | grep '# reads processed:' | awk '{sum+=\\$4} END {print sum}'`);
		chomp($aligned = `cat $bowitem | grep '# reads with at least one reported alignment:' | awk '{sum+=\\$9} END {print sum}'`);
		chomp($uniqAligned = `cat $bowitem | grep '# unique mapped reads:' | awk '{sum+=\\$5} END {print sum}'`);
		## Here we exclude "Number of reads mapped to too many loci" from multimapped reads since in bam file it called as unmapped.
		## Besides, these "too many loci" reads exported as unmapped reads from STAR.
		$RDS_In = int($inputCount);
		$RDS_Multi = int($aligned) -int($uniqAligned);
		if ($RDS_Multi < 0 ){
		    $RDS_Multi = 0;
		}
        $RDS_Uniq = int($uniqAligned);
        $ALGN_T = int($aligned);
		if ($filterArray[$bowCount] eq "Yes"){
            $RDS_After=$RDS_In-$ALGN_T;
        } else {
            $RDS_After=$RDS_In;
        }
    }
    
    print $fh "!{name}\\t$group$RDS_In\\t$RDS_After\\t$RDS_Uniq\\t$RDS_Multi\\t$ALGN_T\\n";
}
close($fh);



'''

}


process Sequential_Mapping_Module_Merge_TSV_Files {

input:
 file tsv from g241_26_outputFileTSV_g241_13.collect()
 val outputFileName from g241_26_name_g241_13.collect()

output:
 file "${name}.tsv"  into g241_13_outputFileTSV_g241_14

errorStrategy 'retry'
maxRetries 3

script:
name = outputFileName[0]
"""    
awk 'FNR==1 && NR!=1 {  getline; } 1 {print} ' *.tsv > ${name}.tsv
"""
}


process Sequential_Mapping_Module_Sequential_Mapping_Short_Summary {

publishDir params.outdir, overwrite: true, mode: 'copy',
	saveAs: {filename ->
	if (filename =~ /sequential_mapping_short_sum.tsv$/) "sequential_mapping_summary/$filename"
	else if (filename =~ /sequential_mapping_detailed_sum.tsv$/) "sequential_mapping_summary/$filename"
}

input:
 file mainSum from g241_13_outputFileTSV_g241_14

output:
 file "sequential_mapping_short_sum.tsv"  into g241_14_outputFileTSV_g_198
 file "sequential_mapping_detailed_sum.tsv"  into g241_14_outputFile

errorStrategy 'retry'
maxRetries 2

shell:
'''
#!/usr/bin/env perl
use List::Util qw[min max];
use File::Basename;
use Getopt::Long;
use Pod::Usage;
use Data::Dumper;

my @header;
my %all_rows;
my @seen_cols_short;
my @seen_cols_detailed;
my $ID_header;

chomp(my $contents = `ls *.tsv`);
my @files = split(/[\\n]+/, $contents);
foreach my $file (@files){
        open IN,"$file";
        my $line1 = <IN>;
        chomp($line1);
        ( $ID_header, my @h) = ( split("\\t", $line1) );
        my $totalHeader = $h[1];
        my $afterFilteringHeader = $h[2];
        my $uniqueHeader = $h[3];
        my $multiHeader = $h[4];
        my $mappedHeader = $h[5];
        push(@seen_cols_short, $totalHeader) unless grep{$_ eq $totalHeader} @seen_cols_short; #Total reads Header
        push(@seen_cols_detailed, $totalHeader) unless grep{$_ eq $totalHeader} @seen_cols_detailed; #Total reads Header

        my $n=0;
        while (my $line=<IN>) {
                
                chomp($line);
                my ( $ID, @fields ) = ( split("\\t", $line) ); 
                #SHORT
                push(@seen_cols_short, $fields[0]) unless grep{$_ eq $fields[0]} @seen_cols_short; #mapped element header
                $all_rows{$ID}{$fields[0]} = $fields[5];#Mapped Reads
                #Grep first line $fields[1] as total reads.
                if (!exists $all_rows{$ID}{$totalHeader}){    
                        $all_rows{$ID}{$totalHeader} = $fields[1];
                } 
                $all_rows{$ID}{$afterFilteringHeader} = $fields[2]; #only use last entry
                #DETAILED
                $uniqueHeadEach = "$fields[0] (${uniqueHeader})";
                $multiHeadEach = "$fields[0] (${multiHeader})";
                $mappedHeadEach = "$fields[0] (${mappedHeader})";
                push(@seen_cols_detailed, $mappedHeadEach) unless grep{$_ eq $mappedHeadEach} @seen_cols_detailed;
                push(@seen_cols_detailed, $uniqueHeadEach) unless grep{$_ eq $uniqueHeadEach} @seen_cols_detailed;
                push(@seen_cols_detailed, $multiHeadEach) unless grep{$_ eq $multiHeadEach} @seen_cols_detailed;
                $all_rows{$ID}{$mappedHeadEach} = $fields[5];
                $all_rows{$ID}{$uniqueHeadEach} = $fields[3];
                $all_rows{$ID}{$multiHeadEach} = $fields[4];
    }
    close IN;
    push(@seen_cols_short, $afterFilteringHeader) unless grep{$_ eq $afterFilteringHeader} @seen_cols_short; #After filtering Header
}


#print Dumper \\%all_rows;
#print Dumper \\%seen_cols_short;

printFiles("sequential_mapping_short_sum.tsv",@seen_cols_short,);
printFiles("sequential_mapping_detailed_sum.tsv",@seen_cols_detailed);


sub printFiles {
    my($summary, @cols_to_print) = @_;
    
    open OUT, ">$summary";
    print OUT join ("\\t", $ID_header,@cols_to_print),"\\n";
    foreach my $key ( keys %all_rows ) { 
        print OUT join ("\\t", $key, (map { $all_rows{$key}{$_} // '' } @cols_to_print)),"\\n";
        }
        close OUT;
}

'''


}

g217_11_outputFileTSV_g_198= g217_11_outputFileTSV_g_198.ifEmpty([""]) 
g241_14_outputFileTSV_g_198= g241_14_outputFileTSV_g_198.ifEmpty([""]) 
g216_10_outputFileTSV_g_198= g216_10_outputFileTSV_g_198.ifEmpty([""]) 
g242_17_outputFileTSV_g_198= g242_17_outputFileTSV_g_198.ifEmpty([""]) 
g218_9_outputFileTSV_g_198= g218_9_outputFileTSV_g_198.ifEmpty([""]) 
g213_11_outputFileTSV_g_198= g213_11_outputFileTSV_g_198.ifEmpty([""]) 
g213_21_outputFileTSV_g_198= g213_21_outputFileTSV_g_198.ifEmpty([""]) 
g213_16_outputFileTSV_g_198= g213_16_outputFileTSV_g_198.ifEmpty([""]) 
g240_22_outFileTSV_g_198= g240_22_outFileTSV_g_198.ifEmpty([""]) 

//* autofill
//* platform
if ($HOSTNAME == "ghpcc06.umassrc.org"){
    $TIME = 30
    $CPU  = 1
    $MEMORY = 10
    $QUEUE = "short"
}
//* platform
//* autofill

process Overall_Summary {

publishDir params.outdir, overwrite: true, mode: 'copy',
	saveAs: {filename ->
	if (filename =~ /overall_summary.tsv$/) "summary/$filename"
}

input:
 file starSum from g217_11_outputFileTSV_g_198
 file sequentialSum from g241_14_outputFileTSV_g_198
 file hisatSum from g216_10_outputFileTSV_g_198
 file rsemSum from g242_17_outputFileTSV_g_198
 file tophatSum from g218_9_outputFileTSV_g_198
 file adapterSum from g213_11_outputFileTSV_g_198
 file trimmerSum from g213_21_outputFileTSV_g_198
 file qualitySum from g213_16_outputFileTSV_g_198
 file kallistoSum from g240_22_outFileTSV_g_198

output:
 file "overall_summary.tsv"  into g_198_outputFileTSV

shell:
'''
#!/usr/bin/env perl
use List::Util qw[min max];
use strict;
use File::Basename;
use Getopt::Long;
use Pod::Usage;
use Data::Dumper;

my @header;
my %all_rows;
my @seen_cols;
my $ID_header;

chomp(my $contents = `ls *.tsv`);
my @rawFiles = split(/[\\n]+/, $contents);
my @files = ();
# order must be in this order for chipseq pipeline: bowtie->dedup
# rsem bam pipeline: dedup->rsem, star->dedup
# riboseq ncRNA_removal->star
my @order = ("adapter_removal","trimmer","quality","extractUMI","extractValid","sequential_mapping","ncRNA_removal","bowtie","star","hisat2","tophat2", "dedup","rsem","kallisto","esat","count");
for ( my $k = 0 ; $k <= $#order ; $k++ ) {
    for ( my $i = 0 ; $i <= $#rawFiles ; $i++ ) {
        if ( $rawFiles[$i] =~ /$order[$k]/ ) {
            push @files, $rawFiles[$i];
        }
    }
}

print Dumper \\@files;
##add rest of the files
for ( my $i = 0 ; $i <= $#rawFiles ; $i++ ) {
    push(@files, $rawFiles[$i]) unless grep{$_ == $rawFiles[$i]} @files;
}
print Dumper \\@files;

##Merge each file according to array order

foreach my $file (@files){
        open IN,"$file";
        my $line1 = <IN>;
        chomp($line1);
        ( $ID_header, my @header) = ( split("\\t", $line1) );
        push @seen_cols, @header;

        while (my $line=<IN>) {
        chomp($line);
        my ( $ID, @fields ) = ( split("\\t", $line) ); 
        my %this_row;
        @this_row{@header} = @fields;

        #print Dumper \\%this_row;

        foreach my $column (@header) {
            if (! exists $all_rows{$ID}{$column}) {
                $all_rows{$ID}{$column} = $this_row{$column}; 
            }
        }   
    }
    close IN;
}

#print for debugging
#print Dumper \\%all_rows;
#print Dumper \\%seen_cols;

#grab list of column headings we've seen, and order them. 
my @cols_to_print = uniq(@seen_cols);
my $summary = "overall_summary.tsv";
open OUT, ">$summary";
print OUT join ("\\t", $ID_header,@cols_to_print),"\\n";
foreach my $key ( keys %all_rows ) { 
    #map iterates all the columns, and gives the value or an empty string. if it's undefined. (prevents errors)
    print OUT join ("\\t", $key, (map { $all_rows{$key}{$_} // '' } @cols_to_print)),"\\n";
}
close OUT;

sub uniq {
    my %seen;
    grep ! $seen{$_}++, @_;
}

'''


}

//* autofill
//* platform
if ($HOSTNAME == "ghpcc06.umassrc.org"){
    $TIME = 30
    $CPU  = 1
    $MEMORY = 10
    $QUEUE = "short"
}
//* platform
//* autofill

process RSEM_module_RSEM_Count {

publishDir params.outdir, overwrite: true, mode: 'copy',
	saveAs: {filename ->
	if (filename =~ /.*.tsv$/) "rsem_summary/$filename"
}

input:
 file rsemOut from g242_22_rsemOut_g242_21.collect()

output:
 file "*.tsv"  into g242_21_outputFile
 val wdir  into g242_21_wdir_g242_14

shell:
wdir="rsem_summary/genes_expression_expected_count.tsv"
if (params.gtf){
	gtf=params.gtf
}
if (params.senseantisense){
	senseantisense = params.senseantisense
}
'''
#!/usr/bin/env perl

my %tf = (
        expected_count => 4,
        tpm => 5,
        fpkm => 6,
    );

my $indir = $ENV{'PWD'};
$outdir = $ENV{'PWD'};

my @sense_antisense = ("", ".forward", ".reverse");
my @gene_iso_ar = ("genes", "isoforms");
my @tpm_fpkm_expectedCount_ar = ("expected_count", "tpm");
for($m = 0; $m <= $#sense_antisense; $m++) {
	my $sa = $sense_antisense[$m];
	for($l = 0; $l <= $#gene_iso_ar; $l++) {
	    my $gene_iso = $gene_iso_ar[$l];
	    for($ll = 0; $ll <= $#tpm_fpkm_expectedCount_ar; $ll++) {
	        my $tpm_fpkm_expectedCount = $tpm_fpkm_expectedCount_ar[$ll];
	
	        opendir D, $indir or die "Could not open $indir\n";
	        my @alndirs = sort { $a cmp $b } grep /^pipe/, readdir(D);
	        closedir D;
	    
	        my @a=();
	        my %b=();
	        my %c=();
	        my $i=0;
	        my $saexist=0;
	        foreach my $d (@alndirs){ 
	            my $dir = "${indir}/$d";
	            print $d."\n";
	            my $libname=$d;
	            $libname=~s/pipe\\.rsem\\.//;
	    
	            $i++;
	            $a[$i]=$libname;
	            if (-e "${dir}/rsem.out$sa.$libname.$gene_iso.results"){
	            	$saexist=1;
	         
		            open IN,"${dir}/rsem.out$sa.$libname.$gene_iso.results";
		            $_=<IN>;
		            while(<IN>)
		            {
		                my @v=split; 
		                $b{$v[0]}{$i}=$v[$tf{$tpm_fpkm_expectedCount}];
		                $c{$v[0]}=$v[1];
		            }
		            close IN;
	            }
	        }
	        if ($saexist==1){
		        my $outfile="${indir}/${gene_iso}${sa}_expression_"."$tpm_fpkm_expectedCount".".tsv";
		        open OUT, ">$outfile";
		        if ($gene_iso ne "isoforms") {
		            print OUT "gene\ttranscript";
		        } else {
		            print OUT "transcript\tgene";
		        }
		    
		        for(my $j=1;$j<=$i;$j++) {
		            print OUT "\t$a[$j]";
		        }
		        print OUT "\n";
		    
		        foreach my $key (keys %b) {
		            print OUT "$key\t$c{$key}";
		            for(my $j=1;$j<=$i;$j++){
		                print OUT "\t$b{$key}{$j}";
		            }
		            print OUT "\n";
		        }
		        close OUT;
		        if ($sa eq ".reverse"){
		        	my $com = "perl !{senseantisense} !{gtf} ${indir}/${gene_iso}.forward_expression_${tpm_fpkm_expectedCount}.tsv ${indir}/${gene_iso}${sa}_expression_${tpm_fpkm_expectedCount}.tsv";
		        	`$com`;
		        	$com = "rm -rf ${indir}/${gene_iso}.forward_expression_${tpm_fpkm_expectedCount}.tsv ${indir}/${gene_iso}${sa}_expression_${tpm_fpkm_expectedCount}.tsv";
		        	`$com`;
		        }
	        }
	    }
	}
}
'''
}

cols = params.RSEM_module_CountData_DE.cols
conds = params.RSEM_module_CountData_DE.conds
cols = convertCommaSepString(cols)
conds = convertCommaSepString(conds)

// convert comma separated string into comma quote and comma separated string
//eg. "a, b, c" -> "a","b","c"
def convertCommaSepString( t ) {
    commaSepList=t.split(",").collect{ '"' + it.trim() + '"'}
    c=commaSepList.toString()
    //remove first and last brackets
    return c.substring(1, c.length()-1)
}

//* autofill
//* platform
if ($HOSTNAME == "ghpcc06.umassrc.org"){
    $TIME = 10
    $CPU  = 1
    $MEMORY = 10
    $QUEUE = "short"
}
//* platform
//* autofill

process RSEM_module_CountData_DE {

publishDir params.outdir, overwrite: true, mode: 'copy',
	saveAs: {filename ->
	if (filename =~ /.*.rmd$/) "rsem_rmarkdown/$filename"
}

input:
 val wdir from g242_21_wdir_g242_14

output:
 file "*.rmd"  into g242_14_rMarkdown

shell:
'''
#!/usr/bin/env perl

my $script = <<'EOF';
## Count data DE analysis

[1]. Reading the data.

The merged the count data table will be read from a web URL to be able 
to run this rmarkdown anywhere. 
In Eukaryotes only a subset of all genes are expressed in 
a given cell. Expression is therefore a bimodal distribution, 
with non-expressed genes having counts that result from experimental 
and biological noise. It is important to filter out the genes 
that are not expressed before doing differential gene expression. 
You can decide which cutoff separates expressed vs non-expressed 
genes by looking your histogram we created.


```{r, echo=FALSE, message=FALSE}
library(debrowser)
library(plotly)
source("https://dolphinnext.umassmed.edu/dist/scripts/funcs.R")
library(RCurl)
url<-{{webpath:"!{wdir}"}}
file <- textConnection(getURL(url)) 
rsem <- read.table(file,sep="\\t", header=TRUE, row.names=1) 
data <- data.frame(rsem[,sapply(rsem, is.numeric)]) 

cols<- c(!{cols})

data <- data[, cols]

h <- hist(log10(rowSums(data)), breaks = as.numeric(100), plot = FALSE) 

plot_ly(x = h$mids, y = h$counts, width = 500, height=300) %>% 
layout( title = "Histogram") %>%
add_bars()
``` 

[2]. All2all scatter plots

To check the reproducibility of biological replicates, we use all2all plots.

```{r, echo=FALSE, message=FALSE}
all2all(data)
``` 

[3]. DESeq ANALYSIS

The goal of Differential gene expression analysis is to find 
genes or transcripts whose difference in expression, when accounting 
for the variance within condition, is higher than expected by chance. 

The first step is to indicate the condition that each column (experiment) 
in the table represent. 
Here we define the correspondence between columns and conditions. 
Make sure the order of the columns matches to your table.

In this case a total sum of 10 counts separates well expressed 
from non-expressed genes. You can change this value and padj value and 
log2FoldChange cutoffs according to your data

```{r, echo=FALSE, message=FALSE}
conds <- factor( c(!{conds}) )
avgall<-cbind(rowSums(data[cols[conds == levels(conds)[1]]])/3, 
              rowSums(data[cols[conds == levels(conds)[2]]])/3)
colnames(avgall)<-c(levels(conds)[1], levels(conds)[2])

gdat<-data.frame(avgall)
de_res <- runDESeq(data, cols, conds,  padj=0.01, log2FoldChange=1, non_expressed_cutoff=10)
overlaid_data <- overlaySig(gdat, de_res$res_selected)
ggplot() +
  geom_point(data=overlaid_data, aes_string(x=levels(conds)[1], y=levels(conds)[2],
                                            colour="Legend"), alpha=6/10, size=3) +
  scale_colour_manual(values=c("All"="darkgrey","Significant"="red"))+
  scale_x_log10() +scale_y_log10()
```

[4]. MA Plot

The Second way to visualize it, we use MA plots.
For MA Plot there is another builtin function that you can use

```{r, echo=FALSE, message=FALSE}
plotMA(de_res$res_detected,ylim=c(-2,2),main="DESeq2");
```

[5]. Volcano Plot

The third way of visualizing the data is making a Volcano Plot.
Here on the x axis you have log2foldChange values and y axis you 
have your -log10 padj values. To see how significant genes are 
distributed. Highlight genes that have an absolute fold change > 2 
and a padj < 0.01

```{r, echo=FALSE, message=FALSE}
volcanoPlot(de_res,  padj=0.01, log2FoldChange=1)
```

[6] Heatmap

The forth way of visualizing the data that is widely used in this 
type of analysis is clustering and Heatmaps.

```{r, echo=FALSE, message=FALSE}
sel_data<-data[rownames(de_res$res_selected),]
norm_data<-getNormalizedMatrix(sel_data, method="TMM")
ld <- log2(norm_data+0.1)
cldt <- scale(t(ld), center=TRUE, scale=TRUE);
cld <- t(cldt)
dissimilarity <- 1 - cor(cld)
distance <- as.dist(dissimilarity)
heatmap.2(cld, Rowv=TRUE,dendrogram="column",
          Colv=TRUE, col=redblue(256),labRow=NA,
          density.info="none",trace="none", cexCol=0.8,
          hclust=function(x) hclust(x,method="complete"),
          distfun=function(x) as.dist((1-cor(t(x)))/2))

```
EOF

open OUT, ">rmark.rmd";
print OUT $script;
close OUT;
'''


}


workflow.onComplete {
println "##Pipeline execution summary##"
println "---------------------------"
println "##Completed at: $workflow.complete"
println "##Duration: ${workflow.duration}"
println "##Success: ${workflow.success ? 'OK' : 'failed' }"
println "##Exit status: ${workflow.exitStatus}"
}
