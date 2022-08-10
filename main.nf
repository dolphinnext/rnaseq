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

def pathChecker(input, path, type){
	recursive = (type == "folder") ? "--recursive" : ""
	cmd = "mkdir -p check && mv ${input} check/. "
	if (!input || input.empty()){
		input = file(path).getName().toString()
		cmd = "mkdir -p check && cd check && ln -s ${path} ${input} && cd .."
		if (path.indexOf('s3:') > -1 || path.indexOf('S3:') >-1){
			cmd = "mkdir -p check && cd check && aws s3 cp ${recursive} ${path} ${workDir}/${input} && ln -s ${workDir}/${input} . && cd .."
		} else if (path.indexOf('/') == -1){
			cmd = ""
		}
}
	return [cmd,input]
}
if (!params.run_FeatureCounts_after_STAR){params.run_FeatureCounts_after_STAR = ""} 
if (!params.run_FeatureCounts_after_Hisat2){params.run_FeatureCounts_after_Hisat2 = ""} 
if (!params.run_FeatureCounts_after_Tophat2){params.run_FeatureCounts_after_Tophat2 = ""} 
if (!params.run_FeatureCounts_after_RSEM){params.run_FeatureCounts_after_RSEM = ""} 
if (!params.run_FeatureCounts_after_Kallisto){params.run_FeatureCounts_after_Kallisto = ""} 
if (!params.mate){params.mate = ""} 
if (!params.reads){params.reads = ""} 
// Stage empty file to be used as an optional input where required
ch_empty_file_1 = file("$baseDir/.emptyfiles/NO_FILE_1", hidden:true)
ch_empty_file_2 = file("$baseDir/.emptyfiles/NO_FILE_2", hidden:true)
ch_empty_file_3 = file("$baseDir/.emptyfiles/NO_FILE_3", hidden:true)
ch_empty_file_4 = file("$baseDir/.emptyfiles/NO_FILE_4", hidden:true)
ch_empty_file_5 = file("$baseDir/.emptyfiles/NO_FILE_5", hidden:true)
ch_empty_file_6 = file("$baseDir/.emptyfiles/NO_FILE_6", hidden:true)
ch_empty_file_7 = file("$baseDir/.emptyfiles/NO_FILE_7", hidden:true)
ch_empty_file_8 = file("$baseDir/.emptyfiles/NO_FILE_8", hidden:true)
ch_empty_file_9 = file("$baseDir/.emptyfiles/NO_FILE_9", hidden:true)
ch_empty_file_10 = file("$baseDir/.emptyfiles/NO_FILE_10", hidden:true)
ch_empty_file_11 = file("$baseDir/.emptyfiles/NO_FILE_11", hidden:true)

Channel.value(params.run_FeatureCounts_after_STAR).set{g_179_run_featureCounts0_g253_125}
Channel.value(params.run_FeatureCounts_after_Hisat2).set{g_188_run_featureCounts0_g252_125}
Channel.value(params.run_FeatureCounts_after_Tophat2).set{g_189_run_featureCounts0_g254_125}
Channel.value(params.run_FeatureCounts_after_RSEM).set{g_203_run_featureCounts0_g251_125}
Channel.value(params.run_FeatureCounts_after_Kallisto).set{g_225_run_featureCounts0_g255_125}
Channel.value(params.mate).into{g_229_mate0_g_127;g_229_mate1_g247_3;g_229_mate0_g247_14;g_229_mate0_g248_36;g_229_mate0_g249_14;g_229_mate0_g250_26;g_229_mate1_g251_82;g_229_mate1_g251_95;g_229_mate0_g251_131;g_229_mate1_g251_133;g_229_mate1_g252_82;g_229_mate1_g252_95;g_229_mate0_g252_131;g_229_mate1_g252_133;g_229_mate1_g253_82;g_229_mate1_g253_95;g_229_mate0_g253_131;g_229_mate1_g253_133;g_229_mate1_g254_82;g_229_mate1_g254_95;g_229_mate0_g254_131;g_229_mate1_g254_133;g_229_mate1_g255_82;g_229_mate1_g255_95;g_229_mate0_g255_131;g_229_mate1_g255_133;g_229_mate1_g256_26;g_229_mate1_g256_30;g_229_mate1_g256_46;g_229_mate1_g257_11;g_229_mate1_g257_16;g_229_mate1_g257_18;g_229_mate1_g257_19;g_229_mate1_g257_20;g_229_mate1_g257_21;g_229_mate1_g257_24;g_229_mate1_g257_23;g_229_mate0_g257_28;g_229_mate0_g257_31;g_229_mate0_g264_31;g_229_mate2_g264_30}
if (params.reads){
Channel
	.fromFilePairs( params.reads , size: params.mate == "single" ? 1 : params.mate == "pair" ? 2 : params.mate == "triple" ? 3 : params.mate == "quadruple" ? 4 : -1 )
	.ifEmpty { error "Cannot find any reads matching: ${params.reads}" }
	.set{g_230_reads1_g257_28}
 } else {  
	g_230_reads1_g257_28 = Channel.empty()
 }


//* params.gtf =  ""  //* @input
//* params.genome =  ""  //* @input
//* params.commondb =  ""  //* @input
//* params.genome_source =  ""  //* @input
//* params.gtf_source =  ""  //* @input
//* params.commondb_source =  ""  //* @input

def downFile(path, task){
    if (path.take(1).indexOf("/") == 0){
      target=path
      if (task.executor == "awsbatch") {
      	a=file(path)
    	fname = a.getName().toString()
    	target = "${workDir}/${fname}"
    	if (!file(target).exists()){
    		a.copyTo(workDir)
    	}
      }
    } else {
      a=file(path)
      fname = a.getName().toString()
      target = "${workDir}/${fname}"
      if (!file(target).exists()){
    		a.copyTo(workDir)
      } 
    }
    return target
}

def getLastName (str){
	if (str.indexOf("/") > -1){
		return  str.substring(str.lastIndexOf('/')+1,str.length())
	} 
	return ""
}

process Check_and_Build_Module_Check_Genome_GTF {


output:
 file "${genomeName}"  into g245_21_genome00_g245_52, g245_21_genome01_g245_54
 file "${gtfName}"  into g245_21_gtfFile10_g245_53, g245_21_gtfFile10_g245_54

when:
params.run_Download_Genomic_Sources == "yes"

script:
gtf_dir = ""
genome_dir = ""
if (params.gtf.indexOf('/') > -1){
	gtf_dir  = params.gtf.substring(0, params.gtf.lastIndexOf('/')) 
}
if (params.genome.indexOf('/') > -1){
	genome_dir  = params.genome.substring(0, params.genome.lastIndexOf('/')) 
}

downGenomePath = ""
downGtfPath = ""
downPathPrefix = (task.executor == "awsbatch") ? 's3:/' : ''
genomeSource = params.genome
gtfSource = params.gtf
genomeSource = !file("${params.genome}").exists() ? params.genome_source : params.genome
downGenomePath = (!file("${params.genome}").exists() || task.executor == "awsbatch" ) ? downPathPrefix + downFile(genomeSource, task) : ""
genomeName = getLastName(genomeSource)
gtfSource = !file("${params.gtf}").exists() ? params.gtf_source : params.gtf
downGtfPath= (!file("${params.gtf}").exists() || task.executor == "awsbatch" ) ? downPathPrefix + downFile(gtfSource, task) : ""
gtfName = getLastName(gtfSource)

"""
if [ ! -e "${params.genome}" ] ; then
    echo "${params.genome} not found"
    
    if [ "${task.executor}" != "awsbatch" ] ; then
    	if [  "${genome_dir}" != "" ] ; then
    		mkdir -p ${genome_dir}
    		cp -n $downGenomePath ${params.genome}
    		ln -s ${params.genome} ${genomeName}
    	else 
    		ln -s $downGenomePath ${genomeName}
    	fi
    else
    	aws s3 cp $downGenomePath ${genomeName}
    fi
else 
	ln -s ${params.genome} ${genomeName}
fi

if [ ! -e "${params.gtf}" ] ; then
    echo "${params.gtf} not found"
    
    if [ "${task.executor}" != "awsbatch" ] ; then
    	if [  "${gtf_dir}" != "" ] ; then
    		mkdir -p ${gtf_dir}
    		cp -n $downGtfPath ${params.gtf}
    		ln -s ${params.gtf} ${gtfName}
    	else 
    		ln -s $downGtfPath ${gtfName}
    	fi
    else
    	aws s3 cp $downGtfPath ${gtfName}
    fi
else 
	ln -s ${params.gtf} ${gtfName}
fi

"""




}

//* params.gtf2bed_path =  ""  //* @input

process Check_and_Build_Module_Check_BED12 {

input:
 file gtf from g245_21_gtfFile10_g245_53

output:
 file "${gtfName}.bed"  into g245_53_bed03_g245_54

when:
params.run_Download_Genomic_Sources == "yes"

script:
gtfName  = gtf.baseName
beddir = ""
if (params.bed.indexOf('/') > -1){
	beddir  = params.bed.substring(0, params.bed.lastIndexOf('/')) 
}
"""

if [ ! -e "${params.bed}" ] ; then
    echo "${params.bed} not found"
    perl ${params.gtf2bed_path} $gtf > ${gtfName}.bed
else 
	cp -n ${params.bed} ${gtfName}.bed
fi
if [ "${beddir}" != "" ] ; then
	mkdir -p ${beddir}
	cp -n ${gtfName}.bed ${params.bed} 
fi
"""




}

//* params.gtf2bed_path =  ""  //* @input
//* params.genome_sizes =  ""  //* @input

process Check_and_Build_Module_Check_chrom_sizes_and_index {

input:
 file genome from g245_21_genome00_g245_52

output:
 file "${genomeName}.chrom.sizes"  into g245_52_genomeSizes02_g245_54

when:
params.run_Download_Genomic_Sources == "yes"

script:
genomeName  = genome.baseName
genome_sizes_dir = ""
if (params.genome_sizes.indexOf('/') > -1){
	genome_sizes_dir  = params.genome_sizes.substring(0, params.genome_sizes.lastIndexOf('/')) 
}

"""
if [ ! -e "${params.genome_sizes}" ] ; then
    echo "${params.genome_sizes} not found"
    cat ${genome} | awk '\$0 ~ ">" {print c; c=0;printf substr(\$1,2,100) "\\t"; } \$0 !~ ">" {c+=length(\$0);} END { print c; }' > ${genomeName}.chrom.sizes
    ##clean first empty line
    sed -i '1{/^\$/d}' ${genomeName}.chrom.sizes
    if [ "${genome_sizes_dir}" != "" ] ; then
    	mkdir -p ${genome_sizes_dir}
		cp -n ${genomeName}.chrom.sizes ${params.genome_sizes} 
	fi
else 
	cp ${params.genome_sizes} ${genomeName}.chrom.sizes
fi

"""




}

g245_21_gtfFile10_g245_54= g245_21_gtfFile10_g245_54.ifEmpty([""]) 
g245_21_genome01_g245_54= g245_21_genome01_g245_54.ifEmpty([""]) 
g245_52_genomeSizes02_g245_54= g245_52_genomeSizes02_g245_54.ifEmpty([""]) 
g245_53_bed03_g245_54= g245_53_bed03_g245_54.ifEmpty([""]) 


process Check_and_Build_Module_check_files {

input:
 file gtf from g245_21_gtfFile10_g245_54
 file genome from g245_21_genome01_g245_54
 file genomeSizes from g245_52_genomeSizes02_g245_54
 file bed from g245_53_bed03_g245_54

output:
 file "*/${gtf}" optional true  into g245_54_gtfFile01_g247_16, g245_54_gtfFile03_g247_14, g245_54_gtfFile01_g250_31, g245_54_gtfFile01_g249_16, g245_54_gtfFile03_g248_36, g245_54_gtfFile01_g248_38, g245_54_gtfFile01_g248_39, g245_54_gtfFile01_g248_32, g245_54_gtfFile03_g251_133, g245_54_gtfFile03_g252_133, g245_54_gtfFile03_g253_133, g245_54_gtfFile03_g254_133, g245_54_gtfFile03_g255_133, g245_54_gtfFile01_g256_47, g245_54_gtfFile01_g264_21
 file "*/${genome}" optional true  into g245_54_genome10_g247_16, g245_54_genome10_g250_31, g245_54_genome10_g249_16, g245_54_genome10_g248_32, g245_54_genome10_g256_47, g245_54_genome10_g264_21
 file "*/${genomeSizes}" optional true  into g245_54_genomeSizes24_g248_36, g245_54_genomeSizes22_g251_131, g245_54_genomeSizes21_g251_142, g245_54_genomeSizes22_g252_131, g245_54_genomeSizes21_g252_142, g245_54_genomeSizes22_g253_131, g245_54_genomeSizes21_g253_142, g245_54_genomeSizes22_g255_131, g245_54_genomeSizes21_g255_142, g245_54_genomeSizes22_g254_131, g245_54_genomeSizes21_g254_142
 file "*/${bed}" optional true  into g245_54_bed31_g251_134, g245_54_bed31_g252_134, g245_54_bed31_g253_134, g245_54_bed31_g254_134, g245_54_bed31_g255_134

script:
(cmd1, gtf) = pathChecker(gtf, params.gtf, "file")
(cmd2, genome) = pathChecker(genome, params.genome, "file")
(cmd3, genomeSizes) = pathChecker(genomeSizes, params.genome_sizes, "file")
(cmd4, bed) = pathChecker(bed, params.bed, "file")
"""
$cmd1
$cmd2
$cmd3
$cmd4
"""
}

build_STAR_index = params.STAR_Module_Check_Build_STAR_Index.build_STAR_index

//* autofill
if ($HOSTNAME == "default"){
    $CPU  = 1
    $MEMORY = 50
}
//* platform
if ($HOSTNAME == "ghpcc06.umassrc.org"){
    $TIME = 3000
    $CPU  = 1
    $MEMORY = 50
    $QUEUE = "long"
}
//* platform
//* autofill

process STAR_Module_Check_Build_STAR_Index {

publishDir params.outdir, mode: 'copy', saveAs: {filename -> if (filename =~ /STARIndex$/) "star_index/$filename"}
input:
 file genome from g245_54_genome10_g264_21
 file gtf from g245_54_gtfFile01_g264_21

output:
 file "STARIndex"  into g264_21_starIndex00_g264_26

when:
build_STAR_index == true && ((params.run_STAR && (params.run_STAR == "yes")) || !params.run_STAR)

script:
star_build_parameters = params.STAR_Module_Check_Build_STAR_Index.star_build_parameters
index_dir = ""
if (params.star_index.indexOf('/') > -1 && params.star_index.indexOf('s3://') < 0){
	index_dir  = file(params.star_index).getParent()
}
newDirName = "STARIndex" 
"""
if [ ! -e "${params.star_index}/SA" ] ; then
    echo "STAR index not found"
    mkdir -p $newDirName 
    STAR --runMode genomeGenerate ${star_build_parameters} --genomeDir $newDirName --genomeFastaFiles ${genome} --sjdbGTFfile ${gtf}
	if [ "${index_dir}" != "" ] ; then
		mkdir -p ${index_dir}
		cp -R $newDirName  ${params.star_index}
	fi
else 
	ln -s ${params.star_index} STARIndex
fi

"""





}

g264_21_starIndex00_g264_26= g264_21_starIndex00_g264_26.ifEmpty([""]) 


if (!((params.run_STAR && (params.run_STAR == "yes")) || !params.run_STAR)){
g264_21_starIndex00_g264_26.set{g264_26_starIndex02_g264_31}
} else {

process STAR_Module_check_STAR_files {

input:
 file star from g264_21_starIndex00_g264_26

output:
 file "*/${star}" optional true  into g264_26_starIndex02_g264_31

when:
(params.run_STAR && (params.run_STAR == "yes")) || !params.run_STAR

script:
(cmd, star) = pathChecker(star, params.star_index, "folder")
"""
$cmd
"""
}
}


download_build_sequential_mapping_indexes = params.Sequential_Mapping_Module_Download_build_sequential_mapping_indexes.download_build_sequential_mapping_indexes

//* autofill
if ($HOSTNAME == "default"){
    $CPU  = 1
    $MEMORY = 50
}
//* platform
if ($HOSTNAME == "ghpcc06.umassrc.org"){
    $TIME = 3000
    $CPU  = 1
    $MEMORY = 50
    $QUEUE = "long"
}
//* platform
//* autofill

process Sequential_Mapping_Module_Download_build_sequential_mapping_indexes {

publishDir params.outdir, mode: 'copy', saveAs: {filename -> if (filename =~ /${commondbName}$/) "commondb/$filename"}
input:
 file genome from g245_54_genome10_g256_47
 file gtf from g245_54_gtfFile01_g256_47

output:
 file "${commondbName}"  into g256_47_commondb00_g256_43
 file "${bowtieIndex}" optional true  into g256_47_bowtieIndex11_g256_43
 file "${bowtie2Index}" optional true  into g256_47_bowtie2index22_g256_43
 file "${starIndex}" optional true  into g256_47_starIndex33_g256_43

when:
download_build_sequential_mapping_indexes == true && params.run_Sequential_Mapping == "yes"

script:
slashCount = params.commondb_source.count("/")
cutDir = slashCount - 3;
commondbSource = !file("${params.commondb}").exists() ? params.commondb_source : params.commondb
commondbName = "commondb"
inputName = file("${params.commondb}").getName().toString()

selectedSeqList = params.Sequential_Mapping_Module_Sequential_Mapping._select_sequence
alignerList = params.Sequential_Mapping_Module_Sequential_Mapping._aligner
genomeIndexes = selectedSeqList.findIndexValues { it ==~ /genome/ }
buildBowtieIndex = alignerList[genomeIndexes].contains("bowtie")
buildBowtie2Index = alignerList[genomeIndexes].contains("bowtie2")
buildSTARIndex = alignerList[genomeIndexes].contains("STAR")

basename = genome.baseName
bowtie2Index = "Bowtie2Index" 
bowtieIndex = "BowtieIndex" 
starIndex = "STARIndex"
bowtie_index_dir = ""
bowtie2_index_dir = ""
star_index_dir = ""
if (params.bowtie_index.indexOf('/') > -1 && params.bowtie_index.indexOf('s3://') < 0){
	bowtie_index_dir  = params.bowtie_index.substring(0, params.bowtie_index.lastIndexOf('/')) 
}
if (params.bowtie2_index.indexOf('/') > -1 && params.bowtie2_index.indexOf('s3://') < 0){
	bowtie2_index_dir  = params.bowtie2_index.substring(0, params.bowtie2_index.lastIndexOf('/')) 
}
if (params.star_index.indexOf('/') > -1 && params.star_index.indexOf('s3://') < 0){
	star_index_dir  = params.star_index.substring(0, params.star_index.lastIndexOf('/')) 
}

"""
if [ ! -e "${params.commondb}" ] ; then
    echo "${params.commondb} not found"
	if [[ "${params.commondb_source}" =~ "s3" ]]; then
		echo "Downloading s3 path from ${params.commondb_source}"
		aws s3 cp --recursive ${params.commondb_source} ${workDir}/${commondbName} && ln -s ${workDir}/${commondbName} ${commondbName}
	else
		echo "Downloading commondb with wget"
		wget --no-check-certificate --secure-protocol=TLSv1 -l inf -nc -nH --cut-dirs=$cutDir -R 'index.html*' -r --no-parent --directory-prefix=\$PWD/${commondbName} ${params.commondb_source}
	fi

else 
	ln -s ${params.commondb} ${commondbName}
fi


if [ "${buildBowtie2Index}" == "true" ]; then
	if [ ! -e "${bowtie2_index_dir}/${basename}.rev.1.bt2" ] ; then
    	echo "${bowtie2_index_dir}/${basename}.rev.1.bt2 Bowtie2 index not found"
    	mkdir -p $bowtie2Index && cp $genome $gtf $bowtie2Index/. && cd $bowtie2Index
    	bowtie2-build ${genome} ${basename}
    	cd ..
    	if [ "${bowtie2_index_dir}" != "" ] ; then
			mkdir -p ${bowtie2_index_dir}
			cp -R -n $bowtie2Index  ${bowtie2_index_dir}
		fi
	else 
		ln -s ${bowtie2_index_dir} $bowtie2Index
	fi
fi

if [ "${buildBowtieIndex}" == "true" ]; then
	if [ ! -e "${bowtie_index_dir}/${basename}.rev.2.ebwt" ] ; then
    	echo "${bowtie_index_dir}/${basename}.rev.2.ebwt Bowtie index not found"
    	mkdir -p $bowtieIndex && cp $genome $gtf $bowtieIndex/. && cd $bowtieIndex
    	bowtie-build ${genome} ${basename}
    	cd ..
    	if [ "${bowtie_index_dir}" != "" ] ; then
			mkdir -p ${bowtie_index_dir}
			cp -R -n $bowtieIndex  ${bowtie_index_dir}
		fi
	else 
		ln -s ${bowtie_index_dir} $bowtieIndex
	fi
fi 

if [ "${buildSTARIndex}" == "true" ]; then
	if [ ! -e "${params.star_index}/SA" ] ; then
    	echo "STAR index not found"
    	mkdir -p $starIndex 
    	STAR --runMode genomeGenerate --genomeDir $starIndex --genomeFastaFiles ${genome} --sjdbGTFfile ${gtf}
		if [ "${star_index_dir}" != "" ] ; then
			mkdir -p ${star_index_dir}
			cp -R $starIndex  ${params.star_index}
		fi
	else 
		ln -s ${params.star_index} $starIndex
	fi
fi
"""




}

//* @style @array:{run_name,run_parameters} @multicolumn:{run_name,run_parameters}

process BAM_Analysis_Module_Kallisto_featureCounts_Prep {

input:
 val run_featureCounts from g_225_run_featureCounts0_g255_125

output:
 val run_params  into g255_125_run_parameters02_g255_133

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

//* @style @array:{run_name,run_parameters} @multicolumn:{run_name,run_parameters}

process BAM_Analysis_Module_Tophat2_featureCounts_Prep {

input:
 val run_featureCounts from g_189_run_featureCounts0_g254_125

output:
 val run_params  into g254_125_run_parameters02_g254_133

when:
run_featureCounts == "yes"

script:
run_name = params.BAM_Analysis_Module_Tophat2_featureCounts_Prep.run_name
run_parameters = params.BAM_Analysis_Module_Tophat2_featureCounts_Prep.run_parameters
sense_antisense = params.BAM_Analysis_Module_Tophat2_featureCounts_Prep.sense_antisense

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

process BAM_Analysis_Module_STAR_featureCounts_Prep {

input:
 val run_featureCounts from g_179_run_featureCounts0_g253_125

output:
 val run_params  into g253_125_run_parameters02_g253_133

when:
run_featureCounts == "yes"

script:
run_name = params.BAM_Analysis_Module_STAR_featureCounts_Prep.run_name
run_parameters = params.BAM_Analysis_Module_STAR_featureCounts_Prep.run_parameters
sense_antisense = params.BAM_Analysis_Module_STAR_featureCounts_Prep.sense_antisense

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

process BAM_Analysis_Module_HISAT2_featureCounts_Prep {

input:
 val run_featureCounts from g_188_run_featureCounts0_g252_125

output:
 val run_params  into g252_125_run_parameters02_g252_133

when:
run_featureCounts == "yes"

script:
run_name = params.BAM_Analysis_Module_HISAT2_featureCounts_Prep.run_name
run_parameters = params.BAM_Analysis_Module_HISAT2_featureCounts_Prep.run_parameters
sense_antisense = params.BAM_Analysis_Module_HISAT2_featureCounts_Prep.sense_antisense

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

process BAM_Analysis_Module_RSEM_featureCounts_Prep {

input:
 val run_featureCounts from g_203_run_featureCounts0_g251_125

output:
 val run_params  into g251_125_run_parameters02_g251_133

when:
run_featureCounts == "yes"

script:
run_name = params.BAM_Analysis_Module_RSEM_featureCounts_Prep.run_name
run_parameters = params.BAM_Analysis_Module_RSEM_featureCounts_Prep.run_parameters
sense_antisense = params.BAM_Analysis_Module_RSEM_featureCounts_Prep.sense_antisense

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

build_Kallisto_index = params.Kallisto_module_Check_Build_Kallisto_Index.build_Kallisto_index

process Kallisto_module_Check_Build_Kallisto_Index {

publishDir params.outdir, mode: 'copy', saveAs: {filename -> if (filename =~ /$index\/.*.idx$/) "kallisto_index/$filename"}
input:
 file genome from g245_54_genome10_g248_32
 file gtf from g245_54_gtfFile01_g248_32

output:
 file "$index/*.idx"  into g248_32_kallisto_index00_g248_31

when:
build_Kallisto_index == true && ((params.run_Kallisto && (params.run_Kallisto == "yes")) || !params.run_Kallisto)

script:
index_dir = ""
if (params.kallisto_index.indexOf('/') > -1 && params.kallisto_index.indexOf('s3://') < 0){
	index_dir = file(params.kallisto_index).getParent()
}
index = "KallistoIndex" 

"""
if [ ! -e "${index_dir}/transcripts.idx" ] ; then
    echo "${index_dir}/transcripts.idx Kallisto index not found"
    
    mkdir -p $index && mv $genome $gtf $index/. && cd $index
    filter_gtf_for_genes_in_genome.py --gtf ${gtf} --fasta ${genome} -o genome_filtered_genes.gtf
    gffread -F -w transcripts_raw.fa -g ${genome} genome_filtered_genes.gtf
    cut -d ' ' -f1 transcripts_raw.fa > transcripts.fa
    gzip transcripts.fa
    kallisto index -i transcripts.idx transcripts.fa.gz
    
    if [ "${index_dir}" != "" ] ; then
    	cd ..
		mkdir -p ${index_dir}
		cp -R -n $index  ${index_dir}
	fi
else 
	ln -s ${index_dir} $index
fi
"""



}

g248_32_kallisto_index00_g248_31= g248_32_kallisto_index00_g248_31.ifEmpty([""]) 


if (!((params.run_Kallisto && (params.run_Kallisto == "yes")) || !params.run_Kallisto)){
g248_32_kallisto_index00_g248_31.set{g248_31_kallisto_index02_g248_36}
} else {

process Kallisto_module_check_kallisto_files {

input:
 file kallisto from g248_32_kallisto_index00_g248_31

output:
 file "*/${kallisto}" optional true  into g248_31_kallisto_index02_g248_36

when:
(params.run_Kallisto && (params.run_Kallisto == "yes")) || !params.run_Kallisto

script:
(cmd, kallisto) = pathChecker(kallisto, params.kallisto_index, "file")
"""
$cmd
"""
}
}


build_Hisat2_index = params.HISAT2_Module_Check_Build_Hisat2_Index.build_Hisat2_index
//* autofill
if ($HOSTNAME == "default"){
    $CPU  = 1
    $MEMORY = 50
}
//* platform
if ($HOSTNAME == "ghpcc06.umassrc.org"){
    $TIME = 3000
    $CPU  = 1
    $MEMORY = 50
    $QUEUE = "long"
}
//* platform
//* autofill

process HISAT2_Module_Check_Build_Hisat2_Index {

publishDir params.outdir, mode: 'copy', saveAs: {filename -> if (filename =~ /$index$/) "hisat2_index/$filename"}
input:
 file genome from g245_54_genome10_g249_16
 file gtf from g245_54_gtfFile01_g249_16

output:
 file "$index"  into g249_16_hisat2Index00_g249_15

when:
build_Hisat2_index == true && ((params.run_HISAT2 && (params.run_HISAT2 == "yes")) || !params.run_HISAT2)

script:
hisat2_build_parameters = params.HISAT2_Module_Check_Build_Hisat2_Index.hisat2_build_parameters
basename = genome.baseName
basenameGTF = gtf.baseName
index_dir = ""
if (params.hisat2_index.indexOf('/') > -1 && params.hisat2_index.indexOf('s3://') < 0){
	index_dir  = file(params.hisat2_index).getParent()
}
index = "Hisat2Index" 

extract_splice_sites = "hisat2_extract_splice_sites.py ${gtf} > ${basenameGTF}.hisat2_splice_sites.txt"
extract_exons = "hisat2_extract_exons.py ${gtf}> ${basenameGTF}.hisat2_exons.txt"
ss = "--ss ${basenameGTF}.hisat2_splice_sites.txt"
exon = "--exon ${basenameGTF}.hisat2_exons.txt"

"""
if [ ! -e "${index_dir}/${basename}.8.ht2" ] ; then
    echo "${index_dir}/${basename}.8.ht2 Hisat2 index not found"
    
    mkdir -p $index && mv $genome $gtf $index/. && cd $index
    $extract_splice_sites
    $extract_exons
    hisat2-build ${hisat2_build_parameters} $ss $exon ${genome} ${basename}
    cd ..
    if [ "${index_dir}" != "" ] ; then
		mkdir -p ${index_dir}
		cp -R -n $index  ${index_dir}
	fi
else 
	ln -s ${index_dir} $index
fi
"""




}

g249_16_hisat2Index00_g249_15= g249_16_hisat2Index00_g249_15.ifEmpty([""]) 


if (!((params.run_HISAT2 && (params.run_HISAT2 == "yes")) || !params.run_HISAT2)){
g249_16_hisat2Index00_g249_15.set{g249_15_hisat2Index02_g249_14}
} else {

process HISAT2_Module_check_Hisat2_files {

input:
 file hisat2 from g249_16_hisat2Index00_g249_15

output:
 file "*/${hisat2}" optional true  into g249_15_hisat2Index02_g249_14

when:
(params.run_HISAT2 && (params.run_HISAT2 == "yes")) || !params.run_HISAT2

script:
(cmd, hisat2) = pathChecker(hisat2, params.hisat2_index, "folder")
"""
$cmd
"""
}
}


build_RSEM_index = params.RSEM_module_Check_Build_Rsem_Index.build_RSEM_index

//* autofill
if ($HOSTNAME == "default"){
    $CPU  = 1
    $MEMORY = 50
}
//* platform
if ($HOSTNAME == "ghpcc06.umassrc.org"){
    $TIME = 3000
    $CPU  = 1
    $MEMORY = 50
    $QUEUE = "long"
}
//* platform
//* autofill

process RSEM_module_Check_Build_Rsem_Index {

publishDir params.outdir, mode: 'copy', saveAs: {filename -> if (filename =~ /${index}$/) "rsem_index/$filename"}
input:
 file genome from g245_54_genome10_g250_31
 file gtf from g245_54_gtfFile01_g250_31

output:
 file "${index}"  into g250_31_rsemIndex00_g250_32

when:
build_RSEM_index == true && ((params.run_RSEM && (params.run_RSEM == "yes")) || !params.run_RSEM)

script:
transcript_to_gene_map = params.RSEM_module_Check_Build_Rsem_Index.transcript_to_gene_map
RSEM_build_parameters = params.RSEM_module_Check_Build_Rsem_Index.RSEM_build_parameters

transcript_to_gene_mapText = ""
if (transcript_to_gene_map?.trim()){
    transcript_to_gene_mapText = "--transcript-to-gene-map " + transcript_to_gene_map
}
RSEM_reference_type = params.RSEM_module_RSEM.RSEM_reference_type
basename = genome.baseName

indexType = ""
index = ""
index_dir = ""
if (RSEM_reference_type == 'bowtie'){
	indexType = "--bowtie "
	index = "RSEM_ref_Bowtie" 
	if (params.rsem_ref_using_bowtie_index.indexOf('/') > -1 && params.rsem_ref_using_bowtie_index.indexOf('s3://') < 0){
		index_dir = params.rsem_ref_using_bowtie_index
	}
} else if (RSEM_reference_type == 'bowtie2'){
	indexType = "--bowtie2 "
	index = "RSEM_ref_Bowtie2" 
	if (params.rsem_ref_using_bowtie2_index.indexOf('/') > -1 && params.rsem_ref_using_bowtie2_index.indexOf('s3://') < 0){
		index_dir = params.rsem_ref_using_bowtie2_index
	}
} else if (RSEM_reference_type == 'star'){
	indexType = "--star "
	index = "RSEM_ref_STAR" 
	if (params.rsem_ref_using_star_index.indexOf('/') > -1 && params.rsem_ref_using_star_index.indexOf('s3://') < 0){
		index_dir = params.rsem_ref_using_star_index
	}
}

"""
if [ ! -e "${index_dir}/${basename}.ti" ] ; then
    echo "${index_dir}/${basename}.ti RSEM index not found"
    
    mkdir -p $index && mv $genome $gtf $index/. && cd $index
    rsem-prepare-reference ${RSEM_build_parameters} --gtf ${gtf} ${transcript_to_gene_mapText} ${indexType} ${genome} ${basename}
    cd ..
    if [ "${index_dir}" != "" ] ; then
		mkdir -p ${index_dir}
		cp -R -n $index  ${index_dir}
	fi
else 
	ln -s ${index_dir} $index
fi
"""

}

build_Bowtie2_index = params.Tophat2_Module_Check_Build_Bowtie2_Index.build_Bowtie2_index
bowtie2_build_parameters = params.Tophat2_Module_Check_Build_Bowtie2_Index.bowtie2_build_parameters

//* autofill
if ($HOSTNAME == "default"){
    $CPU  = 1
    $MEMORY = 30
}
//* platform
if ($HOSTNAME == "ghpcc06.umassrc.org"){
    $TIME = 3000
    $CPU  = 1
    $MEMORY = 30
    $QUEUE = "long"
}
//* platform
//* autofill

process Tophat2_Module_Check_Build_Bowtie2_Index {

publishDir params.outdir, mode: 'copy', saveAs: {filename -> if (filename =~ /$index$/) "tophat2_index/$filename"}
input:
 file genome from g245_54_genome10_g247_16
 file gtf from g245_54_gtfFile01_g247_16

output:
 file "$index"  into g247_16_bowtie2index00_g247_18

when:
build_Bowtie2_index == true && ((params.run_Tophat && (params.run_Tophat == "yes")) || !params.run_Tophat)

script:
basename = genome.baseName
index_dir = ""
if (params.bowtie2_index.indexOf('/') > -1 && params.bowtie2_index.indexOf('s3://') < 0){
	index_dir  = file(params.bowtie2_index).getParent()
}
index = "Bowtie2Index" 

"""
if [ ! -e "${index_dir}/${basename}.rev.1.bt2" ] ; then
    echo "${index_dir}/${basename}.rev.1.bt2 Bowtie2 index not found"
    
    mkdir -p $index && mv $genome $gtf $index/. && cd $index
    bowtie2-build ${bowtie2_build_parameters} ${genome} ${basename}
    if [ "${index_dir}" != "" ] ; then
    	cd ..
		mkdir -p ${index_dir}
		cp -R -n $index  ${index_dir}
	fi
else 
	ln -s ${index_dir} $index
fi
"""



}

g247_16_bowtie2index00_g247_18= g247_16_bowtie2index00_g247_18.ifEmpty([""]) 


if (!((params.run_Tophat && (params.run_Tophat == "yes")) || !params.run_Tophat)){
g247_16_bowtie2index00_g247_18.set{g247_18_bowtie2index02_g247_14}
} else {

process Tophat2_Module_check_Tophat2_files {

input:
 file bowtie2 from g247_16_bowtie2index00_g247_18

output:
 file "*/${bowtie2}" optional true  into g247_18_bowtie2index02_g247_14

when:
(params.run_Tophat && (params.run_Tophat == "yes")) || !params.run_Tophat

script:
(cmd, bowtie2) = pathChecker(bowtie2, params.bowtie2_index, "folder")
"""
$cmd
"""
}
}


g250_31_rsemIndex00_g250_32= g250_31_rsemIndex00_g250_32.ifEmpty([""]) 


if (!((params.run_RSEM && (params.run_RSEM == "yes")) || !params.run_RSEM)){
g250_31_rsemIndex00_g250_32.set{g250_32_rsemIndex02_g250_26}
} else {

process RSEM_module_check_RSEM_files {

input:
 file rsem from g250_31_rsemIndex00_g250_32

output:
 file "*/${rsem}" optional true  into g250_32_rsemIndex02_g250_26

when:
(params.run_RSEM && (params.run_RSEM == "yes")) || !params.run_RSEM

script:
RSEM_reference_type = ''
if (params.RSEM && params.RSEM.RSEM_reference_type){
	RSEM_reference_type = params.RSEM.RSEM_reference_type		
} else if (params.RSEM_module_RSEM && params.RSEM_module_RSEM.RSEM_reference_type){
	RSEM_reference_type = params.RSEM_module_RSEM.RSEM_reference_type
}

if (RSEM_reference_type == 'bowtie'){
	systemInput = params.rsem_ref_using_bowtie_index
} else if (RSEM_reference_type == 'bowtie2'){
	systemInput = params.rsem_ref_using_bowtie2_index
} else if (RSEM_reference_type == 'star'){
	systemInput = params.rsem_ref_using_star_index
}	

	
(cmd, rsem) = pathChecker(rsem, systemInput, "folder")
"""
$cmd
"""
}
}


g256_47_commondb00_g256_43= g256_47_commondb00_g256_43.ifEmpty([""]) 
g256_47_bowtieIndex11_g256_43= g256_47_bowtieIndex11_g256_43.ifEmpty([""]) 
g256_47_bowtie2index22_g256_43= g256_47_bowtie2index22_g256_43.ifEmpty([""]) 
g256_47_starIndex33_g256_43= g256_47_starIndex33_g256_43.ifEmpty([""]) 

//* params.gtf =  ""  //* @input
//* params.genome =  ""  //* @input
//* params.commondb =  ""  //* @input
if (!(params.run_Sequential_Mapping  == "yes")){
g256_47_commondb00_g256_43.into{g256_43_commondb02_g256_46; g256_43_commondb05_g256_44; g256_43_commondb05_g256_45}
g256_47_bowtieIndex11_g256_43.into{g256_43_bowtieIndex13_g256_46; g256_43_bowtieIndex12_g256_44; g256_43_bowtieIndex12_g256_45}
g256_47_bowtie2index22_g256_43.into{g256_43_bowtie2index24_g256_46; g256_43_bowtie2index23_g256_44; g256_43_bowtie2index23_g256_45}
g256_47_starIndex33_g256_43.into{g256_43_starIndex35_g256_46; g256_43_starIndex34_g256_44; g256_43_starIndex34_g256_45}
} else {


process Sequential_Mapping_Module_Check_Sequential_Mapping_Indexes {

input:
 file commondb from g256_47_commondb00_g256_43
 file bowtieIndex from g256_47_bowtieIndex11_g256_43
 file bowtie2Index from g256_47_bowtie2index22_g256_43
 file starIndex from g256_47_starIndex33_g256_43

output:
 file "*/${commondb}" optional true  into g256_43_commondb02_g256_46, g256_43_commondb05_g256_44, g256_43_commondb05_g256_45
 file "*/${bowtieIndex}" optional true  into g256_43_bowtieIndex13_g256_46, g256_43_bowtieIndex12_g256_44, g256_43_bowtieIndex12_g256_45
 file "*/${bowtie2Index}" optional true  into g256_43_bowtie2index24_g256_46, g256_43_bowtie2index23_g256_44, g256_43_bowtie2index23_g256_45
 file "*/${starIndex}" optional true  into g256_43_starIndex35_g256_46, g256_43_starIndex34_g256_44, g256_43_starIndex34_g256_45

when:
params.run_Sequential_Mapping  == "yes"

script:
(cmd1, commondb) = pathChecker(commondb, params.commondb, "folder")
(cmd2, bowtieIndex) = pathChecker(bowtieIndex, params.bowtie_index, "folder")
(cmd3, bowtie2Index) = pathChecker(bowtie2Index, params.bowtie2_index, "folder")
(cmd4, starIndex) = pathChecker(starIndex, params.star_index, "folder")
"""
$cmd1
$cmd2
$cmd3
$cmd4
"""
}
}




if (!((params.run_FastQC && (params.run_FastQC == "yes")))){
g_230_reads1_g257_28.set{g257_28_reads10_g257_18}
g257_28_FastQCout04_g_177 = Channel.empty()
} else {

process Adapter_Trimmer_Quality_Module_FastQC {

publishDir params.outdir, mode: 'copy', saveAs: {filename -> if (filename =~ /.*.(html|zip)$/) "fastqc/$filename"}
input:
 val mate from g_229_mate0_g257_28
 set val(name), file(reads) from g_230_reads1_g257_28

output:
 file '*.{html,zip}'  into g257_28_FastQCout04_g_177
 set val(name), file("reads/*")  into g257_28_reads10_g257_18

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
mkdir reads
mv ${file} reads/.
"""
}
}


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
g257_28_reads10_g257_18.into{g257_18_reads00_g257_23; g257_18_reads01_g257_31}
g257_18_log_file10_g257_11 = Channel.empty()
} else {


process Adapter_Trimmer_Quality_Module_Adapter_Removal {

input:
 set val(name), file(reads) from g257_28_reads10_g257_18
 val mate from g_229_mate1_g257_18

output:
 set val(name), file("reads/*.fastq")  into g257_18_reads00_g257_23, g257_18_reads01_g257_31
 file "*.{fastx,trimmomatic}.log"  into g257_18_log_file10_g257_11

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

runCmd("!{runGzip}");
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
            runCmd("rm -f $targetFile");
            print "INFO: $targetFile deleted.\\n";
        }
    }
}


##Subroutines
sub runCmd {
    my ($com) = @_;
    if ($com eq ""){
		return "";
    }
    my $error = system(@_);
    if   ($error) { die "Command failed: $error $com\\n"; }
    else          { print "Command successful: $com\\n"; }
}
'''

}
}



if (!((params.run_FastQC && params.run_FastQC == "yes" && params.run_Adapter_Removal && params.run_Adapter_Removal == "yes"))){
g257_18_reads01_g257_31.set{g257_31_reads11}
g257_31_FastQCout015_g_177 = Channel.empty()
} else {

process Adapter_Trimmer_Quality_Module_FastQC_after_Adapter_Removal {

input:
 val mate from g_229_mate0_g257_31
 set val(name), file(reads) from g257_18_reads01_g257_31

output:
 file '*.{html,zip}'  into g257_31_FastQCout015_g_177
 set val(name), file("reads/*")  into g257_31_reads11

errorStrategy 'retry'
maxRetries 3

when:
(params.run_FastQC && params.run_FastQC == "yes" && params.run_Adapter_Removal && params.run_Adapter_Removal == "yes")

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
mkdir reads
mv ${file} reads/.
"""
}
}



process Adapter_Trimmer_Quality_Module_Adapter_Removal_Summary {

publishDir params.outdir, mode: 'copy', saveAs: {filename -> if (filename =~ /adapter_removal_detailed_summary.tsv$/) "adapter_removal_detailed_summary/$filename"}
input:
 file logfile from g257_18_log_file10_g257_11.collect()
 val mate from g_229_mate1_g257_11

output:
 file "adapter_removal_summary.tsv"  into g257_11_outputFileTSV05_g_198
 file "adapter_removal_detailed_summary.tsv" optional true  into g257_11_outputFile11

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

//* @style @condition:{single_or_paired_end_reads="single", barcode_pattern1,remove_duplicates_based_on_UMI}, {single_or_paired_end_reads="pair", barcode_pattern1,barcode_pattern2}

if (!(params.run_UMIextract == "yes" || !params.run_UMIextract)){
g257_18_reads00_g257_23.set{g257_23_reads00_g257_19}
g257_23_log_file10_g257_24 = Channel.empty()
} else {


process Adapter_Trimmer_Quality_Module_UMIextract {

input:
 set val(name), file(reads) from g257_18_reads00_g257_23
 val mate from g_229_mate1_g257_23

output:
 set val(name), file("result/*.fastq")  into g257_23_reads00_g257_19
 file "${name}.*.log"  into g257_23_log_file10_g257_24

when:
params.run_UMIextract == "yes" || !params.run_UMIextract

script:
readArray = reads.toString().split(' ')
file2 = ""
file1 =  readArray[0]
if (mate == "pair") {file2 =  readArray[1]}


single_or_paired_end_reads = params.Adapter_Trimmer_Quality_Module_UMIextract.single_or_paired_end_reads
barcode_pattern1 = params.Adapter_Trimmer_Quality_Module_UMIextract.barcode_pattern1
barcode_pattern2 = params.Adapter_Trimmer_Quality_Module_UMIextract.barcode_pattern2
UMIqualityFilterThreshold = params.Adapter_Trimmer_Quality_Module_UMIextract.UMIqualityFilterThreshold
phred = params.Adapter_Trimmer_Quality_Module_UMIextract.phred
remove_duplicates_based_on_UMI = params.Adapter_Trimmer_Quality_Module_UMIextract.remove_duplicates_based_on_UMI

"""
set +e
source activate umi_tools_env 2> /dev/null || true
mkdir result
if [ "${mate}" == "pair" ]; then
umi_tools extract --bc-pattern='${barcode_pattern1}' \
                  --bc-pattern2='${barcode_pattern2}' \
                  --extract-method=regex \
                  --stdin=${file1} \
                  --stdout=result/${name}_R1.fastq \
                  --read2-in=${file2} \
                  --read2-out=result/${name}_R2.fastq\
				  --quality-filter-threshold=${UMIqualityFilterThreshold} \
				  --quality-encoding=phred${phred} \
				  --log=${name}.umitools.log 


else
umi_tools extract --bc-pattern='${barcode_pattern1}' \
                  --log=${name}.umitools.log \
                  --extract-method=regex \
                  --stdin ${file1} \
                  --stdout result/${name}.fastq \
				  --quality-filter-threshold=${UMIqualityFilterThreshold} \
				  --quality-encoding=phred${phred}
	if [ "${remove_duplicates_based_on_UMI}" == "true" ]; then		  
        mv result/${name}.fastq  result/${name}_umitools.fastq 
        ## only checks last part of the underscore splitted header for UMI
        awk '(NR%4==1){name=\$1;header=\$0;len=split(name,umiAr,"_");umi=umiAr[len];} (NR%4==2){total++;if(a[umi]!=1){nondup++;a[umi]=1;  print header;print;getline; print; getline; print;}} END{print FILENAME"\\t"total"\\t"nondup > "${name}.dedup.log"}' result/${name}_umitools.fastq > result/${name}.fastq
        rm result/${name}_umitools.fastq
	fi			  
fi
"""

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
g257_23_reads00_g257_19.set{g257_19_reads00_g257_20}
g257_19_log_file10_g257_21 = Channel.empty()
} else {


process Adapter_Trimmer_Quality_Module_Trimmer {

input:
 set val(name), file(reads) from g257_23_reads00_g257_19
 val mate from g_229_mate1_g257_19

output:
 set val(name), file("reads/*q")  into g257_19_reads00_g257_20
 file "*.log" optional true  into g257_19_log_file10_g257_21

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
 
runCmd("mkdir reads");
runCmd("!{runGzip}");
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
            runCmd("rm -f $targetFile");
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
          runCmd("mv $file reads/.");
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
    if ($com eq ""){
		return "";
    }
    my $error = system(@_);
    if   ($error) { die "Command failed: $error $com\\n"; }
    else          { print "Command successful: $com\\n"; }
}

'''

}
}



process Adapter_Trimmer_Quality_Module_Trimmer_Removal_Summary {

input:
 file logfile from g257_19_log_file10_g257_21.collect()
 val mate from g_229_mate1_g257_21

output:
 file "trimmer_summary.tsv"  into g257_21_outputFileTSV06_g_198

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
g257_19_reads00_g257_20.set{g257_20_reads00_g256_46}
g257_20_log_file10_g257_16 = Channel.empty()
} else {


process Adapter_Trimmer_Quality_Module_Quality_Filtering {

input:
 set val(name), file(reads) from g257_19_reads00_g257_20
 val mate from g_229_mate1_g257_20

output:
 set val(name), file("reads/*q")  into g257_20_reads00_g256_46
 file "*.{fastx,trimmomatic}_quality.log" optional true  into g257_20_log_file10_g257_16

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
 
runCmd("mkdir reads unpaired");
runCmd("!{runGzip}");
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
        runCmd("mv !{file1} !{file2} reads/.");
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
            runCmd("rm -f $targetFile");
            print "INFO: $targetFile deleted.\\n";
        }
    }
}

##Subroutines
sub runCmd {
    my ($com) = @_;
    if ($com eq ""){
		return "";
    }
    my $error = system(@_);
    if   ($error) { die "Command failed: $error $com\\n"; }
    else          { print "Command successful: $com\\n"; }
}


'''

}
}


g256_43_commondb02_g256_46= g256_43_commondb02_g256_46.ifEmpty([""]) 
g256_43_bowtieIndex13_g256_46= g256_43_bowtieIndex13_g256_46.ifEmpty([""]) 
g256_43_bowtie2index24_g256_46= g256_43_bowtie2index24_g256_46.ifEmpty([""]) 
g256_43_starIndex35_g256_46= g256_43_starIndex35_g256_46.ifEmpty([""]) 

//* params.bowtie_index =  ""  //* @input
//* params.bowtie2_index =  ""  //* @input
//* params.star_index =  ""  //* @input

//both bowtie and bowtie2 indexes located in same path
bowtieIndexes = [rRNA:  "commondb/rRNA/rRNA", ercc:  "commondb/ercc/ercc", miRNA: "commondb/miRNA/miRNA", tRNA:  "commondb/tRNA/tRNA", piRNA: "commondb/piRNA/piRNA", snRNA: "commondb/snRNA/snRNA", rmsk:  "commondb/rmsk/rmsk"]
genomeIndexes = [bowtie: "BowtieIndex", bowtie2: "Bowtie2Index", STAR: "STARIndex"]


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

selectSequenceList = []
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
selectSequenceList = (_select_sequence)

mappingList = mapList.join(" ") // convert into space separated format in order to use in bash for loop
paramsList = paramList.join(",") // convert into comma separated format in order to use in as array in bash
alignersList = alignerList.join(",") 
filtersList = filterList.join(",") 
indexesList = indexList.join(",") 
senseList = senseList.join(",")
selectSequencesList = selectSequenceList.join(",")

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
g257_20_reads00_g256_46.into{g256_46_reads01_g_127; g256_46_reads01_g248_36; g256_46_reads01_g250_26}
g256_46_bowfiles10_g256_26 = Channel.empty()
g256_46_bowfiles15_g_177 = Channel.empty()
g256_46_bam_file20_g256_44 = Channel.empty()
g256_46_bam_file50_g256_45 = Channel.empty()
g256_46_bam_index31_g256_44 = Channel.empty()
g256_46_bam_index61_g256_45 = Channel.empty()
g256_46_filter42_g256_26 = Channel.empty()
g256_46_log_file70_g256_30 = Channel.empty()
g256_46_bigWig_file88 = Channel.empty()
} else {


process Sequential_Mapping_Module_Sequential_Mapping {

publishDir params.outdir, mode: 'copy', saveAs: {filename -> if (filename =~ /.*\/.*_sorted.bam$/) "sequential_mapping/$filename"}
publishDir params.outdir, mode: 'copy', saveAs: {filename -> if (filename =~ /.*\/.*_sorted.bam.bai$/) "sequential_mapping/$filename"}
publishDir params.outdir, mode: 'copy', saveAs: {filename -> if (filename =~ /.*\/.*_duplicates_stats.log$/) "sequential_mapping/$filename"}
publishDir params.outdir, mode: 'copy', saveAs: {filename -> if (filename =~ /.*\/.*.bw$/) "sequential_mapping/$filename"}
input:
 set val(name), file(reads) from g257_20_reads00_g256_46
 val mate from g_229_mate1_g256_46
 file commondb from g256_43_commondb02_g256_46
 file bowtie_index from g256_43_bowtieIndex13_g256_46
 file bowtie2_index from g256_43_bowtie2index24_g256_46
 file star_index from g256_43_starIndex35_g256_46

output:
 set val(name), file("final_reads/*q")  into g256_46_reads01_g_127, g256_46_reads01_g248_36, g256_46_reads01_g250_26
 set val(name), file("bowfiles/?*") optional true  into g256_46_bowfiles10_g256_26, g256_46_bowfiles15_g_177
 file "*/*_sorted.bam" optional true  into g256_46_bam_file20_g256_44
 file "*/*_sorted.bam.bai" optional true  into g256_46_bam_index31_g256_44
 val filtersList  into g256_46_filter42_g256_26
 file "*/*_sorted.dedup.bam" optional true  into g256_46_bam_file50_g256_45
 file "*/*_sorted.dedup.bam.bai" optional true  into g256_46_bam_index61_g256_45
 file "*/*_duplicates_stats.log" optional true  into g256_46_log_file70_g256_30
 file "*/*.bw" optional true  into g256_46_bigWig_file88

errorStrategy 'retry'

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
create_bigWig = params.Sequential_Mapping_Module_Sequential_Mapping.create_bigWig
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
    IFS=',' read -r -a selectSeqListAr <<< "${selectSequencesList}"
    IFS=',' read -r -a paramsListAr <<< "${paramsList}" #create comma separated array 
    IFS=',' read -r -a filtersListAr <<< "${filtersList}"
    IFS=',' read -r -a indexesListAr <<< "${indexesList}"
    IFS=',' read -r -a alignersListAr <<< "${alignersList}"
    IFS=',' read -r -a senseListAr <<< "${senseList}"
    wrkDir=\$(pwd)
    startDir=\$(pwd)
    for rna_set in ${mappingList}
    do
        ((k++))
        printf -v k2 "%02d" "\$k" #turn into two digit format
        mkdir -p \${rna_set}/unmapped
        cd \$rna_set
        ## create link of the target file to prevent "too many symlinks error"
        for r in \${startDir}/\${prev}/*; do
            targetRead=\$(readlink -e \$r)
            rname=\$(basename \$r)
            echo "INFO: ln -s \$targetRead \$rname"
            ln -s \$targetRead \$rname
        done
        basename=""
        genomeDir="\${startDir}/\${indexesListAr[\$k-1]}"
        
        if [ "\${selectSeqListAr[\$k-1]}" == "genome" ]; then
        	wrkDir="\${startDir}"
        	if [ "\${alignersListAr[\$k-1]}" == "bowtie" ]; then
        		basename="/\$(basename \${startDir}/\${indexesListAr[\$k-1]}/*.rev.1.ebwt | cut -d. -f1)"
        	elif [ "\${alignersListAr[\$k-1]}" == "bowtie2" ]; then
        		basename="/\$(basename \${startDir}/\${indexesListAr[\$k-1]}/*.rev.1.bt2 | cut -d. -f1)"
        	elif [ "\${alignersListAr[\$k-1]}" == "STAR" ]; then
        		basename="/\$(basename \${startDir}/\${indexesListAr[\$k-1]}/*.gtf | cut -d. -f1)"
        	fi
        elif [ "\${selectSeqListAr[\$k-1]}" == "custom" ] ; then
        	wrkDir=""
        fi
        echo "INFO: basename: \$basename"
        echo "INFO: genomeDir: \$genomeDir"
        echo "INFO: check bowtie index: \${wrkDir}/\${indexesListAr[\$k-1]}\${basename}.rev.1.ebwt"
        echo "INFO: check bowtie2 index: \${wrkDir}/\${indexesListAr[\$k-1]}\${basename}.1.bt2"
        echo "INFO: check star index: \${genomeDir}/SAindex"
        if [ -e "\${wrkDir}/\${indexesListAr[\$k-1]}\${basename}.rev.1.ebwt" -o -e "\${wrkDir}/\${indexesListAr[\$k-1]}\${basename}.1.bt2" -o  -e "\${wrkDir}/\${indexesListAr[\$k-1]}\${basename}.fa"  -o  -e "\${wrkDir}/\${indexesListAr[\$k-1]}\${basename}.fasta"  -o  -e "\${genomeDir}/SAindex" ]; then
            if [ -e "\${wrkDir}/\${indexesListAr[\$k-1]}\${basename}.fa" ] ; then
                fasta=\${wrkDir}/\${indexesListAr[\$k-1]}\${basename}.fa
            elif [ -e "\${wrkDir}/\${indexesListAr[\$k-1]}\${basename}.fasta" ] ; then
                fasta=\${wrkDir}/\${indexesListAr[\$k-1]}\${basename}.fasta
            fi
            echo "INFO: fasta: \$fasta"
            if [ -e "\${wrkDir}/\${indexesListAr[\$k-1]}\${basename}.1.bt2" -a "\${alignersListAr[\$k-1]}" == "bowtie2" ] ; then
                echo "INFO: \${wrkDir}/\${indexesListAr[\$k-1]}\${basename}.1.bt2 Bowtie2 index found."
            elif [ -e "\${wrkDir}/\${indexesListAr[\$k-1]}\${basename}.1.ebwt" -a "\${alignersListAr[\$k-1]}" == "bowtie" ] ; then
                echo "INFO: \${wrkDir}/\${indexesListAr[\$k-1]}\${basename}.1.ebwt Bowtie index found."
            elif [ -e "\$genomeDir/SAindex" -a "\${alignersListAr[\$k-1]}" == "STAR" ] ; then
                echo "INFO: \$genomeDir/SAindex STAR index found."
            elif [ -e "\${wrkDir}/\${indexesListAr[\$k-1]}\${basename}.fa" -o  -e "\${wrkDir}/\${indexesListAr[\$k-1]}\${basename}.fasta" ] ; then
                if [ "\${alignersListAr[\$k-1]}" == "bowtie2" ]; then
                    bowtie2-build \$fasta \${wrkDir}/\${indexesListAr[\$k-1]}\${basename}
                elif [ "\${alignersListAr[\$k-1]}" == "STAR" ]; then
                    if [ -e "\${wrkDir}/\${indexesListAr[\$k-1]}\${basename}.gtf" ]; then
                        STAR --runMode genomeGenerate --genomeDir \$genomeDir --genomeFastaFiles \$fasta --sjdbGTFfile \${wrkDir}/\${indexesListAr[\$k-1]}\${basename}.gtf --genomeSAindexNbases 5
                    else
                        echo "WARNING: \${wrkDir}/\${indexesListAr[\$k-1]}\${basename}.gtf not found. STAR index is not generated."
                    fi
                elif [ "\${alignersListAr[\$k-1]}" == "bowtie" ]; then
                    bowtie-build \$fasta \${wrkDir}/\${indexesListAr[\$k-1]}\${basename}
                fi
            fi
                
            if [ "${mate}" == "pair" ]; then
                if [ "\${alignersListAr[\$k-1]}" == "bowtie2" ]; then
                    bowtie2 \${paramsListAr[\$k-1]} -x \${wrkDir}/\${indexesListAr[\$k-1]}\${basename} --no-unal --un-conc unmapped/${name}.unmapped.fastq -1 ${name}.1.fastq -2 ${name}.2.fastq --al-conc ${name}.fq.mapped -S \${rna_set}_${name}_alignment.sam 2>&1 | tee \${k2}_${name}.bow_\${rna_set}
                elif [ "\${alignersListAr[\$k-1]}" == "STAR" ]; then
                    STAR \${paramsListAr[\$k-1]}  --genomeDir \$genomeDir --readFilesIn ${name}.1.fastq ${name}.2.fastq --outSAMtype SAM  --outFileNamePrefix ${name}.star --outReadsUnmapped Fastx
                    mv ${name}.starAligned.out.sam \${rna_set}_${name}_alignment.sam
                    mv ${name}.starUnmapped.out.mate1 unmapped/${name}.unmapped.1.fastq
                    mv ${name}.starUnmapped.out.mate2 unmapped/${name}.unmapped.2.fastq
                    mv ${name}.starLog.final.out \${k2}_${name}.star_\${rna_set}
                elif [ "\${alignersListAr[\$k-1]}" == "bowtie" ]; then
                    bowtie \${paramsListAr[\$k-1]}   \${wrkDir}/\${indexesListAr[\$k-1]}\${basename}  --un  unmapped/${name}.unmapped.fastq -1 ${name}.1.fastq -2 ${name}.2.fastq -S  \${rna_set}_${name}_alignment.sam 2>&1 | tee \${k2}_${name}.bow1_\${rna_set}  
                    mv unmapped/${name}.unmapped_1.fastq unmapped/${name}.unmapped.1.fastq
                    mv unmapped/${name}.unmapped_2.fastq unmapped/${name}.unmapped.2.fastq
                fi
            else
                if [ "\${alignersListAr[\$k-1]}" == "bowtie2" ]; then
                    bowtie2 \${paramsListAr[\$k-1]} -x \${wrkDir}/\${indexesListAr[\$k-1]}\${basename} --no-unal --un  unmapped/${name}.unmapped.fastq -U ${name}.fastq --al ${name}.fq.mapped -S \${rna_set}_${name}_alignment.sam 2>&1 | tee \${k2}_${name}.bow_\${rna_set}  
                elif [ "\${alignersListAr[\$k-1]}" == "STAR" ]; then
                    STAR \${paramsListAr[\$k-1]}  --genomeDir \$genomeDir --readFilesIn ${name}.fastq --outSAMtype SAM  --outFileNamePrefix ${name}.star --outReadsUnmapped Fastx
                    mv ${name}.starAligned.out.sam \${rna_set}_${name}_alignment.sam
                    mv ${name}.starUnmapped.out.mate1 unmapped/${name}.unmapped.fastq
                    mv ${name}.starLog.final.out \${k2}_${name}.star_\${rna_set}
                elif [ "\${alignersListAr[\$k-1]}" == "bowtie" ]; then
                    bowtie \${paramsListAr[\$k-1]}  \${wrkDir}/\${indexesListAr[\$k-1]}\${basename}  --un  unmapped/${name}.unmapped.fastq  ${name}.fastq  -S \${rna_set}_${name}_alignment.sam 2>&1 | tee \${k2}_${name}.bow1_\${rna_set}  
                    
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
            
            if [ "${create_bigWig}" == "yes" ]; then
				echo "INFO: creating genome.sizes file"
				cat \$fasta | awk '\$0 ~ ">" {print c; c=0;printf substr(\$0,2,100) "\\t"; } \$0 !~ ">" {c+=length(\$0);} END { print c; }' > \${rna_set}.chrom.sizes && sed -i '1{/^\$/d}' \${rna_set}.chrom.sizes
				echo "INFO: creating bigWig file"
				bedtools genomecov -split -bg -ibam \${rna_set}@${name}_sorted.bam -g \${rna_set}.chrom.sizes > \${rna_set}@${name}.bg && wigToBigWig -clip -itemsPerSlot=1 \${rna_set}@${name}.bg \${rna_set}.chrom.sizes \${rna_set}@${name}.bw 
			fi
			
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
                if [ "${create_bigWig}" == "yes" ]; then
					echo "INFO: creating bigWig file for sense antisense bam"
					bedtools genomecov -split -bg -ibam \${rna_set}@${name}_sense_sorted.bam -g \${rna_set}.chrom.sizes > \${rna_set}@${name}_sense.bg && wigToBigWig -clip -itemsPerSlot=1 \${rna_set}@${name}_sense.bg \${rna_set}.chrom.sizes \${rna_set}@${name}_sense.bw 
					bedtools genomecov -split -bg -ibam \${rna_set}@${name}_antisense_sorted.bam -g \${rna_set}.chrom.sizes > \${rna_set}@${name}_antisense.bg && wigToBigWig -clip -itemsPerSlot=1 \${rna_set}@${name}_antisense.bg \${rna_set}.chrom.sizes \${rna_set}@${name}_antisense.bw 
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
                if [ "${create_bigWig}" == "yes" ]; then
					echo "INFO: creating bigWig file for dedup bam"
					bedtools genomecov -split -bg -ibam \${rna_set}@${name}_sorted.dedup.bam -g \${rna_set}.chrom.sizes > \${rna_set}@${name}_dedup.bg && wigToBigWig -clip -itemsPerSlot=1 \${rna_set}@${name}_dedup.bg \${rna_set}.chrom.sizes \${rna_set}@${name}_dedup.bw 
				fi
                
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
	                if [ "${create_bigWig}" == "yes" ]; then
						echo "INFO: creating bigWig file for sense antisense bam"
						bedtools genomecov -split -bg -ibam \${rna_set}@${name}_sense_sorted.dedup.bam -g \${rna_set}.chrom.sizes > \${rna_set}@${name}_sense_sorted.dedup.bg && wigToBigWig -clip -itemsPerSlot=1 \${rna_set}@${name}_sense_sorted.dedup.bg \${rna_set}.chrom.sizes \${rna_set}@${name}_sense_sorted.dedup.bw
						bedtools genomecov -split -bg -ibam \${rna_set}@${name}_antisense_sorted.dedup.bam -g \${rna_set}.chrom.sizes > \${rna_set}@${name}_antisense_sorted.dedup.bg && wigToBigWig -clip -itemsPerSlot=1 \${rna_set}@${name}_antisense_sorted.dedup.bg \${rna_set}.chrom.sizes \${rna_set}@${name}_antisense_sorted.dedup.bw
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
            echo "WARNING: \${startDir}/\${indexesListAr[\$k-1]}\${basename} Mapping skipped. File not found."
            cd unmapped 
            ln -s \${startDir}/\${rna_set}/*fastq .
            cd ..
            cd ..
        fi
        
        if [ "\${filtersListAr[\$k-1]}" == "Yes" ]; then
            prev=\${rna_set}/unmapped
        fi
    done
    cd final_reads && ln -s \${startDir}/\${prev}/* .
else 
    mv ${reads} final_reads/.
fi
"""

}
}


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

publishDir params.outdir, mode: 'copy', saveAs: {filename -> if (filename =~ /pipe.rsem.${name}$/) "rsem/$filename"}
input:
 val mate from g_229_mate0_g250_26
 set val(name), file(reads) from g256_46_reads01_g250_26
 file rsemIndex from g250_32_rsemIndex02_g250_26

output:
 file "pipe.rsem.${name}"  into g250_26_rsemOut00_g250_17, g250_26_rsemOut00_g250_21, g250_26_rsemOut01_g_177
 set val(name), file("pipe.rsem.*/*.genome.bam") optional true  into g250_26_bam_file11
 set val(name), file("pipe.rsem.*/*.bam") optional true  into g250_26_mapped_reads22
 file "pipe.rsem.*/*.genome.bam" optional true  into g250_26_genome_bam30_g250_23

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


refType = ""
if (RSEM_reference_type == "star"){
    refType = "--star"
} else if (RSEM_reference_type == "bowtie2"){
    refType = "--bowtie2"
} else if (RSEM_reference_type == "bowtie"){
    refType = ""
}
"""
basename=\$(basename ${rsemIndex}/*.ti | cut -d. -f1)
rsemRef=${rsemIndex}/\${basename}
$runGzip
mkdir -p pipe.rsem.${name}

if [ "${mate}" == "pair" ]; then
    echo "rsem-calculate-expression ${refType} ${RSEM_parameters} ${genome_BamText} ${noBamText} --paired-end ${file1} ${file2} \${rsemRef} pipe.rsem.${name}/rsem.out.${name}"
    rsem-calculate-expression ${refType} ${RSEM_parameters} ${genome_BamText} ${noBamText} --paired-end ${file1} ${file2} \${rsemRef} pipe.rsem.${name}/rsem.out.${name}
	if [ "${sense_antisense}" == "Yes" ]; then
		 rsem-calculate-expression ${refType} ${RSEM_parameters} ${genome_BamText} ${noBamText} --forward-prob 1 --paired-end ${file1} ${file2} \${rsemRef} pipe.rsem.${name}/rsem.out.forward.${name}
		 rsem-calculate-expression ${refType} ${RSEM_parameters} ${genome_BamText} ${noBamText} --forward-prob 0 --paired-end ${file1} ${file2} \${rsemRef} pipe.rsem.${name}/rsem.out.reverse.${name}
	fi
else
    echo "rsem-calculate-expression ${refType} ${RSEM_parameters} ${genome_BamText} ${noBamText} --calc-ci ${file1} \${rsemRef} pipe.rsem.${name}/rsem.out.${name}"
    rsem-calculate-expression ${refType} ${RSEM_parameters} ${genome_BamText} ${noBamText} --calc-ci ${file1} \${rsemRef} pipe.rsem.${name}/rsem.out.${name}
	if [ "${sense_antisense}" == "Yes" ]; then
	    rsem-calculate-expression ${refType} ${RSEM_parameters} ${genome_BamText} ${noBamText} --forward-prob 1 --calc-ci ${file1} \${rsemRef} pipe.rsem.${name}/rsem.out.forward.${name}
    	rsem-calculate-expression ${refType} ${RSEM_parameters} ${genome_BamText} ${noBamText} --forward-prob 0 --calc-ci ${file1} \${rsemRef} pipe.rsem.${name}/rsem.out.reverse.${name}
	fi
fi

if [ -e "pipe.rsem.${name}/rsem.out.${name}.genome.bam" ] ; then
    mv pipe.rsem.${name}/rsem.out.${name}.genome.bam pipe.rsem.${name}/rsem.out.${name}.genome.unsorted.bam
    samtools sort -o pipe.rsem.${name}/rsem.out.${name}.genome.bam pipe.rsem.${name}/rsem.out.${name}.genome.unsorted.bam
    rm pipe.rsem.${name}/rsem.out.${name}.genome.unsorted.bam
fi
"""

}


process RSEM_module_file_to_set_conversion_for_bam {

input:
 file bam from g250_26_genome_bam30_g250_23.flatten()

output:
 set val(name),file(bam)  into g250_23_bam_file00_g251_121, g250_23_bam_file01_g251_131, g250_23_bam_file00_g251_133, g250_23_bam_file00_g251_134, g250_23_bam_file00_g251_142

script:
name = bam.baseName
"""
echo "done"	
"""
}



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

process BAM_Analysis_Module_RSEM_UCSC_BAM2BigWig_converter {

publishDir params.outdir, mode: 'copy', saveAs: {filename -> if (filename =~ /.*.bw$/) "bigwig_rsem/$filename"}
input:
 set val(name), file(bam) from g250_23_bam_file00_g251_142
 file genomeSizes from g245_54_genomeSizes21_g251_142

output:
 file "*.bw"  into g251_142_outputFileBw00

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
bedtools genomecov -split -bg -ibam ${nameFinal} -g ${genomeSizes} > ${name}.bg 
wigToBigWig -clip -itemsPerSlot=1 ${name}.bg ${genomeSizes} ${name}.bw 
"""
}

//* params.bed =  ""  //* @input

process BAM_Analysis_Module_RSEM_RSeQC {

publishDir params.outdir, mode: 'copy', saveAs: {filename -> if (filename =~ /result\/.*.out$/) "rseqc_rsem/$filename"}
input:
 set val(name), file(bam) from g250_23_bam_file00_g251_134
 file bed from g245_54_bed31_g251_134

output:
 file "result/*.out"  into g251_134_outputFileOut00_g251_95, g251_134_outputFileOut07_g_177

when:
(params.run_RSeQC && (params.run_RSeQC == "yes")) || !params.run_RSeQC

script:
"""
mkdir result
read_distribution.py  -i ${bam} -r ${bed}> result/RSeQC.${name}.out
"""
}


process BAM_Analysis_Module_RSEM_RSeQC_Summary {

publishDir params.outdir, mode: 'copy', saveAs: {filename -> if (filename =~ /.*.tsv$/) "rseqc_summary_rsem/$filename"}
input:
 file rseqcOut from g251_134_outputFileOut00_g251_95.collect()
 val mate from g_229_mate1_g251_95

output:
 file "*.tsv"  into g251_95_outputFileTSV00

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

//* params.gtf =  ""  //* @input


process BAM_Analysis_Module_RSEM_featureCounts {

publishDir params.outdir, mode: 'copy', saveAs: {filename -> if (filename =~ /.*$/) "featureCounts_after_rsem/$filename"}
input:
 set val(name), file(bam) from g250_23_bam_file00_g251_133
 val paired from g_229_mate1_g251_133
 each run_params from g251_125_run_parameters02_g251_133
 file gtf from g245_54_gtfFile03_g251_133

output:
 file "*"  into g251_133_outputFileTSV00_g251_117

script:
pairText = ""
if (paired == "pair"){
    pairText = "-p"
}

run_name = run_params["run_name"] 
run_parameters = run_params["run_parameters"] 

"""
featureCounts ${pairText} ${run_parameters} -a ${gtf} -o ${name}@${run_name}@fCounts.txt ${bam}
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

process BAM_Analysis_Module_RSEM_summary_featureCounts {

publishDir params.outdir, mode: 'copy', saveAs: {filename -> if (filename =~ /.*_featureCounts.tsv$/) "featureCounts_after_rsem_summary/$filename"}
publishDir params.outdir, mode: 'copy', saveAs: {filename -> if (filename =~ /.*_featureCounts.sum.tsv$/) "featureCounts_after_rsem_details/$filename"}
input:
 file featureCountsOut from g251_133_outputFileTSV00_g251_117.collect()

output:
 file "*_featureCounts.tsv"  into g251_117_outputFile00
 file "*_featureCounts.sum.tsv"  into g251_117_outFileTSV11

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

igv_extention_factor = params.BAM_Analysis_Module_RSEM_IGV_BAM2TDF_converter.igv_extention_factor
igv_window_size = params.BAM_Analysis_Module_RSEM_IGV_BAM2TDF_converter.igv_window_size

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

process BAM_Analysis_Module_RSEM_IGV_BAM2TDF_converter {

publishDir params.outdir, mode: 'copy', saveAs: {filename -> if (filename =~ /.*.tdf$/) "igvtools_rsem/$filename"}
input:
 val mate from g_229_mate0_g251_131
 set val(name), file(bam) from g250_23_bam_file01_g251_131
 file genomeSizes from g245_54_genomeSizes22_g251_131

output:
 file "*.tdf"  into g251_131_outputFileOut00

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
igvtools count -w ${igv_window_size} -e ${igv_extention_factor} ${pairedText} ${nameFinal} ${name}.tdf ${genomeSizes}
"""
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

process BAM_Analysis_Module_RSEM_Picard {

input:
 set val(name), file(bam) from g250_23_bam_file00_g251_121

output:
 file "*_metrics"  into g251_121_outputFileOut00_g251_82
 file "results/*.pdf"  into g251_121_outputFilePdf12_g251_82

when:
(params.run_Picard_CollectMultipleMetrics && (params.run_Picard_CollectMultipleMetrics == "yes")) || !params.run_Picard_CollectMultipleMetrics

script:
"""
picard CollectMultipleMetrics OUTPUT=${name}_multiple.out VALIDATION_STRINGENCY=LENIENT INPUT=${bam}
mkdir results && java -jar ${params.pdfbox_path} PDFMerger *.pdf results/${name}_multi_metrics.pdf
"""
}


process BAM_Analysis_Module_RSEM_Picard_Summary {

publishDir params.outdir, mode: 'copy', saveAs: {filename -> if (filename =~ /.*.tsv$/) "picard_summary_rsem/$filename"}
publishDir params.outdir, mode: 'copy', saveAs: {filename -> if (filename =~ /results\/.*.pdf$/) "picard_summary_pdf_rsem/$filename"}
input:
 file picardOut from g251_121_outputFileOut00_g251_82.collect()
 val mate from g_229_mate1_g251_82
 file picardPdf from g251_121_outputFilePdf12_g251_82.collect()

output:
 file "*.tsv"  into g251_82_outputFileTSV00
 file "results/*.pdf"  into g251_82_outputFilePdf11

shell:
'''
#!/usr/bin/env perl
use List::Util qw[min max];
use strict;
use File::Basename;
use Getopt::Long;
use Pod::Usage; 
use Data::Dumper;

runCommand("mkdir results && mv *.pdf results/. ");

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


sub runCommand {
	my ($com) = @_;
	if ($com eq ""){
		return "";
    }
    my $error = system(@_);
	if   ($error) { die "Command failed: $error $com\\n"; }
    else          { print "Command successful: $com\\n"; }
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

publishDir params.outdir, mode: 'copy', saveAs: {filename -> if (filename =~ /.*.tsv$/) "rsem_summary/$filename"}
input:
 file rsemOut from g250_26_rsemOut00_g250_21.collect()

output:
 file "*.tsv"  into g250_21_outputFile00
 val wdir  into g250_21_wdir10_g250_14

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

publishDir params.outdir, mode: 'copy', saveAs: {filename -> if (filename =~ /.*.rmd$/) "rsem_rmarkdown/$filename"}
input:
 val wdir from g250_21_wdir10_g250_14

output:
 file "*.rmd" optional true  into g250_14_rMarkdown00

shell:
cols = params.RSEM_module_CountData_DE.cols
conds = params.RSEM_module_CountData_DE.conds
colsC = convertCommaSepString(cols)
condsC = convertCommaSepString(conds)
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

cols<- c(!{colsC})

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
conds <- factor( c(!{condsC}) )
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

if ("!{cols}" !~ /^ *$/) {
  open OUT, ">rmark.rmd";
  print OUT $script;
  close OUT;
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
 file rsemDir from g250_26_rsemOut00_g250_17.collect()

output:
 file "rsem_alignment_sum.tsv"  into g250_17_outputFileTSV03_g_198

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

publishDir params.outdir, mode: 'copy', saveAs: {filename -> if (filename =~ /kallisto_${name}$/) "kallisto/$filename"}
input:
 val mate from g_229_mate0_g248_36
 set val(name), file(reads) from g256_46_reads01_g248_36
 file kallisto_index from g248_31_kallisto_index02_g248_36
 file gtf from g245_54_gtfFile03_g248_36
 file genome_sizes from g245_54_genomeSizes24_g248_36

output:
 file "kallisto_${name}"  into g248_36_outputDir00_g248_22, g248_36_outputDir00_g248_38, g248_36_outputDir013_g_177
 set val(name), file("*.bam") optional true  into g248_36_bam_file10_g255_121, g248_36_bam_file11_g255_131, g248_36_bam_file10_g255_133, g248_36_bam_file10_g255_134, g248_36_bam_file10_g255_142
 file "*.bam.bai" optional true  into g248_36_bam_bai22

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

genomebamText = (genomebam.toString() != "false") ? "--genomebam --gtf _${gtf} --chromosomes ${genome_sizes}" : ""
fragment_lengthText = (fragment_length.toString() != "" && single_or_paired_end_reads.toString() == "single") ? "-l ${fragment_length}" : ""
standard_deviationText = (standard_deviation.toString() != "" && single_or_paired_end_reads.toString() == "single") ? "-s ${standard_deviation}" : ""
"""
$runGzip
if [[ \$(awk '{print \$3}' ${gtf} | grep -c transcript) -le 1 ]]; then
    echo "transcript entries are not found in gtf file. gffread will add transcript entries."
    gffread -E --keep-genes ${gtf} -T -o- >_${gtf} 2>gffread.log
else
    ln -s ${gtf} _${gtf}
fi


mkdir -p kallisto_${name}
if [ "${mate}" == "pair" ]; then
    kallisto quant ${kallisto_parameters} -i ${kallisto_index} ${genomebamText} -o kallisto_${name} ${file1} ${file2} > kallisto_${name}/kallisto.log 2>&1
else
    kallisto quant --single ${kallisto_parameters} ${fragment_lengthText} ${standard_deviationText}  -i ${kallisto_index} ${genomebamText} -o kallisto_${name} ${file1} > kallisto_${name}/kallisto.log 2>&1
fi


if ls kallisto_${name}/*.bam 1> /dev/null 2>&1; then
    mv kallisto_${name}/*.bam  ${name}.bam
fi
if ls kallisto_${name}/*.bam.bai 1> /dev/null 2>&1; then
    mv kallisto_${name}/*.bam.bai  ${name}.bam.bai
fi

if [ -f kallisto_${name}/abundance.tsv ]; then
   mv kallisto_${name}/abundance.tsv  kallisto_${name}/abundance_isoforms.tsv
fi



"""

}



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

publishDir params.outdir, mode: 'copy', saveAs: {filename -> if (filename =~ /.*.bw$/) "bigwig_kallisto/$filename"}
input:
 set val(name), file(bam) from g248_36_bam_file10_g255_142
 file genomeSizes from g245_54_genomeSizes21_g255_142

output:
 file "*.bw"  into g255_142_outputFileBw00

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
bedtools genomecov -split -bg -ibam ${nameFinal} -g ${genomeSizes} > ${name}.bg 
wigToBigWig -clip -itemsPerSlot=1 ${name}.bg ${genomeSizes} ${name}.bw 
"""
}

//* params.bed =  ""  //* @input

process BAM_Analysis_Module_Kallisto_RSeQC {

publishDir params.outdir, mode: 'copy', saveAs: {filename -> if (filename =~ /result\/.*.out$/) "rseqc_kallisto/$filename"}
input:
 set val(name), file(bam) from g248_36_bam_file10_g255_134
 file bed from g245_54_bed31_g255_134

output:
 file "result/*.out"  into g255_134_outputFileOut00_g255_95, g255_134_outputFileOut014_g_177

when:
(params.run_RSeQC && (params.run_RSeQC == "yes")) || !params.run_RSeQC

script:
"""
mkdir result
read_distribution.py  -i ${bam} -r ${bed}> result/RSeQC.${name}.out
"""
}


process BAM_Analysis_Module_Kallisto_RSeQC_Summary {

publishDir params.outdir, mode: 'copy', saveAs: {filename -> if (filename =~ /.*.tsv$/) "rseqc_summary_kallisto/$filename"}
input:
 file rseqcOut from g255_134_outputFileOut00_g255_95.collect()
 val mate from g_229_mate1_g255_95

output:
 file "*.tsv"  into g255_95_outputFileTSV00

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

//* params.gtf =  ""  //* @input


process BAM_Analysis_Module_Kallisto_featureCounts {

publishDir params.outdir, mode: 'copy', saveAs: {filename -> if (filename =~ /.*$/) "featureCounts_after_kallisto/$filename"}
input:
 set val(name), file(bam) from g248_36_bam_file10_g255_133
 val paired from g_229_mate1_g255_133
 each run_params from g255_125_run_parameters02_g255_133
 file gtf from g245_54_gtfFile03_g255_133

output:
 file "*"  into g255_133_outputFileTSV00_g255_117

script:
pairText = ""
if (paired == "pair"){
    pairText = "-p"
}

run_name = run_params["run_name"] 
run_parameters = run_params["run_parameters"] 

"""
featureCounts ${pairText} ${run_parameters} -a ${gtf} -o ${name}@${run_name}@fCounts.txt ${bam}
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

publishDir params.outdir, mode: 'copy', saveAs: {filename -> if (filename =~ /.*_featureCounts.tsv$/) "featureCounts_after_Kallisto_summary/$filename"}
publishDir params.outdir, mode: 'copy', saveAs: {filename -> if (filename =~ /.*_featureCounts.sum.tsv$/) "featureCounts_after_Kallisto_details/$filename"}
input:
 file featureCountsOut from g255_133_outputFileTSV00_g255_117.collect()

output:
 file "*_featureCounts.tsv"  into g255_117_outputFile00
 file "*_featureCounts.sum.tsv"  into g255_117_outFileTSV11

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

igv_extention_factor = params.BAM_Analysis_Module_Kallisto_IGV_BAM2TDF_converter.igv_extention_factor
igv_window_size = params.BAM_Analysis_Module_Kallisto_IGV_BAM2TDF_converter.igv_window_size

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

publishDir params.outdir, mode: 'copy', saveAs: {filename -> if (filename =~ /.*.tdf$/) "igvtools_kallisto/$filename"}
input:
 val mate from g_229_mate0_g255_131
 set val(name), file(bam) from g248_36_bam_file11_g255_131
 file genomeSizes from g245_54_genomeSizes22_g255_131

output:
 file "*.tdf"  into g255_131_outputFileOut00

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
igvtools count -w ${igv_window_size} -e ${igv_extention_factor} ${pairedText} ${nameFinal} ${name}.tdf ${genomeSizes}
"""
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
 set val(name), file(bam) from g248_36_bam_file10_g255_121

output:
 file "*_metrics"  into g255_121_outputFileOut00_g255_82
 file "results/*.pdf"  into g255_121_outputFilePdf12_g255_82

when:
(params.run_Picard_CollectMultipleMetrics && (params.run_Picard_CollectMultipleMetrics == "yes")) || !params.run_Picard_CollectMultipleMetrics

script:
"""
picard CollectMultipleMetrics OUTPUT=${name}_multiple.out VALIDATION_STRINGENCY=LENIENT INPUT=${bam}
mkdir results && java -jar ${params.pdfbox_path} PDFMerger *.pdf results/${name}_multi_metrics.pdf
"""
}


process BAM_Analysis_Module_Kallisto_Picard_Summary {

publishDir params.outdir, mode: 'copy', saveAs: {filename -> if (filename =~ /.*.tsv$/) "picard_summary_kallisto/$filename"}
publishDir params.outdir, mode: 'copy', saveAs: {filename -> if (filename =~ /results\/.*.pdf$/) "picard_summary_pdf_kallisto/$filename"}
input:
 file picardOut from g255_121_outputFileOut00_g255_82.collect()
 val mate from g_229_mate1_g255_82
 file picardPdf from g255_121_outputFilePdf12_g255_82.collect()

output:
 file "*.tsv"  into g255_82_outputFileTSV00
 file "results/*.pdf"  into g255_82_outputFilePdf11

shell:
'''
#!/usr/bin/env perl
use List::Util qw[min max];
use strict;
use File::Basename;
use Getopt::Long;
use Pod::Usage; 
use Data::Dumper;

runCommand("mkdir results && mv *.pdf results/. ");

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


sub runCommand {
	my ($com) = @_;
	if ($com eq ""){
		return "";
    }
    my $error = system(@_);
	if   ($error) { die "Command failed: $error $com\\n"; }
    else          { print "Command successful: $com\\n"; }
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

process Kallisto_module_Kallisto_transcript_to_gene_count {

input:
 file outDir from g248_36_outputDir00_g248_38
 file gtf from g245_54_gtfFile01_g248_38

output:
 file newoutDir  into g248_38_outputDir00_g248_39

shell:
newoutDir = "genes_" + outDir
'''
#!/usr/bin/env perl
use strict;
use Getopt::Long;
use IO::File;
use Data::Dumper;

my $gtf_file = "!{gtf}";
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

publishDir params.outdir, mode: 'copy', saveAs: {filename -> if (filename =~ /.*.tsv$/) "kallisto_count/$filename"}
input:
 file kallistoOut from g248_38_outputDir00_g248_39.collect()
 file gtf from g245_54_gtfFile01_g248_39

output:
 file "*.tsv"  into g248_39_outputFile00

shell:
'''
#!/usr/bin/env perl
use Data::Dumper;
use strict;

### Parse gtf file
my $gtf_file = "!{gtf}";
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
 file kallistoDir from g248_36_outputDir00_g248_22.collect()

output:
 file "kallisto_alignment_sum.tsv"  into g248_22_outFileTSV09_g_198

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
    my $unique;
    
    # eg. [quant] processed 24,788 reads, 19,238 reads pseudoaligned
	chomp($total   = `cat ${dir}/kallisto.log | grep 'pseudoaligned' | sed 's/,//g' | awk '{sum+=\\$3} END {print sum}'`);
	chomp($aligned = `cat ${dir}/kallisto.log | grep 'pseudoaligned' | sed 's/,//g' | awk '{sum+=\\$5} END {print sum}'`);
	chomp($unique = `cat ${dir}/run_info.json | grep 'n_unique' | sed 's/,//g' | sed 's/"//g' | sed 's/://g' | awk '{sum+=\\$2} END {print sum}'`);
    $tsv{$libname}=[$libname, $total];
    push(@{$tsv{$libname}}, $aligned);
    push(@{$tsv{$libname}}, $unique);
}

push(@headers, "Sample");
push(@headers, "Total Reads");
push(@headers, "Pseudoaligned Reads (Kallisto)");
push(@headers, "Uniquely Mapped Reads (Kallisto)");

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
g256_46_reads01_g_127.into{g_127_reads01_g247_14; g_127_reads01_g249_14; g_127_reads01_g264_31}
} else {


process SplitFastq {

input:
 val mate from g_229_mate0_g_127
 set val(name), file(reads) from g256_46_reads01_g_127.map(flatPairsClosure).splitFastq(splitFastqParams).map(groupPairsClosure)

output:
 set val(name), file("split/*q")  into g_127_reads01_g247_14, g_127_reads01_g249_14, g_127_reads01_g264_31

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
 val mate from g_229_mate0_g249_14
 set val(name), file(reads) from g_127_reads01_g249_14
 file hisat2index from g249_15_hisat2Index02_g249_14

output:
 set val(name), file("${newName}.bam")  into g249_14_mapped_reads00_g249_13
 set val(name), file("${newName}.align_summary.txt")  into g249_14_outputFileTxt10_g249_2
 set val(name), file("${newName}.flagstat.txt")  into g249_14_outputFileOut22

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
basename=\$(basename ${hisat2index}/*.8.ht2 | cut -d. -f1)
$runGzip
if [ "${mate}" == "pair" ]; then
    hisat2 ${HISAT2_parameters} -x ${hisat2index}/\${basename} -1 ${file1} -2 ${file2} -S ${newName}.sam &> ${newName}.align_summary.txt
else
    hisat2 ${HISAT2_parameters} -x ${hisat2index}/\${basename} -U ${file1} -S ${newName}.sam &> ${newName}.align_summary.txt
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

publishDir params.outdir, mode: 'copy', saveAs: {filename -> if (filename =~ /.*_sorted.*bam$/) "hisat2/$filename"}
input:
 set val(oldname), file(bamfiles) from g249_14_mapped_reads00_g249_13.groupTuple()

output:
 set val(oldname), file("${oldname}.bam")  into g249_13_merged_bams00
 set val(oldname), file("*_sorted*bai")  into g249_13_bam_index11
 set val(oldname), file("*_sorted*bam")  into g249_13_sorted_bam20_g252_121, g249_13_sorted_bam21_g252_131, g249_13_sorted_bam20_g252_133, g249_13_sorted_bam20_g252_134, g249_13_sorted_bam20_g252_142

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

process BAM_Analysis_Module_HISAT2_UCSC_BAM2BigWig_converter {

publishDir params.outdir, mode: 'copy', saveAs: {filename -> if (filename =~ /.*.bw$/) "bigwig_hisat2/$filename"}
input:
 set val(name), file(bam) from g249_13_sorted_bam20_g252_142
 file genomeSizes from g245_54_genomeSizes21_g252_142

output:
 file "*.bw"  into g252_142_outputFileBw00

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
bedtools genomecov -split -bg -ibam ${nameFinal} -g ${genomeSizes} > ${name}.bg 
wigToBigWig -clip -itemsPerSlot=1 ${name}.bg ${genomeSizes} ${name}.bw 
"""
}

//* params.bed =  ""  //* @input

process BAM_Analysis_Module_HISAT2_RSeQC {

publishDir params.outdir, mode: 'copy', saveAs: {filename -> if (filename =~ /result\/.*.out$/) "rseqc_hisat2/$filename"}
input:
 set val(name), file(bam) from g249_13_sorted_bam20_g252_134
 file bed from g245_54_bed31_g252_134

output:
 file "result/*.out"  into g252_134_outputFileOut00_g252_95, g252_134_outputFileOut09_g_177

when:
(params.run_RSeQC && (params.run_RSeQC == "yes")) || !params.run_RSeQC

script:
"""
mkdir result
read_distribution.py  -i ${bam} -r ${bed}> result/RSeQC.${name}.out
"""
}


process BAM_Analysis_Module_HISAT2_RSeQC_Summary {

publishDir params.outdir, mode: 'copy', saveAs: {filename -> if (filename =~ /.*.tsv$/) "rseqc_summary_hisat2/$filename"}
input:
 file rseqcOut from g252_134_outputFileOut00_g252_95.collect()
 val mate from g_229_mate1_g252_95

output:
 file "*.tsv"  into g252_95_outputFileTSV00

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

//* params.gtf =  ""  //* @input


process BAM_Analysis_Module_HISAT2_featureCounts {

publishDir params.outdir, mode: 'copy', saveAs: {filename -> if (filename =~ /.*$/) "featureCounts_after_hisat2/$filename"}
input:
 set val(name), file(bam) from g249_13_sorted_bam20_g252_133
 val paired from g_229_mate1_g252_133
 each run_params from g252_125_run_parameters02_g252_133
 file gtf from g245_54_gtfFile03_g252_133

output:
 file "*"  into g252_133_outputFileTSV00_g252_117

script:
pairText = ""
if (paired == "pair"){
    pairText = "-p"
}

run_name = run_params["run_name"] 
run_parameters = run_params["run_parameters"] 

"""
featureCounts ${pairText} ${run_parameters} -a ${gtf} -o ${name}@${run_name}@fCounts.txt ${bam}
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

process BAM_Analysis_Module_HISAT2_summary_featureCounts {

publishDir params.outdir, mode: 'copy', saveAs: {filename -> if (filename =~ /.*_featureCounts.tsv$/) "featureCounts_after_hisat2_summary/$filename"}
publishDir params.outdir, mode: 'copy', saveAs: {filename -> if (filename =~ /.*_featureCounts.sum.tsv$/) "featureCounts_after_hisat2_details/$filename"}
input:
 file featureCountsOut from g252_133_outputFileTSV00_g252_117.collect()

output:
 file "*_featureCounts.tsv"  into g252_117_outputFile00
 file "*_featureCounts.sum.tsv"  into g252_117_outFileTSV11

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

igv_extention_factor = params.BAM_Analysis_Module_HISAT2_IGV_BAM2TDF_converter.igv_extention_factor
igv_window_size = params.BAM_Analysis_Module_HISAT2_IGV_BAM2TDF_converter.igv_window_size

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

process BAM_Analysis_Module_HISAT2_IGV_BAM2TDF_converter {

publishDir params.outdir, mode: 'copy', saveAs: {filename -> if (filename =~ /.*.tdf$/) "igvtools_hisat2/$filename"}
input:
 val mate from g_229_mate0_g252_131
 set val(name), file(bam) from g249_13_sorted_bam21_g252_131
 file genomeSizes from g245_54_genomeSizes22_g252_131

output:
 file "*.tdf"  into g252_131_outputFileOut00

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
igvtools count -w ${igv_window_size} -e ${igv_extention_factor} ${pairedText} ${nameFinal} ${name}.tdf ${genomeSizes}
"""
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

process BAM_Analysis_Module_HISAT2_Picard {

input:
 set val(name), file(bam) from g249_13_sorted_bam20_g252_121

output:
 file "*_metrics"  into g252_121_outputFileOut00_g252_82
 file "results/*.pdf"  into g252_121_outputFilePdf12_g252_82

when:
(params.run_Picard_CollectMultipleMetrics && (params.run_Picard_CollectMultipleMetrics == "yes")) || !params.run_Picard_CollectMultipleMetrics

script:
"""
picard CollectMultipleMetrics OUTPUT=${name}_multiple.out VALIDATION_STRINGENCY=LENIENT INPUT=${bam}
mkdir results && java -jar ${params.pdfbox_path} PDFMerger *.pdf results/${name}_multi_metrics.pdf
"""
}


process BAM_Analysis_Module_HISAT2_Picard_Summary {

publishDir params.outdir, mode: 'copy', saveAs: {filename -> if (filename =~ /.*.tsv$/) "picard_summary_hisat2/$filename"}
publishDir params.outdir, mode: 'copy', saveAs: {filename -> if (filename =~ /results\/.*.pdf$/) "picard_summary_pdf_hisat2/$filename"}
input:
 file picardOut from g252_121_outputFileOut00_g252_82.collect()
 val mate from g_229_mate1_g252_82
 file picardPdf from g252_121_outputFilePdf12_g252_82.collect()

output:
 file "*.tsv"  into g252_82_outputFileTSV00
 file "results/*.pdf"  into g252_82_outputFilePdf11

shell:
'''
#!/usr/bin/env perl
use List::Util qw[min max];
use strict;
use File::Basename;
use Getopt::Long;
use Pod::Usage; 
use Data::Dumper;

runCommand("mkdir results && mv *.pdf results/. ");

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


sub runCommand {
	my ($com) = @_;
	if ($com eq ""){
		return "";
    }
    my $error = system(@_);
	if   ($error) { die "Command failed: $error $com\\n"; }
    else          { print "Command successful: $com\\n"; }
}
'''

}


process HISAT2_Module_HISAT2_Summary {

input:
 set val(name), file(alignSum) from g249_14_outputFileTxt10_g249_2.groupTuple()

output:
 file "*.tsv"  into g249_2_outputFile00_g249_10
 val "hisat2_alignment_sum"  into g249_2_name11_g249_10

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

publishDir params.outdir, mode: 'copy', saveAs: {filename -> if (filename =~ /${name}.tsv$/) "hisat2_summary/$filename"}
input:
 file tsv from g249_2_outputFile00_g249_10.collect()
 val outputFileName from g249_2_name11_g249_10.collect()

output:
 file "${name}.tsv"  into g249_10_outputFileTSV02_g_198

errorStrategy 'retry'
maxRetries 3

script:
name = outputFileName[0]
"""    
awk 'FNR==1 && NR!=1 {  getline; } 1 {print} ' *.tsv > ${name}.tsv
"""
}

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
 val mate from g_229_mate0_g247_14
 set val(name), file(reads) from g_127_reads01_g247_14
 file bowtie2_index from g247_18_bowtie2index02_g247_14
 file gtf from g245_54_gtfFile03_g247_14

output:
 set val(name), file("${newName}.bam")  into g247_14_mapped_reads00_g247_13
 set val(name), file("${newName}_unmapped.bam")  into g247_14_unmapped_reads11
 set val(name), file("${newName}_align_summary.txt")  into g247_14_summary20_g247_3, g247_14_summary20_g_177

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
basename=\$(basename ${bowtie2_index}/*.rev.1.bt2 | cut -d. -f1)
$runGzip
if [ "${mate}" == "pair" ]; then
    tophat2 ${tophat2_parameters}  --keep-tmp -G ${gtf} -o . ${bowtie2_index}/\${basename} $file
else
    tophat2 ${tophat2_parameters}  --keep-tmp -G ${gtf} -o . ${bowtie2_index}/\${basename} $file
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

publishDir params.outdir, mode: 'copy', saveAs: {filename -> if (filename =~ /.*_sorted.*bai$/) "tophat2/$filename"}
publishDir params.outdir, mode: 'copy', saveAs: {filename -> if (filename =~ /.*_sorted.*bam$/) "tophat2/$filename"}
input:
 set val(oldname), file(bamfiles) from g247_14_mapped_reads00_g247_13.groupTuple()

output:
 set val(oldname), file("${oldname}.bam")  into g247_13_merged_bams00
 set val(oldname), file("*_sorted*bai")  into g247_13_bam_index11
 set val(oldname), file("*_sorted*bam")  into g247_13_sorted_bam20_g254_121, g247_13_sorted_bam21_g254_131, g247_13_sorted_bam20_g254_133, g247_13_sorted_bam20_g254_134, g247_13_sorted_bam20_g254_142

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

process BAM_Analysis_Module_Tophat2_UCSC_BAM2BigWig_converter {

publishDir params.outdir, mode: 'copy', saveAs: {filename -> if (filename =~ /.*.bw$/) "bigwig_tophat/$filename"}
input:
 set val(name), file(bam) from g247_13_sorted_bam20_g254_142
 file genomeSizes from g245_54_genomeSizes21_g254_142

output:
 file "*.bw"  into g254_142_outputFileBw00

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
bedtools genomecov -split -bg -ibam ${nameFinal} -g ${genomeSizes} > ${name}.bg 
wigToBigWig -clip -itemsPerSlot=1 ${name}.bg ${genomeSizes} ${name}.bw 
"""
}

//* params.bed =  ""  //* @input

process BAM_Analysis_Module_Tophat2_RSeQC {

publishDir params.outdir, mode: 'copy', saveAs: {filename -> if (filename =~ /result\/.*.out$/) "rseqc_tophat2/$filename"}
input:
 set val(name), file(bam) from g247_13_sorted_bam20_g254_134
 file bed from g245_54_bed31_g254_134

output:
 file "result/*.out"  into g254_134_outputFileOut00_g254_95, g254_134_outputFileOut06_g_177

when:
(params.run_RSeQC && (params.run_RSeQC == "yes")) || !params.run_RSeQC

script:
"""
mkdir result
read_distribution.py  -i ${bam} -r ${bed}> result/RSeQC.${name}.out
"""
}


process BAM_Analysis_Module_Tophat2_RSeQC_Summary {

publishDir params.outdir, mode: 'copy', saveAs: {filename -> if (filename =~ /.*.tsv$/) "rseqc_summary_tophat2/$filename"}
input:
 file rseqcOut from g254_134_outputFileOut00_g254_95.collect()
 val mate from g_229_mate1_g254_95

output:
 file "*.tsv"  into g254_95_outputFileTSV00

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

//* params.gtf =  ""  //* @input


process BAM_Analysis_Module_Tophat2_featureCounts {

publishDir params.outdir, mode: 'copy', saveAs: {filename -> if (filename =~ /.*$/) "featureCounts_after_Tophat2/$filename"}
input:
 set val(name), file(bam) from g247_13_sorted_bam20_g254_133
 val paired from g_229_mate1_g254_133
 each run_params from g254_125_run_parameters02_g254_133
 file gtf from g245_54_gtfFile03_g254_133

output:
 file "*"  into g254_133_outputFileTSV00_g254_117

script:
pairText = ""
if (paired == "pair"){
    pairText = "-p"
}

run_name = run_params["run_name"] 
run_parameters = run_params["run_parameters"] 

"""
featureCounts ${pairText} ${run_parameters} -a ${gtf} -o ${name}@${run_name}@fCounts.txt ${bam}
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

process BAM_Analysis_Module_Tophat2_summary_featureCounts {

publishDir params.outdir, mode: 'copy', saveAs: {filename -> if (filename =~ /.*_featureCounts.tsv$/) "featureCounts_after_Tophat2_summary/$filename"}
publishDir params.outdir, mode: 'copy', saveAs: {filename -> if (filename =~ /.*_featureCounts.sum.tsv$/) "featureCounts_after_Tophat2_details/$filename"}
input:
 file featureCountsOut from g254_133_outputFileTSV00_g254_117.collect()

output:
 file "*_featureCounts.tsv"  into g254_117_outputFile00
 file "*_featureCounts.sum.tsv"  into g254_117_outFileTSV11

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

igv_extention_factor = params.BAM_Analysis_Module_Tophat2_IGV_BAM2TDF_converter.igv_extention_factor
igv_window_size = params.BAM_Analysis_Module_Tophat2_IGV_BAM2TDF_converter.igv_window_size

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

process BAM_Analysis_Module_Tophat2_IGV_BAM2TDF_converter {

publishDir params.outdir, mode: 'copy', saveAs: {filename -> if (filename =~ /.*.tdf$/) "igvtools_tophat2/$filename"}
input:
 val mate from g_229_mate0_g254_131
 set val(name), file(bam) from g247_13_sorted_bam21_g254_131
 file genomeSizes from g245_54_genomeSizes22_g254_131

output:
 file "*.tdf"  into g254_131_outputFileOut00

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
igvtools count -w ${igv_window_size} -e ${igv_extention_factor} ${pairedText} ${nameFinal} ${name}.tdf ${genomeSizes}
"""
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

process BAM_Analysis_Module_Tophat2_Picard {

input:
 set val(name), file(bam) from g247_13_sorted_bam20_g254_121

output:
 file "*_metrics"  into g254_121_outputFileOut00_g254_82
 file "results/*.pdf"  into g254_121_outputFilePdf12_g254_82

when:
(params.run_Picard_CollectMultipleMetrics && (params.run_Picard_CollectMultipleMetrics == "yes")) || !params.run_Picard_CollectMultipleMetrics

script:
"""
picard CollectMultipleMetrics OUTPUT=${name}_multiple.out VALIDATION_STRINGENCY=LENIENT INPUT=${bam}
mkdir results && java -jar ${params.pdfbox_path} PDFMerger *.pdf results/${name}_multi_metrics.pdf
"""
}


process BAM_Analysis_Module_Tophat2_Picard_Summary {

publishDir params.outdir, mode: 'copy', saveAs: {filename -> if (filename =~ /.*.tsv$/) "picard_summary_tophat2/$filename"}
publishDir params.outdir, mode: 'copy', saveAs: {filename -> if (filename =~ /results\/.*.pdf$/) "picard_summary_pdf_tophat2/$filename"}
input:
 file picardOut from g254_121_outputFileOut00_g254_82.collect()
 val mate from g_229_mate1_g254_82
 file picardPdf from g254_121_outputFilePdf12_g254_82.collect()

output:
 file "*.tsv"  into g254_82_outputFileTSV00
 file "results/*.pdf"  into g254_82_outputFilePdf11

shell:
'''
#!/usr/bin/env perl
use List::Util qw[min max];
use strict;
use File::Basename;
use Getopt::Long;
use Pod::Usage; 
use Data::Dumper;

runCommand("mkdir results && mv *.pdf results/. ");

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


sub runCommand {
	my ($com) = @_;
	if ($com eq ""){
		return "";
    }
    my $error = system(@_);
	if   ($error) { die "Command failed: $error $com\\n"; }
    else          { print "Command successful: $com\\n"; }
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
 set val(name), file(alignSum) from g247_14_summary20_g247_3.groupTuple()
 val mate from g_229_mate1_g247_3

output:
 set val(name), file("${name}_tophat_sum.tsv")  into g247_3_report00_g247_9
 val "tophat2_alignment_sum"  into g247_3_name11_g247_9

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
 file tsv from g247_3_report00_g247_9.collect()
 val outputFileName from g247_3_name11_g247_9.collect()

output:
 file "${name}.tsv"  into g247_9_outputFileTSV04_g_198

errorStrategy 'retry'
maxRetries 3

script:
name = outputFileName[0]
"""    
awk 'FNR==1 && NR!=1 {  getline; } 1 {print} ' *.tsv > ${name}.tsv
"""
}

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
 val mate from g_229_mate0_g264_31
 set val(name), file(reads) from g_127_reads01_g264_31
 file star_index from g264_26_starIndex02_g264_31

output:
 set val(name), file("${newName}Log.final.out")  into g264_31_outputFileOut00_g264_18
 set val(name), file("${newName}.flagstat.txt")  into g264_31_outputFileTxt11
 set val(name), file("${newName}Log.out")  into g264_31_logOut21_g264_18
 set val(name), file("${newName}.bam")  into g264_31_mapped_reads30_g264_30
 set val(name), file("${newName}SJ.out.tab")  into g264_31_outputFileTab43_g264_18
 set val(name), file("${newName}Log.progress.out")  into g264_31_progressOut52_g264_18
 set val(name), file("${newName}Aligned.toTranscriptome.out.bam") optional true  into g264_31_transcriptome_bam60_g264_15
 val sense_antisense  into g264_31_sense_antisense71_g264_30

errorStrategy 'retry'

when:
(params.run_STAR && (params.run_STAR == "yes")) || !params.run_STAR

script:
params_STAR = params.STAR_Module_Map_STAR.params_STAR
sense_antisense = params.STAR_Module_Map_STAR.sense_antisense

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
STAR ${params_STAR}  --genomeDir ${star_index} --readFilesIn $file --outFileNamePrefix ${newName}
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
 set val(oldname), file(bamfiles) from g264_31_transcriptome_bam60_g264_15.groupTuple()

output:
 set val(oldname), file("${oldname}.bam")  into g264_15_merged_bams00
 set val(oldname), file("*_sorted*bai")  into g264_15_bam_index11
 set val(oldname), file("*_sorted*bam")  into g264_15_sorted_bam22

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


process Sequential_Mapping_Module_Sequential_Mapping_Bam_Dedup_count {

publishDir params.outdir, mode: 'copy', saveAs: {filename -> if (filename =~ /.*.counts.tsv$/) "Sequential_Mapping_Bam_dedup_count/$filename"}
input:
 file bam from g256_46_bam_file50_g256_45.collect()
 file index from g256_46_bam_index61_g256_45.collect()
 file bowtie_index from g256_43_bowtieIndex12_g256_45
 file bowtie2_index from g256_43_bowtie2index23_g256_45
 file star_index from g256_43_starIndex34_g256_45
 file commondb from g256_43_commondb05_g256_45

output:
 file "*.counts.tsv"  into g256_45_outputFileTSV00

shell:
mappingListQuoteSep = mapList.collect{ '"' + it + '"'}.join(",")
rawIndexList = indexList.collect{ '"' + it + '"'}.join(",")
selectSeqListQuote = selectSequenceList.collect{ '"' + it + '"'}.join(",")
alignerListQuote = alignerList.collect{ '"' + it + '"'}.join(",")
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
my @selectSeqList = (!{selectSeqListQuote});
my @alignerList = (!{alignerListQuote});

my %indexHash;
my %selectSeqHash;
my %alignerHash;
my $dedup = "";
@indexHash{@mappingList} = @rawIndexList;
@selectSeqHash{@mappingList} = @selectSeqList;
@alignerHash{@mappingList} = @alignerList;

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
	   
	   my $prefix = ${indexHash}{$key};
	   my $selectedSeq = ${selectSeqHash}{$key};
	   my $aligner = ${alignerHash}{$key};
	   if ($selectedSeq eq "genome"){
	   	  if ($aligner eq "bowtie"){
	   		$basename = `basename $prefix/*.rev.1.ebwt | cut -d. -f1`;
	   	  } elsif ($aligner eq "bowtie2"){
	   		$basename = `basename $prefix/*.rev.1.bt2 | cut -d. -f1`;
	   	  } elsif ($aligner eq "STAR"){
	   	    $basename = `basename $prefix/*.gtf | cut -d. -f1`;
	   	  }
	   	  $basename =~ s|\\s*$||;
	   	  $prefix = $prefix."/".$basename;
	   }
	   
		unless (-e $prefix.".bed") {
            print "2: bed not found run makeBed\\n";
                if (-e $prefix.".fa") {
                    makeBed($prefix.".fa", $key, $prefix.".bed");
                } elsif(-e $prefix.".fasta"){
                    makeBed($prefix.".fasta", $key, $prefix.".bed");
                }
        }
	    
		my $com =  "bedtools multicov $par -bams $bamFiles -bed ".$prefix.".bed > $key${dedup}${sense_antisense}.counts.tmp\\n";
        print $com;
        `$com`;
        my $iniResColumn = int(countColumn($prefix.".bed")) + 1;
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
    $name=~s/\\r//g;
    print OUT "$name\\t1\\t".length($seq)."\\t$name\\t0\\t+\\n" if (length($seq)>0); 
    close OUT;
}

'''
}


process Sequential_Mapping_Module_Sequential_Mapping_Bam_count {

publishDir params.outdir, mode: 'copy', saveAs: {filename -> if (filename =~ /.*.counts.tsv$/) "Sequential_Mapping_Bam_count/$filename"}
input:
 file bam from g256_46_bam_file20_g256_44.collect()
 file index from g256_46_bam_index31_g256_44.collect()
 file bowtie_index from g256_43_bowtieIndex12_g256_44
 file bowtie2_index from g256_43_bowtie2index23_g256_44
 file star_index from g256_43_starIndex34_g256_44
 file commondb from g256_43_commondb05_g256_44

output:
 file "*.counts.tsv"  into g256_44_outputFileTSV00

shell:
mappingListQuoteSep = mapList.collect{ '"' + it + '"'}.join(",")
rawIndexList = indexList.collect{ '"' + it + '"'}.join(",")
selectSeqListQuote = selectSequenceList.collect{ '"' + it + '"'}.join(",")
alignerListQuote = alignerList.collect{ '"' + it + '"'}.join(",")
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
my @selectSeqList = (!{selectSeqListQuote});
my @alignerList = (!{alignerListQuote});

my %indexHash;
my %selectSeqHash;
my %alignerHash;
my $dedup = "";
@indexHash{@mappingList} = @rawIndexList;
@selectSeqHash{@mappingList} = @selectSeqList;
@alignerHash{@mappingList} = @alignerList;

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
	   
	   my $prefix = ${indexHash}{$key};
	   my $selectedSeq = ${selectSeqHash}{$key};
	   my $aligner = ${alignerHash}{$key};
	   if ($selectedSeq eq "genome"){
	   	  if ($aligner eq "bowtie"){
	   		$basename = `basename $prefix/*.rev.1.ebwt | cut -d. -f1`;
	   	  } elsif ($aligner eq "bowtie2"){
	   		$basename = `basename $prefix/*.rev.1.bt2 | cut -d. -f1`;
	   	  } elsif ($aligner eq "STAR"){
	   	    $basename = `basename $prefix/*.gtf | cut -d. -f1`;
	   	  }
	   	  $basename =~ s|\\s*$||;
	   	  $prefix = $prefix."/".$basename;
	   }
	   
		unless (-e $prefix.".bed") {
            print "2: bed not found run makeBed\\n";
                if (-e $prefix.".fa") {
                    makeBed($prefix.".fa", $key, $prefix.".bed");
                } elsif(-e $prefix.".fasta"){
                    makeBed($prefix.".fasta", $key, $prefix.".bed");
                }
        }
	    
		my $com =  "bedtools multicov $par -bams $bamFiles -bed ".$prefix.".bed > $key${dedup}${sense_antisense}.counts.tmp\\n";
        print $com;
        `$com`;
        my $iniResColumn = int(countColumn($prefix.".bed")) + 1;
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
    $name=~s/\\r//g;
    print OUT "$name\\t1\\t".length($seq)."\\t$name\\t0\\t+\\n" if (length($seq)>0); 
    close OUT;
}

'''
}


process Sequential_Mapping_Module_Deduplication_Summary {

publishDir params.outdir, mode: 'copy', saveAs: {filename -> if (filename =~ /deduplication_summary.tsv$/) "sequential_mapping_summary/$filename"}
input:
 file flagstat from g256_46_log_file70_g256_30.collect()
 val mate from g_229_mate1_g256_30

output:
 file "deduplication_summary.tsv"  into g256_30_outputFileTSV00

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
 set val(name), file(bowfile) from g256_46_bowfiles10_g256_26
 val mate from g_229_mate1_g256_26
 val filtersList from g256_46_filter42_g256_26

output:
 file '*.tsv'  into g256_26_outputFileTSV00_g256_13
 val "sequential_mapping_sum"  into g256_26_name11_g256_13

errorStrategy 'retry'
maxRetries 2

shell:
'''
#!/usr/bin/env perl
open(my \$fh, '>', '!{name}.tsv');
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
 file tsv from g256_26_outputFileTSV00_g256_13.collect()
 val outputFileName from g256_26_name11_g256_13.collect()

output:
 file "${name}.tsv"  into g256_13_outputFileTSV00_g256_14

errorStrategy 'retry'
maxRetries 3

script:
name = outputFileName[0]
"""    
awk 'FNR==1 && NR!=1 {  getline; } 1 {print} ' *.tsv > ${name}.tsv
"""
}


process Sequential_Mapping_Module_Sequential_Mapping_Short_Summary {

publishDir params.outdir, mode: 'copy', saveAs: {filename -> if (filename =~ /sequential_mapping_short_sum.tsv$/) "sequential_mapping_summary/$filename"}
publishDir params.outdir, mode: 'copy', saveAs: {filename -> if (filename =~ /sequential_mapping_detailed_sum.tsv$/) "sequential_mapping_summary/$filename"}
input:
 file mainSum from g256_13_outputFileTSV00_g256_14

output:
 file "sequential_mapping_short_sum.tsv"  into g256_14_outputFileTSV01_g_198
 file "sequential_mapping_detailed_sum.tsv"  into g256_14_outputFile11

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


process Adapter_Trimmer_Quality_Module_Quality_Filtering_Summary {

input:
 file logfile from g257_20_log_file10_g257_16.collect()
 val mate from g_229_mate1_g257_16

output:
 file "quality_filter_summary.tsv"  into g257_16_outputFileTSV07_g_198

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


process Adapter_Trimmer_Quality_Module_Umitools_Summary {

input:
 file logfile from g257_23_log_file10_g257_24.collect()
 val mate from g_229_mate1_g257_24

output:
 file "umitools_summary.tsv"  into g257_24_outputFileTSV08_g_198

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
    if ($file =~ /(.*)\\.umitools\\.log/){
        $file =~ /(.*)\\.umitools\\.log/;
        my $mapper   = "umitools";
        my $name = $1;    ##sample name
        push( @header, $mapper );
        my $in;
        my $out;
        my $dedupout;
        chomp( $in =`cat $file | grep 'INFO Input Reads:' | awk '{sum=\\$6} END {print sum}'` );
        chomp( $out =`cat $file | grep 'INFO Reads output:' | awk '{sum=\\$6} END {print sum}'` );
        my $deduplog = $name.".dedup.log";
        $headerHash{$mapOrder} = $mapper;
        if (-e $deduplog) {
            print "dedup log found\\n";
            chomp( $dedupout =`cat $deduplog | grep '$name' | awk '{sum=\\$3} END {print sum}'` );
            $tsv{$name}{$mapper} = [ $in, $out, $dedupout];
            $headerText{$mapOrder} = [ "Total Reads", "Reads After Umiextract", "Reads After Deduplication" ]; 
        } else {
            $tsv{$name}{$mapper} = [ $in, $out ];
            $headerText{$mapOrder} = [ "Total Reads", "Reads After Umiextract" ]; 
        }
        
        
    }
}

my @mapOrderArray = ( keys %headerHash );
my @sortedOrderArray = sort { $a <=> $b } @mapOrderArray;

my $summary          = "umitools_summary.tsv";
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


process STAR_Module_STAR_Summary {

publishDir params.outdir, mode: 'copy', saveAs: {filename -> if (filename =~ /.*.(out|tab)$/) "star/$filename"}
input:
 set val(name), file(alignSum) from g264_31_outputFileOut00_g264_18.groupTuple()
 set val(name), file(LogOut) from g264_31_logOut21_g264_18.groupTuple()
 set val(name), file(progressOut) from g264_31_progressOut52_g264_18.groupTuple()
 set val(name), file(TabOut) from g264_31_outputFileTab43_g264_18.groupTuple()

output:
 file "*.tsv"  into g264_18_outputFile00_g264_11
 file "*.{out,tab}"  into g264_18_logOut12_g_177
 val "star_alignment_sum"  into g264_18_name21_g264_11

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
my $name = '!{name}';

# merge output files 
runCommand('cat !{alignSum} >'."${name}_Merged_Log.final.out");
runCommand('cat !{LogOut} >'."${name}_Merged_Log.out");
runCommand('cat !{progressOut} >'."${name}_Merged_Log.progress.out");
runCommand('cat !{TabOut} >'."${name}_Merged_SJ.out.tab");

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

sub runCommand {
    my ($com) = @_;
    if ($com eq ""){
		return "";
    }
    my $error = system(@_);
    if   ($error) { die "Command failed: $error $com\\n"; }
    else          { print "Command successful: $com\\n"; }
}
'''

}


process STAR_Module_merge_tsv_files_with_same_header {

input:
 file tsv from g264_18_outputFile00_g264_11.collect()
 val outputFileName from g264_18_name21_g264_11.collect()

output:
 file "${name}.tsv"  into g264_11_outputFileTSV00_g_198

errorStrategy 'retry'
maxRetries 3

script:
name = outputFileName[0]
"""    
awk 'FNR==1 && NR!=1 {  getline; } 1 {print} ' *.tsv > ${name}.tsv
"""
}

g264_11_outputFileTSV00_g_198= g264_11_outputFileTSV00_g_198.ifEmpty([""]) 
g256_14_outputFileTSV01_g_198= g256_14_outputFileTSV01_g_198.ifEmpty([""]) 
g249_10_outputFileTSV02_g_198= g249_10_outputFileTSV02_g_198.ifEmpty([""]) 
g250_17_outputFileTSV03_g_198= g250_17_outputFileTSV03_g_198.ifEmpty([""]) 
g247_9_outputFileTSV04_g_198= g247_9_outputFileTSV04_g_198.ifEmpty([""]) 
g257_11_outputFileTSV05_g_198= g257_11_outputFileTSV05_g_198.ifEmpty([""]) 
g257_21_outputFileTSV06_g_198= g257_21_outputFileTSV06_g_198.ifEmpty([""]) 
g257_16_outputFileTSV07_g_198= g257_16_outputFileTSV07_g_198.ifEmpty([""]) 
g257_24_outputFileTSV08_g_198= g257_24_outputFileTSV08_g_198.ifEmpty([""]) 
g248_22_outFileTSV09_g_198= g248_22_outFileTSV09_g_198.ifEmpty([""]) 

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

publishDir params.outdir, mode: 'copy', saveAs: {filename -> if (filename =~ /overall_summary.tsv$/) "summary/$filename"}
input:
 file starSum from g264_11_outputFileTSV00_g_198
 file sequentialSum from g256_14_outputFileTSV01_g_198
 file hisatSum from g249_10_outputFileTSV02_g_198
 file rsemSum from g250_17_outputFileTSV03_g_198
 file tophatSum from g247_9_outputFileTSV04_g_198
 file adapterSum from g257_11_outputFileTSV05_g_198
 file trimmerSum from g257_21_outputFileTSV06_g_198
 file qualitySum from g257_16_outputFileTSV07_g_198
 file umiSum from g257_24_outputFileTSV08_g_198
 file kallistoSum from g248_22_outFileTSV09_g_198

output:
 file "overall_summary.tsv" optional true  into g_198_outputFileTSV00

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
if (scalar @rawFiles == 0){
    exit;
}
my @files = ();
# order must be in this order for chipseq pipeline: bowtie->dedup
# rsem bam pipeline: dedup->rsem, star->dedup
# riboseq ncRNA_removal->star
my @order = ("adapter_removal","trimmer","quality","extractUMI","extractValid","tRAX","sequential_mapping","ncRNA_removal","bowtie","star","hisat2","tophat2", "dedup","rsem","kallisto","esat","count");
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
    $TIME = 2000
    $CPU  = 1
    $MEMORY = 8
    $QUEUE = "long"
}
//* platform
//* autofill

process STAR_Module_Merge_Bam_and_create_sense_antisense {

publishDir params.outdir, mode: 'copy', saveAs: {filename -> if (filename =~ /.*_sorted.bam.bai$/) "star/$filename"}
input:
 set val(oldname), file(bamfiles) from g264_31_mapped_reads30_g264_30.groupTuple()
 val sense_antisense from g264_31_sense_antisense71_g264_30
 val mate from g_229_mate2_g264_30

output:
 file "*_sorted.bam.bai"  into g264_30_bam_index00
 file "*_sorted.bam"  into g264_30_bamFile10_g264_29

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

if [ "!{sense_antisense}" == "yes" ]; then
    if [ "!{mate}" == "pair" ]; then
    	echo "INFO: paired end sense antisense separation"
    	samtools view -f 65 -b !{oldname}_sorted.bam > !{oldname}_forward.bam
		samtools view -F 16 -b !{oldname}_forward.bam > !{oldname}_sense.bam && samtools sort -o !{oldname}_sense_sorted.bam !{oldname}_sense.bam && samtools index !{oldname}_sense_sorted.bam
		samtools view -f 16 -b !{oldname}_forward.bam > !{oldname}_antisense.bam  && samtools sort -o !{oldname}_antisense_sorted.bam !{oldname}_antisense.bam && samtools index !{oldname}_antisense_sorted.bam
		rm !{oldname}_forward.bam
	else
	    echo "INFO: single end sense antisense separation"
	    samtools view -F 16 -b !{oldname}_sorted.bam >!{oldname}_sense.bam && samtools sort -o !{oldname}_sense_sorted.bam !{oldname}_sense.bam && samtools index !{oldname}_sense_sorted.bam
	    samtools view -f 16 -b !{oldname}_sorted.bam >!{oldname}_antisense.bam  && samtools sort -o !{oldname}_antisense_sorted.bam !{oldname}_antisense.bam && samtools index !{oldname}_antisense_sorted.bam
	fi
fi
'''
}


process STAR_Module_file_to_set_conversion_for_bam {

publishDir params.outdir, mode: 'copy', saveAs: {filename -> if (filename =~ /bam$/) "star/$filename"}
input:
 file bam from g264_30_bamFile10_g264_29.flatten()

output:
 set val(name),file(bam)  into g264_29_bam_file00_g253_121, g264_29_bam_file01_g253_131, g264_29_bam_file00_g253_133, g264_29_bam_file00_g253_134, g264_29_bam_file00_g253_142

script:
name = bam.baseName
"""
echo "done"	
"""
}



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

process BAM_Analysis_Module_STAR_UCSC_BAM2BigWig_converter {

publishDir params.outdir, mode: 'copy', saveAs: {filename -> if (filename =~ /.*.bw$/) "bigwig_star/$filename"}
input:
 set val(name), file(bam) from g264_29_bam_file00_g253_142
 file genomeSizes from g245_54_genomeSizes21_g253_142

output:
 file "*.bw"  into g253_142_outputFileBw00

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
bedtools genomecov -split -bg -ibam ${nameFinal} -g ${genomeSizes} > ${name}.bg 
wigToBigWig -clip -itemsPerSlot=1 ${name}.bg ${genomeSizes} ${name}.bw 
"""
}

//* params.bed =  ""  //* @input

process BAM_Analysis_Module_STAR_RSeQC {

publishDir params.outdir, mode: 'copy', saveAs: {filename -> if (filename =~ /result\/.*.out$/) "rseqc_star/$filename"}
input:
 set val(name), file(bam) from g264_29_bam_file00_g253_134
 file bed from g245_54_bed31_g253_134

output:
 file "result/*.out"  into g253_134_outputFileOut00_g253_95, g253_134_outputFileOut08_g_177

when:
(params.run_RSeQC && (params.run_RSeQC == "yes")) || !params.run_RSeQC

script:
"""
mkdir result
read_distribution.py  -i ${bam} -r ${bed}> result/RSeQC.${name}.out
"""
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

publishDir params.outdir, mode: 'copy', saveAs: {filename -> if (filename =~ /multiqc_report.html$/) "multiqc/$filename"}
input:
 file "tophat/*" from g247_14_summary20_g_177.flatten().toList()
 file "rsem/*" from g250_26_rsemOut01_g_177.flatten().toList()
 file "star/*" from g264_18_logOut12_g_177.flatten().toList()
 file "fastqc/*" from g257_28_FastQCout04_g_177.flatten().toList()
 file "sequential_mapping/*" from g256_46_bowfiles15_g_177.flatten().toList()
 file "rseqc_tophat/*" from g254_134_outputFileOut06_g_177.flatten().toList()
 file "rseqc_rsem/*" from g251_134_outputFileOut07_g_177.flatten().toList()
 file "rseqc_star/*" from g253_134_outputFileOut08_g_177.flatten().toList()
 file "rseqc_hisat/*" from g252_134_outputFileOut09_g_177.flatten().toList()
 file "kallisto/*" from g248_36_outputDir013_g_177.flatten().toList()
 file "rseqc_kallisto/*" from g255_134_outputFileOut014_g_177.flatten().toList()
 file "after_adapter_removal/*" from g257_31_FastQCout015_g_177.flatten().toList()

output:
 file "multiqc_report.html" optional true  into g_177_outputHTML00

errorStrategy 'ignore'

script:
multiqc_parameters = params.MultiQC.multiqc_parameters
"""
multiqc ${multiqc_parameters} -e general_stats -d -dd 2 .
"""

}


process BAM_Analysis_Module_STAR_RSeQC_Summary {

publishDir params.outdir, mode: 'copy', saveAs: {filename -> if (filename =~ /.*.tsv$/) "rseqc_summary_star/$filename"}
input:
 file rseqcOut from g253_134_outputFileOut00_g253_95.collect()
 val mate from g_229_mate1_g253_95

output:
 file "*.tsv"  into g253_95_outputFileTSV00

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

//* params.gtf =  ""  //* @input


process BAM_Analysis_Module_STAR_featureCounts {

publishDir params.outdir, mode: 'copy', saveAs: {filename -> if (filename =~ /.*$/) "featureCounts_after_STAR/$filename"}
input:
 set val(name), file(bam) from g264_29_bam_file00_g253_133
 val paired from g_229_mate1_g253_133
 each run_params from g253_125_run_parameters02_g253_133
 file gtf from g245_54_gtfFile03_g253_133

output:
 file "*"  into g253_133_outputFileTSV00_g253_117

script:
pairText = ""
if (paired == "pair"){
    pairText = "-p"
}

run_name = run_params["run_name"] 
run_parameters = run_params["run_parameters"] 

"""
featureCounts ${pairText} ${run_parameters} -a ${gtf} -o ${name}@${run_name}@fCounts.txt ${bam}
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

process BAM_Analysis_Module_STAR_summary_featureCounts {

publishDir params.outdir, mode: 'copy', saveAs: {filename -> if (filename =~ /.*_featureCounts.tsv$/) "featureCounts_after_STAR_summary/$filename"}
publishDir params.outdir, mode: 'copy', saveAs: {filename -> if (filename =~ /.*_featureCounts.sum.tsv$/) "featureCounts_after_STAR_details/$filename"}
input:
 file featureCountsOut from g253_133_outputFileTSV00_g253_117.collect()

output:
 file "*_featureCounts.tsv"  into g253_117_outputFile00
 file "*_featureCounts.sum.tsv"  into g253_117_outFileTSV11

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

igv_extention_factor = params.BAM_Analysis_Module_STAR_IGV_BAM2TDF_converter.igv_extention_factor
igv_window_size = params.BAM_Analysis_Module_STAR_IGV_BAM2TDF_converter.igv_window_size

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

process BAM_Analysis_Module_STAR_IGV_BAM2TDF_converter {

publishDir params.outdir, mode: 'copy', saveAs: {filename -> if (filename =~ /.*.tdf$/) "igvtools_star/$filename"}
input:
 val mate from g_229_mate0_g253_131
 set val(name), file(bam) from g264_29_bam_file01_g253_131
 file genomeSizes from g245_54_genomeSizes22_g253_131

output:
 file "*.tdf"  into g253_131_outputFileOut00

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
igvtools count -w ${igv_window_size} -e ${igv_extention_factor} ${pairedText} ${nameFinal} ${name}.tdf ${genomeSizes}
"""
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

process BAM_Analysis_Module_STAR_Picard {

input:
 set val(name), file(bam) from g264_29_bam_file00_g253_121

output:
 file "*_metrics"  into g253_121_outputFileOut00_g253_82
 file "results/*.pdf"  into g253_121_outputFilePdf12_g253_82

when:
(params.run_Picard_CollectMultipleMetrics && (params.run_Picard_CollectMultipleMetrics == "yes")) || !params.run_Picard_CollectMultipleMetrics

script:
"""
picard CollectMultipleMetrics OUTPUT=${name}_multiple.out VALIDATION_STRINGENCY=LENIENT INPUT=${bam}
mkdir results && java -jar ${params.pdfbox_path} PDFMerger *.pdf results/${name}_multi_metrics.pdf
"""
}


process BAM_Analysis_Module_STAR_Picard_Summary {

publishDir params.outdir, mode: 'copy', saveAs: {filename -> if (filename =~ /.*.tsv$/) "picard_summary_star/$filename"}
publishDir params.outdir, mode: 'copy', saveAs: {filename -> if (filename =~ /results\/.*.pdf$/) "picard_summary_pdf_star/$filename"}
input:
 file picardOut from g253_121_outputFileOut00_g253_82.collect()
 val mate from g_229_mate1_g253_82
 file picardPdf from g253_121_outputFilePdf12_g253_82.collect()

output:
 file "*.tsv"  into g253_82_outputFileTSV00
 file "results/*.pdf"  into g253_82_outputFilePdf11

shell:
'''
#!/usr/bin/env perl
use List::Util qw[min max];
use strict;
use File::Basename;
use Getopt::Long;
use Pod::Usage; 
use Data::Dumper;

runCommand("mkdir results && mv *.pdf results/. ");

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


sub runCommand {
	my ($com) = @_;
	if ($com eq ""){
		return "";
    }
    my $error = system(@_);
	if   ($error) { die "Command failed: $error $com\\n"; }
    else          { print "Command successful: $com\\n"; }
}
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
