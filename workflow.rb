require 'rbbt-util'
require 'rbbt/workflow'

Misc.add_libdir if __FILE__ == $0

Workflow.require_workflow "VEP"
module PVacSeq
  extend Workflow

  Rbbt.claim Rbbt.software.opt.IEDB[".source"]["IEDB_MHC_I.tar.gz"], :proc do  |filename|
    raise "Download from http://downloads.iedb.org/ and place in #{filename}"
  end
  Rbbt.claim Rbbt.software.opt.IEDB.mhc_i, :proc do |directory|
    source = Rbbt.software.opt.IEDB[".source"]["IEDB_MHC_I.tar.gz"].produce.find
    Misc.untar(source, directory + '/..')
    Misc.in_dir directory do
      CMD.cmd_log("python src/configure.py")
    end
  end

  MHCFLURRY_DOWNLOADS_DIR=Rbbt.share.databases["mhcflurry"].find
  IEDB_INSTALL_DIR=Rbbt.software.opt.IEDB.mhc_i.produce.find

  extension :vcf
  task :add_sample_info => :text do
    TSV.traverse step(:analysis), :into => :stream, :type => :array, :bar => "Preparing VEP VCF" do |line|
      if line =~ /^#/
        if line =~ /^#CHR/
          line + "\tFORMAT\tSAMPLE"
        else
          line
        end
      else
        next unless line =~ /CSQ=/
        format = "GT"
        sample = "1/1"
        line = line
        line + "\t" + format + "\t" + sample
      end
    end
  end

  dep VEP, :analysis, :args_VEP => "--format vcf --pick --symbol --terms SO --plugin Downstream --plugin Wildtype --tsl", :compute => :produce
  dep :add_sample_info do |jobname,options,dependencies|
    vep = dependencies.flatten.select{|dep| dep.task_name == :analysis && dep.workflow.to_s == "VEP"}.first
    if options[:vcf_file]
      nil
    else
      {:input => options, :jobname => jobname}
    end
  end
  input :alleles, :array, "Alleles to query"
  task :analysis => :tsv do |alleles|
    target = file('output')
    alleles = alleles.reject{|a| a =~ /---/}.collect{|a| a =~ /^HLA-/ ? a : "HLA-" << a }.collect{|a| a.split(":").values_at(0,1) * ":"}.uniq

    cpus = config(:cpus, :PVacSeq, :pvacseq, :neo_epitopes, :default => 3)

    if dependencies.select{|dep| dep.task_name == :add_sample_info}.any?
      file = step(:add_sample_info).path
    else
      file = step(:analysis).path
    end

    samples = TSV.open(file).fields[8..-1]

    if samples.length == 1
      normal_sample = nil
      tumor_sample = samples.first
    else
      normal_sample = samples.select{|s| s =~ /(normal|control|blood|gDNA)/i}.first
      if normal_sample
        tumor_sample = (samples - [normal_sample]).last
      else
        tumor_sample = samples.last
      end
    end

    CMD.cmd_log("env MHCFLURRY_DOWNLOADS_DIR=#{MHCFLURRY_DOWNLOADS_DIR} \
pvacseq run '#{file}' '#{tumor_sample}' #{alleles * ","} \
MHCflurry MHCnuggetsI MHCnuggetsII NNalign NetMHC NetMHCIIpan NetMHCcons NetMHCpan PickPocket SMM SMMPMBEC SMMalign \
'#{target}' \
--iedb-install-directory #{IEDB_INSTALL_DIR}/.. \
#{normal_sample ? "--normal-sample-name #{normal_sample}" : ""} \
-t #{cpus} \
--netmhc-stab \
--pass-only \
-e 8,9,10 --binding-threshold 1000000")

    files = file('output').glob("MHC_Class_*/#{tumor_sample}.filtered.tsv")
    tsv = files.shift.tsv :header_hash => "", :merge => true
    if files.any?
      tsv = tsv.attach files.first.tsv(:header_hash => "", :merge => true)
    end

    tsv.to_s nil, false, true
  end

end
