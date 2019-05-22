require 'rbbt-util'
require 'rbbt/workflow'

Misc.add_libdir if __FILE__ == $0

#require 'rbbt/sources/MODULE'

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

  dep VEP, :analysis, :args_VEP => "--format vcf --vcf --offline --cache --pick --symbol --terms SO --plugin Downstream --plugin Wildtype --tsl"
  extension :vcf
  task :prepare => :text do
    TSV.traverse step(:analysis), :into => :stream, :type => :array, :bar => "Preparing VEP VCF" do |line|
      if line =~ /^#/
        if line =~ /^#CHR/
          line + "\tFORMAT\tSAMPLE"
        elsif line =~ /WildtypeProtein/
          line.sub!('WildtypeProtein','WildtypeProtein|DownstreamProtein|ProteinLengthChange')
        else
          line
        end
      else
        next unless line =~ /CSQ=/
        #next unless %w(A C T G).include? line.split("\t")[4]
        format = "GT"
        sample = "1/1"
        line = line + "||" 
        line + "\t" + format + "\t" + sample
      end
    end
  end

  extension :vcf
  dep_task :_prepare, VEP, :analysis, :args_VEP => "--pick --symbol --terms SO --plugin Downstream --plugin Wildtype"

  dep :prepare
  input :alleles, :array, "Alleles to query"
  task :analysis => :tsv do |alleles|
    target = file('output')
    alleles = alleles.reject{|a| a =~ /---/}.collect{|a| a =~ /^HLA-/ ? a : "HLA-" << a }.collect{|a| a.split(":").values_at(0,1) * ":"}.uniq

    cpus = config(:cpus, :PVacSeq, :pvacseq, :neo_epitopes, :default => 3)

    CMD.cmd_log("env MHCFLURRY_DOWNLOADS_DIR=#{MHCFLURRY_DOWNLOADS_DIR} \
pvacseq run '#{step(:prepare).produce.path}' Test #{alleles * ","} \
MHCflurry MHCnuggetsI MHCnuggetsII NNalign NetMHC PickPocket SMM SMMPMBEC SMMalign \
'#{target}' \
--iedb-install-directory #{IEDB_INSTALL_DIR}/.. \
-t #{cpus} \
-e 8,9,10 --binding-threshold 1000000")

    files = file('output').glob('MHC_Class_*/Test.filtered.tsv')
    tsv = files.shift.tsv :header_hash => "", :merge => true
    if files.any?
      tsv = tsv.attach files.first.tsv(:header_hash => "", :merge => true)
    end

    tsv.to_s nil, false, true
  end

end

#require 'MODULE/tasks/basic.rb'

#require 'rbbt/knowledge_base/MODULE'
#require 'rbbt/entity/MODULE'

