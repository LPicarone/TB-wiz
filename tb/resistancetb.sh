#ruby your_ruby_script.rb sample

def leggi_matrice(nome)
  m=[]
  open(nome).readlines.each{ |linea|
    vettore=linea.chomp.split(",")
    vett_interi=[]
    vettore.each{|n| vett_interi.push(n) }
    m.push(vett_interi)
  }
  return m
end


def leggi_vcf(nome)
  m=[]
  open(nome).readlines.each{ |linea|
    vettore=linea.chomp.split("	")
    vett_interi=[]
    vettore.each{|n| vett_interi.push(n) }
    m.push(vett_interi)
  }
  return m
end

def vcf2(nome)
  m=[]
  open(nome).readlines.each{ |linea|
    m.push(linea)
  }
  return m
end


def resistance(nome)
  output_file = "#{nome}.resistance.txt"
  $stdout.reopen(output_file, "w")
annotated=vcf2(nome)
n=0
dbresistance=leggi_matrice("/home/TBwiz/TB-wiz-main/tb/resistancetb.txt")
dbresistance.each{|i|
  gene=i[0]
  mutation=i[1]
  if gene != "-"
    annotated.each { |snp|
      if snp.include?(gene) && snp.include?(mutation)
        array=snp.split("|")
        index=array.index(gene)
        index2=index+6
        index3=index+7
        if array[index]==gene && array[index2]==mutation
          splitting=snp.split(" ")
	  found=splitting[9].split(":")
          puts "Resistance: #{i[2]}"
          puts "Position: #{splitting[1]}"
          puts "Reference allele: #{splitting[3]}"
          puts "Alternate allele: #{splitting[4]}"
	  puts "Depth of variant-supporting bases: #{found[5]}"
          puts "Variant allele frequency #{found[6]} \n"
          puts "Gene: #{i[0]}"
          puts "Mutation:#{i[1]}"
          puts "\n"
          n+=1
        end
      if array[index]==gene && array[index3]==mutation
        splitting=snp.split(" ")
    	found=splitting[9].split(":")
        puts "Resistance: #{i[2]}"
        puts "Position: #{splitting[1]}"
        puts "Reference allele: #{splitting[3]}"
        puts "Alternate allele: #{splitting[4]}"
	puts "Depth of variant-supporting bases: #{found[5]}"
        puts "Variant allele frequency #{found[6]} \n"
        puts "Gene: #{i[0]}"
        puts "Mutation:#{i[1]}"
        puts "\n"
	n+=1
  end


  end
}
  end
}
if n==0
  puts "No resistance found."
end
$stdout.reopen(STDOUT)
end

resistance(ARGV[0])
