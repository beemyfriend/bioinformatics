module DNA
	include RNA

	#return number of DNA neucleotides in order of ACGT
	def A_C_G_T(string)
		count = { 	"A" => 0,
					"C" => 0,
					"G" => 0,
					"T" => 0}

		string.split("").each { |a| count[a] += 1}

		 "#{count["A"]} #{count["C"]} #{count["G"]} #{count["T"]}"
	end

	#Return reverse complement of a DNA string
	def complement(string)
		string.reverse.gsub(/[ATCG]/, "A" => "T", "T" => "A", "C" => "G", "G" =>"C" )
	end

	#find a DNA stand's transcribed RNA string
	def to_RNA(string)
		 string.gsub("T", "U")
	end

	#find's the frequency of 'c' and 'g' in a given string
	def GC_Content(string)
		gc = string.split("").select { |a| a == "G" or a == "C"}
		(( gc.length).to_f / (string.strip.length).to_f) * 100 
	end

	#return the location of every motif within a dna string. First symbol of dna string is location 1.
	def motif_location(string, motif)
		locations = []
		(0..(string.length - motif.length)).each { |a| string[a...(a+motif.length)] == motif ? locations << (a+1) : next }
		locations.join(" ")
	end

	#intron's are segements of a gene not used for translation into protein
	def remove_intron(string, intron)
		string.gsub(/#{intron}/, "")
	end

	#hamming distance:  distance between two strings having the same length
	#is the minimum number of symbol substitutions required to transform one string into the other 
	def hamm(x, y)
		count = 0
		(0...x.length).each { |a| x[a] != y[a] ? count += 1 : next}
		count
	end

	#transition/transversion ration
	#transition: point mutation of purine for purine or of a prymadine for another prymadine (a to g) or (c to t)
	#transversions: substitution of a purine for a prymadine
	def trans_ratio(x, y)
		itions = 0.0
		versions = 0.0
		
		transitions = {
			"G" => "A",
			"A" => "G",
			"C" => "T",
			"T" => "C"
		}

		(0...x.length).each do |a|
			if x[a] != y[a]
				transitions[x[a]] == y[a] ? itions += 1 : versions += 1
			end
		end

		itions/versions
	end

	#return longest common substring for an array of strings
	def longest_common_substring(array_of_strings)
		long = 1
		
		common = {}
		common_copy = {}

		(0...(array_of_strings[0].length - long)).each do |a|
				array_of_strings.all? {|b| b.include?(array_of_strings[0][a..(a + long)])} ? common[a] = array_of_strings[0][a..(a+long)] : next 
		end

		until common == {}
			long += 1
			common_copy = {}
			common.each {|k,v| common_copy[k] = v}
			common = {}
			common_copy.each do |k,v|
				if (array_of_strings[0][k+long] != nil) and (array_of_strings.all? {|a| a.include?(v+(array_of_strings[0][k+long]))} ) 
					common[k] = v + array_of_strings[0][k+long]
				end
			end
		end 

		common_copy
	end

	#return a single collection of indices in which the symbols of the motif appear as a subsequence of  the main string
	#subsequence is a collection of symbols contained in order(though not necessarily contiguously)
	def spliced_motif_index(main, motif)
		indices = []
		starter = 0

		(0...motif.length).each do |a|
			indices << (main[starter..-1].split("").index(motif[a]).to_i + 1 + starter)
			starter = indices[-1] 
		end
		indices.join(" ")
	end

	#Return every distinct candidate protein string that can be translated from open reading frames of a string
	#A ORF is a substring which starts from the start codon and ends by a stop codon, without any other stop codons inbetween
	def protein_candidates(data)
		complements = [to_RNA(data), to_RNA(complement(data))]
		candidate_protein = []

		complements.each do |i|
			codons = []
			index = []

			i.split("").each_cons(3) {|a|  codons << a}

			codons.map! {|a| a.join("")}

			(0...codons.length).each {|a| codons[a] == "AUG" ? index << a : next}

			index.each {|a|  candidate_protein << rna_to_protein(i[a..-1])}
		end

		candidate_protein.uniq.compact
	end
end

Module Fasta

	#parse through a fasta formatted file
	def parse(fasta_file)
		counter = 0
		matrix = []

		File.open(fasta_file, "r") do |f|
			f.each_line do |a|
				if a.include?(">")
					counter += 1
					matrix[counter-1] = ""
				else
					matrix[counter-1] += a.chomp
				end
			end
		end
		matrix
	end
end

module RNA
	Codon_table = {
	"UUU" => "F", "CUU" => "L", "AUU" => "I", "GUU" => "V",
	"UUC" => "F", "CUC" => "L", "AUC" => "I", "GUC" => "V",
	"UUA" => "L", "CUA" => "L", "AUA" => "I", "GUA" => "V",
	"UUG" => "L", "CUG" => "L", "AUG" => "M", 'GUG' => 'V',
	'UCU' => 'S', 'CCU' => 'P', 'ACU' => 'T', 'GCU' => 'A',
	'UCC' => 'S', "CCC" => "P", 'ACC' => 'T', 'GCC' => 'A',
	'UCA' => 'S', 'CCA' => 'P', 'ACA' => 'T', 'GCA' => 'A',
	"UCG" => 'S', 'CCG' => 'P', 'ACG' => 'T', 'GCG' => 'A',
	'UAU' => 'Y', 'CAU' => 'H', 'AAU' => 'N', 'GAU' => 'D',
	'UAC' => 'Y', 'CAC' => 'H', 'AAC' => 'N', 'GAC' => 'D',
	"UAA" => 'Stop', 'CAA' => 'Q', 'AAA' => 'K', 'GAA' => 'E',
	'UAG' => 'Stop', 'CAG' => 'Q', 'AAG' => 'K', 'GAG' => 'E',
	'UGU' => 'C', 'CGU' => 'R', 'AGU' => 'S', 'GGU' => 'G',
	'UGC' => 'C', 'CGC' => 'R', 'AGC' => 'S', 'GGC' => 'G',
	'UGA' => 'Stop', 'CGA' => 'R', 'AGA' => 'R', 'GGA' => 'G',
	'UGG' => 'W', 'CGG' => 'R', 'AGG' => 'R', 'GGG' => 'G' 
	}

	#The total number of different RNA strings from which the protein could have been translated, modulo 1,000,000
	def potential_mRNA_strings(string)
		pot_string = {}

		Codon_table.each do |k, v|
			 pot_string.has_key?(v) ? pot_string[v] += 1 : pot_string[v]=1
		end

		mrna = string.split("") << 'Stop'
		answer = mrna.map! {|a| pot_string[a]}.inject(:*)
		modanswer = answer % 1000000
	end

	#return the protein string encoded by the RNA string
	def rna_to_protein(string)
		protein = []
		triplets = string.split("").each_slice(3).to_a.map { |a| a.join }
		
		(0...triplets.length).each do |a|
		 	if (Codon_table[triplets[a-1]] != 'Stop') or (a ==0) 
		 		protein << Codon_table[triplets[a]] 
		 	else
		 		break
		 	end
		 end

		if protein.include?('Stop') 
			protein -= ['Stop']
			protein.join("") 
		end
	end

	Mass_table = {
		'A'  => 71.03711,
		'C'  => 103.00919,
		'D'  => 115.02694,
		'E'  => 129.04259,
		'F'  => 147.06841,
		'G'  => 57.02146,
		'H'  => 137.05891,
		'I'  => 113.08406,
		'K'  => 128.09496,
		'L'  => 113.08406,
		'M'  => 131.04049,
		'N'  => 114.04293,
		'P'  => 97.05276,
		'Q'  => 128.05858,
		'R'  => 156.10111,
		'S'  => 87.03203,
		'T'  => 101.04768,
		'V'  => 99.06841,
		'W'  => 186.07931,
		'Y'  => 163.06333 
	}

	#The total monoisiotopic mass of a protein string
	def monoisotopic_mass(string)
		answer = string.split("").map {|a| Mass_table[a]}.inject(:+)
	end
end

module Tree
	def vertex_edges(array_of_edges)
		vertex_edges = {}
		array_of_edges.each do |a|
			vertex_edges.has_key?(a[0]) ? vertex_edges[a[0]] << a[1] : vertex_edges[a[0]] = [a[1]]
			vertex_edges.has_key?(a[1]) ? vertex_edges[a[1]] << a[0] : vertex_edges[a[1]] = [a[0]]
		end
		vertex_edges
	end

	def leaves(vertex_edges)
		leaves = []
		vertex_edges.each do |k,v|
			if (v.length == 1) and (vertex_edges[k].length >= v.length)
				leaves << k
			end
		end
		leaves
	end

	def prune_leaves(vertex_edges)
		tree = vertex_edges.clone
		leaves = leaves(tree)
		leaves.each do |a|
			tree.delete(a)
			tree.each do |k,v|
				tree[k] = v - [a]
			end
		end
		tree
	end

	def single_edge_trees(vertex_edges)
		tree = vertex_edges.clone
		single_edge_trees = []
		tree.each do |k,v|
			if (tree[k].length == 1) 
				if (tree[v[0]].length == 1)
					single_edge_trees <<  [k,v[0]]
				end
			end
		end
		single_edge_trees
	end

	def min_number_edges(total_vertices)
		total_vertices - 1
	end

	def unrooted_internal_vertices(number_leaves)
		number_leaves - 2
	end
end
