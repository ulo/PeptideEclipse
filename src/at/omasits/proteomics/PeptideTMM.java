package at.omasits.proteomics;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.zip.GZIPInputStream;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Option;
import org.apache.commons.cli.OptionBuilder;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.apache.commons.cli.PosixParser;

import com.google.common.base.Joiner;
import com.google.common.base.Splitter;
import com.google.common.collect.Lists;

/**
 * PeptideTMM matches proteomics peptide identifications to transmembrane regions as annotated in UniProt.
 * 
 * @author Ulrich Omasits
 * @date 22.07.2013
 */
public class PeptideTMM {
	
	final static Joiner tabJoiner = Joiner.on('\t').useForNull("");
	
	@SuppressWarnings("serial")
	public static Options options = new Options() {{
		addOption( CLIUtils.createArgOption("in", "prot.xml|prot.xls", "ProteinProphet file (can be gzipped)", true, false) );
		addOption( CLIUtils.createArgOption("up", "uniprot.dat", "UniProt knowledgebase (can be gzipped)", true, true) );
	}};
	
	public static void printUsageAndExit() {
		new HelpFormatter().printHelp("java -jar "+PeptideTMM.class.getSimpleName()+".jar", "PeptideTransMembraneMapper by Ulrich Omasits", options, null, true);
		System.exit(0);
	}
	
	public static void main(String[] args) throws FileNotFoundException, IOException {
				
		// parse the command line arguments
		CommandLine cli = null;
		try {
			cli = new PosixParser().parse( options, args );
		} catch (ParseException e) {
			System.err.println(e.getMessage());
			printUsageAndExit();
		}
		
		File fileIn = new File(cli.getOptionValue("in"));
		File[] filesUniProt = CLIUtils.getFileArray(cli, "up");
		
		
		// set up internal data structures
		Map<String, Set<String>> peptidesByProtein = new HashMap<String, Set<String>>();
		Map<String, Map<String, Integer>> tmAAsPerPeptideAndProtein = new HashMap<String, Map<String, Integer>>();
		Map<String, Integer> tmCountPerProtein = new HashMap<String, Integer>();
		Map<String, Integer> tmAAsPerProtein = new HashMap<String, Integer>();
		Map<String, Integer> coverageTmPerProtein = new HashMap<String, Integer>();
		Map<String, Integer> coveragePerProtein = new HashMap<String, Integer>();
		Map<String, Integer> proteinLength = new HashMap<String, Integer>();
		
		
		/*
		 *  read in input file
		 */
		boolean isXML = false;
		System.err.println("reading input file "+fileIn.getName());
		int c_peps = 0;
		int c_entries = 0;
		BufferedReader in = reader(fileIn);
		String line = in.readLine();
		
		if (line.startsWith("<?xml")) {
		// XML input
			isXML = true;
			while((line=in.readLine())!=null) {
				if (line.startsWith("<protein_group ")) {
					while (! (line=in.readLine().trim()).equals("</protein_group>")) { // read to end of protein group
						if (line.startsWith("<protein ")) {
							String prot = extractProteinID(getXMLattribute(line, "protein_name"));
							List<String> peps = Lists.newArrayList(Splitter.on('+').split(getXMLattribute(line, "unique_stripped_peptides")));
							if (!peptidesByProtein.containsKey(prot))
								peptidesByProtein.put(prot, new HashSet<String>());
							for (String pep : peps) {
								c_entries++;
								if (peptidesByProtein.get(prot).add(extractPeptideSequenceIL(pep)))
									c_peps++;
							}
							break; // finish with this protein group after first protein
						}
					}
				}
			}
		} else {
		// XLS (TSV) input
			List<String> header = Lists.newArrayList(Splitter.on('\t').split(line));
			int i_protein = header.indexOf("protein");
			int i_peptide = header.indexOf("peptide sequence");
			while((line=in.readLine())!=null) {
				if (line.length()==0)
					continue;
				c_entries++;
				List<String> elems = Lists.newArrayList(Splitter.on('\t').split(line));
				String prot = extractProteinID( elems.get(i_protein) );
				String pep = extractPeptideSequenceIL( elems.get(i_peptide) );
				if (!peptidesByProtein.containsKey(prot))
					peptidesByProtein.put(prot, new HashSet<String>());
				if (peptidesByProtein.get(prot).add(pep))
					c_peps++;
			}
		}
		in.close();
		System.err.println("loaded "+c_entries+" entries corresponding to a total of "+peptidesByProtein.size()+" proteins with "+c_peps+" peptides");
		
		
		/*
		 * parse swissprot and trembl files and match peptides and topologies
		 */
		int c_prots_found = 0;
		int c_prots_found_TM = 0;
		int c_peps_found = 0;
		int c_peps_found_TM = 0;
		for (File file : filesUniProt) {
			System.err.println("looking for protein sequence and transmembrane regions in "+file.getName());
			in = reader(file);
			while((line=in.readLine())!=null) {
				if (line.startsWith("AC ")) {
					List<String> prots = new ArrayList<String>();
					for (String acc : Splitter.on(';').omitEmptyStrings().split(line.substring(5)))
						if (peptidesByProtein.containsKey(acc)) prots.add(acc);
					if (prots.size()>0) {
						// read in entry
						c_prots_found++;
						StringBuilder seqBuilder = new StringBuilder();
						Set<Integer> tmAAs = new HashSet<Integer>();
						int c_tmRegions = 0;
						while(!(line=in.readLine()).equals("//")) {
							if (line.startsWith("     ")) {
								seqBuilder.append(line.replaceAll(" ", ""));
							} else if (line.startsWith("FT   TRANSMEM ")) {
								c_tmRegions++;
								List<String> tmElems = Lists.newArrayList(Splitter.on(' ').omitEmptyStrings().split(line));
								int from = Integer.parseInt(tmElems.get(2));
								int to = Integer.parseInt(tmElems.get(3));
								for (int i=from; i<=to; i++)
									tmAAs.add(i);
							}
						}
						String seq = seqBuilder.toString().toUpperCase().replace('I', 'L');
						
						// iterate all matched protein accessions by this entry
						for (String prot : prots) {
							tmAAsPerProtein.put(prot, tmAAs.size());
							tmCountPerProtein.put(prot, c_tmRegions);
							if (c_tmRegions > 0)
								c_prots_found_TM++;
							proteinLength.put(prot, seq.length());
							
							// query all peptides
							Set<Integer> pepAAs = new HashSet<Integer>();
							Set<Integer> pepAAsTM = new HashSet<Integer>();
							tmAAsPerPeptideAndProtein.put(prot, new HashMap<String, Integer>());
							for (String pep : peptidesByProtein.get(prot)) {
								int pepFrom = seq.indexOf(pep) + 1;
								if (pepFrom==0) {
									System.err.println("WARN: could not find peptide "+pep+" in protein "+prot+", skipping this peptide...");
									continue;
								}
								c_peps_found++;
								int pepTo = pepFrom + pep.length() - 1;
								int c_pepTmAAs = 0;
								for (int i=pepFrom; i<=pepTo; i++) {
									pepAAs.add(i);
									if (tmAAs.contains(i)) {
										pepAAsTM.add(i);
										c_pepTmAAs++;
									}
								}
								tmAAsPerPeptideAndProtein.get(prot).put(pep, c_pepTmAAs);
								if (c_pepTmAAs>0)
									c_peps_found_TM++;
							}
							
							coveragePerProtein.put(prot, pepAAs.size());
							coverageTmPerProtein.put(prot, pepAAsTM.size());
						}
					}
				}
			}
			in.close();
		}
		System.err.println("could match "+c_prots_found+" proteins ("+c_prots_found_TM+" of them having TM regions) with "+c_peps_found+" peptides ("+c_peps_found_TM+" of them intersecting with a TM region)");
		
		/*
		 * generate output
		 */
		System.err.println("writing output file...");
		in = reader(fileIn);
		if (isXML) {
			println("group_number","group_probability",
					"protein_name","protein_peptides","protein_pct_spectrum_ids","protein_percent_coverage",
					"peptide_sequence","peptide_sequence_modified","peptide_charge","peptide_nsp_adjusted_probability","peptide_n_enzymatic_termini","peptide_calc_neutral_pep_mass",
					"proteinID","protein_transmemRegions","protein_AA","protein_AA_covered","protein_AA_transmem","protein_AA_transmem_covered","peptide_AA","peptide_AA_transmem");
			while((line=in.readLine())!=null) {
				if (line.startsWith("<protein_group ")) {
					List<String> groupElems = Lists.newArrayList(getXMLattribute(line, "group_number"), getXMLattribute(line, "probability"));
					while (! (line=in.readLine().trim()).equals("</protein_group>")) { // read to end of protein group
						if (line.startsWith("<protein ")) {
							List<String> protElems = Lists.newArrayList(getXMLattribute(line, "protein_name"), getXMLattribute(line, "total_number_peptides"),
									getXMLattribute(line, "pct_spectrum_ids"), getXMLattribute(line, "percent_coverage"));
							String prot = extractProteinID(getXMLattribute(line, "protein_name"));
							while (! (line=in.readLine().trim()).equals("</protein>")) { // read to end of protein
								if (line.startsWith("<peptide ")) {
									List<String> pepElems = Lists.newArrayList(getXMLattribute(line, "peptide_sequence"), getXMLattribute(line, "peptide_sequence"),
											getXMLattribute(line, "charge"), getXMLattribute(line, "nsp_adjusted_probability"), 
											getXMLattribute(line, "n_enzymatic_termini"), getXMLattribute(line, "calc_neutral_pep_mass"));
									String pep = extractPeptideSequenceIL(getXMLattribute(line, "peptide_sequence"));
									line = in.readLine().trim();
									if (line.startsWith("<modification_info "))
										pepElems.set(1, getXMLattribute(line, "modified_peptide"));
									
									List<String> elems = new ArrayList<String>();
									elems.addAll(groupElems);
									elems.addAll(protElems);
									elems.addAll(pepElems);
									
									if (proteinLength.containsKey(prot)) {
										elems.add(prot);
										if (tmCountPerProtein.containsKey(prot))
											elems.add(""+tmCountPerProtein.get(prot));
										else
											elems.add("");
										elems.add(""+proteinLength.get(prot));
										elems.add(""+coveragePerProtein.get(prot));				
									} else {
										elems.add("");
										elems.add("");
										elems.add("");
										elems.add("");
									}
									if (tmCountPerProtein.containsKey(prot) && tmCountPerProtein.get(prot)>0) {
										elems.add(""+tmAAsPerProtein.get(prot));
										elems.add(""+coverageTmPerProtein.get(prot));					
									} else {
										elems.add("");
										elems.add("");
									}
									
									if (tmAAsPerPeptideAndProtein.containsKey(prot) && tmAAsPerPeptideAndProtein.get(prot).containsKey(pep)) {
										elems.add(""+pep.length());
										elems.add(""+tmAAsPerPeptideAndProtein.get(prot).get(pep));
									} else {
										elems.add("");
										elems.add("");
									}
									
									println(elems);
								}
							}
							break; // finish with this protein group after first protein
						}
					}
				}
			}
		} else {
			/*
			 *  add columns to each peptide in input file
			 */
			List<String> header = Lists.newArrayList(Splitter.on('\t').split(in.readLine()));
			int i_protein = header.indexOf("protein");
			int i_peptide = header.indexOf("peptide sequence");
			header.remove(header.size()-1);
			header.add("peptide_AA");
			header.add("peptide_AA_transmem");
			header.add("proteinID");
			header.add("protein_transmemRegions");
			header.add("protein_AA");
			header.add("protein_AA_covered");
			header.add("protein_AA_transmem");
			header.add("protein_AA_transmem_covered");
			println(header);
			while((line=in.readLine())!=null) {
				if (line.length()==0)
					continue;
				List<String> elems = Lists.newArrayList(Splitter.on('\t').split(line));
				elems.remove(elems.size()-1); // remove last (empty) columns
				String prot = extractProteinID( elems.get(i_protein) );
				String pep = extractPeptideSequenceIL( elems.get(i_peptide) );
				
				if (tmAAsPerPeptideAndProtein.containsKey(prot) && tmAAsPerPeptideAndProtein.get(prot).containsKey(pep)) {
					elems.add(""+pep.length());
					elems.add(""+tmAAsPerPeptideAndProtein.get(prot).get(pep));
				} else {
					elems.add("");
					elems.add("");
				}
				if (proteinLength.containsKey(prot)) {
					elems.add(prot);
					if (tmCountPerProtein.containsKey(prot))
						elems.add(""+tmCountPerProtein.get(prot));
					else
						elems.add("");
					elems.add(""+proteinLength.get(prot));
					elems.add(""+coveragePerProtein.get(prot));				
				} else {
					elems.add("");
					elems.add("");
					elems.add("");
					elems.add("");
				}
				if (tmCountPerProtein.containsKey(prot) && tmCountPerProtein.get(prot)>0) {
					elems.add(""+tmAAsPerProtein.get(prot));
					elems.add(""+coverageTmPerProtein.get(prot));					
				} else {
					elems.add("");
					elems.add("");
				}
				println(elems);
			}
		}
		in.close();
		System.err.println("Done!");
	}
	
	public static String extractProteinID(String elem) {
		String id = Splitter.on(',').split( elem ).iterator().next();
		if (id.toLowerCase().startsWith("sp|") || id.toLowerCase().startsWith("tr|"))
			id = id.split("\\|")[1];
		if (id.matches(".*\\-\\d+"))
			id = id.substring(0,id.lastIndexOf('-'));
		return id;
	}
	
	public static String extractPeptideSequenceIL(String elem) {
		return elem.toUpperCase().replaceAll("\\[[0123456789\\.]+\\]", "").replace('I', 'L');
	}


	public static void println(Object... values) {
		System.out.println(tabJoiner.join(values));
	}
	
	public static void println(Iterable<?> values) {
		System.out.println(tabJoiner.join(values));
	}
	
	
	public static class CLIUtils {
		public static Option createArgOption(String option, String arg, String description, boolean required, boolean multiple) {
			if (multiple)
				OptionBuilder.hasArgs();
			else
				OptionBuilder.hasArg();
			OptionBuilder.withArgName(arg);
			OptionBuilder.isRequired(required);
			OptionBuilder.withDescription(description);
			return OptionBuilder.create(option);
		}
		
		public static File getFileOption(CommandLine cli, String option, boolean errorOnExistance) {
			if (cli.hasOption(option)) {
				File f = new File(cli.getOptionValue(option));
				if (errorOnExistance && f.exists()) {
					System.err.println("ERROR: " + f.getName() + " already exists.");
					System.exit(0);
				}
				return f;
			}
			return null;
		}
		
		public static File[] getFileArray(CommandLine cli, String option) {
			if (!cli.hasOption(option))
				return null;
			String[] arr = cli.getOptionValues(option);
			File[] files = new File[arr.length];
			for (int i=0; i<arr.length; i++)
				files[i] = new File(arr[i]);
			return files;
		}
		
	}
	
	public static BufferedReader reader(File file) throws FileNotFoundException, IOException {
		if (file.getName().endsWith(".gz") || file.getName().endsWith(".tgz"))
			return new BufferedReader(new InputStreamReader(new GZIPInputStream(new FileInputStream(file))));
		else
			return new BufferedReader(new FileReader(file));
	}
	
	public static String getXMLattribute(String string, String attribute) {
		int begin;
		if ((begin = string.indexOf(attribute+"=\"")) != -1)
			return string.substring(begin+attribute.length()+2, string.indexOf("\"",begin+attribute.length()+2));
		else if ((begin = string.indexOf(attribute+"='")) != -1)
			return string.substring(begin+attribute.length()+2, string.indexOf("'",begin+attribute.length()+2));
		else if ((begin = string.indexOf(attribute+" = \"")) != -1)
			return string.substring(begin+attribute.length()+4, string.indexOf("\"",begin+attribute.length()+4));
		else if ((begin = string.indexOf(attribute+" = '")) != -1)
			return string.substring(begin+attribute.length()+4, string.indexOf("'",begin+attribute.length()+4));
		else if ((begin = string.indexOf(attribute+" =\"")) != -1)
			return string.substring(begin+attribute.length()+3, string.indexOf("\"",begin+attribute.length()+3));
		else if ((begin = string.indexOf(attribute+" ='")) != -1)
			return string.substring(begin+attribute.length()+3, string.indexOf("'",begin+attribute.length()+3));
		else if ((begin = string.indexOf(attribute+"= \"")) != -1)
			return string.substring(begin+attribute.length()+3, string.indexOf("\"",begin+attribute.length()+3));
		else if ((begin = string.indexOf(attribute+"= '")) != -1)
			return string.substring(begin+attribute.length()+3, string.indexOf("'",begin+attribute.length()+3));
		else
			return null;
	}
}
