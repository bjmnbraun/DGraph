package com.dgraph;

import java.util.ArrayDeque;
import java.util.LinkedHashMap;
import java.util.SortedMap;
import java.util.TreeMap;

public class FastaUtil {
	private static LinkedHashMap<String, ArrayDeque<String>> parseFasta(String[] lines){
		LinkedHashMap<String, ArrayDeque<String>> fastaData = new LinkedHashMap<String, ArrayDeque<String>>();
		String currentHeader = null;
		for(String line : lines) {
			if (line.startsWith(">")) {
				currentHeader = line;
			}
			if (currentHeader == null) {
				throw new RuntimeException("Invalid fasta file. First line did not start with >");
			}
			ArrayDeque<String> _fastaData = fastaData.get(currentHeader);
			if (_fastaData == null) {
				_fastaData = new ArrayDeque<String>();
				fastaData.put(line, _fastaData);
			}
			_fastaData.add(line);
		}
		return fastaData;
	}
	private static TreeMap<String, ArrayDeque<String>> parseFastaOrdered(String[] lines){
		TreeMap<String, ArrayDeque<String>> fastaData = new TreeMap<String, ArrayDeque<String>>();
		String currentHeader = null;
		for(String line : lines) {
			if (line.startsWith(">")) {
				currentHeader = line;
			}
			if (currentHeader == null) {
				throw new RuntimeException("Invalid fasta file. First line did not start with >");
			}
			ArrayDeque<String> _fastaData = fastaData.get(currentHeader);
			if (_fastaData == null) {
				_fastaData = new ArrayDeque<String>();
				fastaData.put(line, _fastaData);
			}
			_fastaData.add(line);
		}
		return fastaData;
	}
	public static void reorderFasta(String[] lines, String[] ordering_reference) {
		TreeMap<String, ArrayDeque<String>> data = parseFastaOrdered(lines);
		LinkedHashMap<String, ArrayDeque<String>> reference = parseFasta(ordering_reference);
		//Iteration is over the order of keys in the reference by virtue of LinkedHashMap guarantees
		for(String key : reference.keySet()) {
			String keyPrefix = key.substring(0,Math.min(key.length(),30));
			SortedMap<String, ArrayDeque<String>> tailMap = data.tailMap(keyPrefix);
			String matchingKey = tailMap.firstKey();
			if (!matchingKey.startsWith(keyPrefix)) {
				throw new RuntimeException("Did not find a matching fasta sequence for reference header prefix: "+ keyPrefix);
			}
			//Print the reference header, not the matching in the original.
			System.out.println(key);
			data.get(matchingKey).removeFirst();
			for(String line : data.get(matchingKey)) {
				System.out.println(line);
			}
			data.remove(matchingKey);
		}
	}
}
