package com.dgraph;

import java.util.ArrayDeque;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import processing.core.PApplet;

public class ClustalwUtil {
	public static void clustalwAnalyze(String[] lines){
		ArrayDeque<Object[]> validLines = new ArrayDeque<Object[]>(100);
		boolean inAlignment = false;
		
		Pattern pairwise_score_re = Pattern.compile("Sequences \\((\\d+):(\\d+)\\) Aligned. Score: ([\\d.]+)");
		for(int k = 0; k < lines.length; k++){
			String line = lines[k];
			if (line.contains("Start of Pairwise alignments")){
				inAlignment = true;
			}
			if (line.contains("Start of Multiple Alignment")){
				inAlignment = false;
			}
			if (!inAlignment){
				continue;
			}
			Matcher matcher = pairwise_score_re.matcher(line);
			if (matcher.matches()) {
				validLines.add(new Object[] {
						Integer.parseInt(matcher.group(1)),
						Integer.parseInt(matcher.group(2)),
						Double.parseDouble(matcher.group(3)),
				});
			}
		}
		int Ncomp = validLines.size();
		int s = (int)((1+PApplet.sqrt(1+8*Ncomp))/2); 
		double[][] toFill = new double[s][s];
		for(int x = 0; x < s; x++){
			toFill[x][x]=0;
		}
		for(Object[] row : validLines) {
			double value = 1/Math.max(0.001f,(Double)row[2] * (Double)row[2]);
			toFill[(Integer)row[0]-1][(Integer)row[1]-1] = value;
			toFill[(Integer)row[1]-1][(Integer)row[0]-1] = value;
		}
		for(int y = 0; y < s; y++){
			for(int x = 0; x < s; x++){
				System.out.print(toFill[y][x]+(x < s - 1 ? "\t" : ""));
			}
			System.out.println();
		}
	}
}
