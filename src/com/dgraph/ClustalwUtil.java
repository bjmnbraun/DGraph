package com.dgraph;

import java.util.ArrayList;

import processing.core.PApplet;

public class ClustalwUtil {
	public static void clustalwAnalyze(String[] lines){
		ArrayList validLines = new ArrayList(30000);
		int[] vals = new int[3];
		boolean inAlignment = false;
		for(int k = 0; k < lines.length; k++){
			if (lines[k].contains("Start of Pairwise alignments")){
				inAlignment = true;
			}
			if (lines[k].contains("Start of Multiple Alignment")){
				inAlignment = false;
			}
			if (!inAlignment){
				continue;
			}
			String[] doLine = lines[k].split("[^\\d]+");
			int got = 0;
			for(int p = 0; p < doLine.length; p++){
				if (!doLine[p].matches("\\s*")){
					vals[got]=new Integer(doLine[p]);
					got++;
					if (got==3){
						got = 0;
						validLines.add(vals);
						vals = new int[3];
					}
				}
			}
		}
		int Ncomp = validLines.size();
		int s = (int)((1+PApplet.sqrt(1+8*Ncomp))/2); 
		double[][] toFill = new double[s][s];
		for(int x = 0; x < s; x++){
			toFill[x][x]=0;
		}
		for(int k = 0; k < Ncomp; k++){
			vals = (int[])validLines.get(k);
			double value = 1/Math.max(0.001f,vals[2] * vals[2]);
			toFill[vals[0]-1][vals[1]-1] = value;
			toFill[vals[1]-1][vals[0]-1] = value;
		}
		for(int y = 0; y < s; y++){
			for(int x = 0; x < s; x++){
				System.out.print(toFill[y][x]+(x < s - 1 ? "\t" : ""));
			}
			System.out.println();
		}
	}
}
