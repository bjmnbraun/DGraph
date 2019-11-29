package com.pdscore;

import java.util.Arrays;

//Code for computing PD (peptide distance) scores
public class PDScore {
	private static final char[] list1 = new char[]{
			'A',	'R',	'N',	'D',	'C',	'Q',	'E',	'G',	'H',	
			'I',	'L',	'K',	'M',	'F',	'P',	'S',	'T',	'W',	'Y',	'V',	'-'};
	private static final int[] list;
	static {
		list = new int[128];
		for(int k = 0; k < list1.length; k++){
			list[list1[k]]=k;
		}
	}
	private static final double[][] atts = new double[][]{
		{
			0.354,	7.573,	11.294,	13.42,	-5.846,	6.599,	9.788,	9.655,	1.019,	-15.634,	-11.825,	10.762,	-10.585, -14.571, 7.662, 8.813,	3.012,	-13.11,	-6.245,	-12.135, 0,  }
		,
		{
			3.762,	-10.135,	1.067,	-1.6,	4.885,	-5.166,	-7.861,	15.778,	-4.969,	1.993,	0.505,	-9.517,	-3.959,	-0.646,	8.029,	6.682,	4.127,	-5.222,	-1.6,	3.818,	0,  }
		,
		{
			-11.036,	2.486,	2.718,	-0.325,	1.626,	-0.697,	-7.318,	-0.558,	0.953,	-2.045,	-6.157,	-1.022,	-3.601,	1.673,	9.456,	-0.348,	-0.348,	9.038,	9.874,	-4.345,	0,  }
		,
		{
			-0.649,	-4.291,	1.963,	3.742,	9.397,	0.582,	2.611,	0.299,	4.657,	-3.243,	-4.557,	-5.405,	5.339,	-0.033,	-3.576,	-1.131,	-2.195,	1.38,	-1.597,	-3.26,	0,  }
		,
		{
			2.828,	-5.687,	-0.859,	2.437,	-5.843,	-1.75,	4.734,	1.656,	-0.328,	-1.672,	3.219,	-0.422,	1.203,	3.25,	6,	-3.062,	-4.281,	4.64,	-1.422,	-4.672,	0,  }
	};
	public static double pdVal(char one, int att){
		int whichChar = list[(int)one];
		/*
			   for(int k = 0; k < list.length; k++){
			   if (one == list[k]){
			   whichChar = k;
			   }
			   }
		 */
		return atts[att][whichChar];
	}
	public static double PDScore(String a, String b){
		a = a.trim().toUpperCase();
		b = b.trim().toUpperCase();
		//PRE: a shorter than B.
		String smaller = a.length()<b.length()?a:b;
		String bigger = a.length()<b.length()?b:a;
		a = smaller;
		b = bigger;
		//Done.
		double toMin = Double.MAX_VALUE;
		int k = 0;
		for(k = 0; k <= bigger.length()-smaller.length(); k++){
			toMin = Math.min(toMin, PDScore0(a,b.substring(k,k+a.length())));
		}
		return toMin;
	}
	public static double PDScore0(String a, String b){
		double totSum = 0;
		double pval1 = 0;
		double pval2 = 0;
		int k = 0;
		for(k = 0; k < b.length(); k++){
			double partSum = 0;
			for(int att = 0; att < 5; att++){
				pval1 = pdVal(a.charAt(k),att);
				pval2 = pdVal(b.charAt(k),att);
				pval1 = pval1-pval2;
				partSum = partSum + pval1*pval1;
			}
			totSum += Math.sqrt(partSum);
		}
		return totSum/a.length();
	}
	//Computes the PD between two amino acids a and b. 
	public static double _PD(char a, char b){
		double partSum = 0;
		for(int att = 0; att < 5; att++){
			double pval1 = pdVal(a,att);
			double pval2 = pdVal(b,att);
			pval1 = pval1-pval2;
			partSum = partSum + pval1*pval1;
		}
		return Math.sqrt(partSum);
	}
	/**
	 * Computes the PD2 score between sequences A and B. Both may be of variable lengths.
	 **/
	public static double PDScoreWindowed(String a, String b, int wSize){
		if (a.length() > b.length()){
			String tmp = a;
			a = b;
			b = tmp;
		}
		//Shrink window size to the length of the shortest sequence.
		if (wSize > a.length()) {
			wSize = a.length();
		}
		
		//Keep track of the "best" PDScore * wSize for each window of A:
		//Note that we don't do the division by wSize until the last line of this function.
		double[] bestMatches = new double[a.length() - wSize + 1];
		Arrays.fill(bestMatches, Double.MAX_VALUE);
		
		//For all possible ways to place the first wSize window of a onto b
		for(int shift = 0; shift <= b.length() - wSize; shift++) {
			//Initialize a sliding window
			double slidingWindowPD = 0;
			int off = 0;
			//See above for why this is never greater than a's length
			for(; off < wSize - 1; off++) {
				slidingWindowPD += _PD(a.charAt(off), b.charAt(shift + off));
			}
			//For all wSize windows in a, which we are matching to b shift-many characters later.
			//The window ends on character "off".
			for(; off < a.length() && off + shift < b.length();)  {
				//Integrate a[off]
				slidingWindowPD += _PD(a.charAt(off), b.charAt(shift + off));
				//Advance off++
				off++;
				//At this point we have a full window of wSize PD values in slidingWindowPD ending before off.
				//i.e. the segment of a: [off - wSize, off) with some window in b
				double oldBest = bestMatches[off - wSize];
				if (slidingWindowPD < oldBest) {
					bestMatches[off - wSize] = slidingWindowPD;
				}
				//Advance and subtract off the first character in the window 
				slidingWindowPD -= _PD(a.charAt(off - wSize), b.charAt(shift + off - wSize));
			}
		}
		//Compute the average of bestMatches
		double sum = 0;
		for(double x : bestMatches) {
			sum += x;
		}
		//bestMatches are in fact wSize * PDScore, so handle that here:
		return sum / bestMatches.length / wSize;
	}


	public static void pdScoreMat(String[] lines){
		for(int k = 0; k < lines.length; k++){
			for(int p = 0; p < lines.length; p++){
				double val = PDScore(lines[p],lines[k]);
				System.out.printf("%.3f ",val);
			}
			System.out.println();
			System.err.println(k+" out of "+lines.length);
		}
	}


}
