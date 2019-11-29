package com.dgraph;

import java.util.ArrayList;
import java.util.Scanner;

public class JalviewUtil {
	private static class jalViewParsingState{
		private static int SCORE_LINE_STATE = 1,
					LENGTH_OF_ALIGNMENT_LINE_STATE = SCORE_LINE_STATE + 1,
					SEQUENCE_A_STATE = LENGTH_OF_ALIGNMENT_LINE_STATE + 1,
					SEQUENCE_B_STATE = SEQUENCE_A_STATE + 2; //unused
		private int state = 0;
		private void hasScore() {
			if (state != 0) throw new RuntimeException("Got score unexpectedly.");
			state = SCORE_LINE_STATE;
		}
		private void hasSequence() {
			if (state == SCORE_LINE_STATE) {
				state = SEQUENCE_A_STATE;
			} else if (state == SEQUENCE_A_STATE) {
				state = 0;
			} else {
				throw new RuntimeException("Got 'sequence' line unexpectedly.");
			}
		}
	};

	//Modifies / gets the index (in order of occurrence in the file) of the given sequence name
	public static void jalviewPairwiseAlignmentScores(String[] lines) {
		jalViewParsingState state = new jalViewParsingState();
		double score = 0;
		int alignmentLength = 0;
		String sequenceA = null, sequenceB = null;
		ArrayList<Object[]> scores = new ArrayList<Object[]>();
		for(String line : lines) {
			Scanner in = new Scanner(line);
			if (!in.hasNext()) {
				continue;
			}
			String firstToken = in.next();
			if (firstToken.equals("Score")) {
				if (in.next().equals("=")) {
					score = in.nextDouble();
					state.hasScore();
				}
			} else if (firstToken.equals("Length")) {
				if (in.next().equals("of") && in.next().equals("alignment") && in.next().equals("=")) {
					alignmentLength = in.nextInt();
				}
			} else if (firstToken.equals("Sequence")) {
				String restOfLine = in.nextLine();
				//Replace the "/d-dddd" business Jalview appends
				//System.err.println(restOfLine);
				String sequenceName = restOfLine.substring(0, restOfLine.indexOf('/')).trim();
				//System.err.println(sequenceName);
				if (state.state == jalViewParsingState.SCORE_LINE_STATE) {
					sequenceA = sequenceName;
				}
				if (state.state == jalViewParsingState.SEQUENCE_A_STATE) {
					sequenceB = sequenceName;
					//The score is a sum over all positions in the alignment of sum similarity matrix (i.e. BLOSUM)
					//Divide by length of alignment to compute a length-independent metric
					double distance = Math.max(score,1)/alignmentLength;
					//Need to make large scores
					//Square to emphasize that a 2x similarity score should lead to a 1/4 reduction in distance,
					//not a 1/2 reduction
					distance = 1/(distance*distance);
					System.out.println(sequenceA+"\t"+sequenceB+"\t"+distance);
				}
				state.hasSequence();
			}
			//System.err.println(state.state);
		}
		/*
		int Ncomp = validLines.size();
		int s = (int)((1+PApplet.sqrt(1+8*Ncomp))/2); 
		float[][] toFill = new float[s][s];
		for(int x = 0; x < s; x++){
			toFill[x][x]=0;
		}
		for(int k = 0; k < Ncomp; k++){
			vals = (int[])validLines.get(k);
			float value = PApplet.min(20,1/PApplet.pow(vals[2]/40f,2));
			toFill[vals[0]-1][vals[1]-1] = value;
			toFill[vals[1]-1][vals[0]-1] = value;
		}
		for(int y = 0; y < s; y++){
			for(int x = 0; x < s; x++){
				System.out.print(toFill[y][x]+" ");
			}
			System.out.println();
		}
		*/
	}

}
