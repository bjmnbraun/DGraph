package com.dgraph;

import java.io.PrintWriter;

import processing.core.PApplet;

public class PointCloudUtil {
    //An algorithm that computes a flip + rotation transforming an input set of 2D points to "match" another set
	//Note that this is destructive of toTransform.
	//Minimizes RMSD between the rotated set (out) and target.
	public static void transformToMatch2D(float[][] target, float[][] toTransform, float[][] out, PrintWriter dbg){
		//Uses out as temporary space until we dedicate a final transformation
		float bestMag = Float.MAX_VALUE;
		int needsAFlip = 0;
		float andTheta = 0;
		for(int doIt = 0; doIt < 2; doIt++){
			float theta = calcRotation2D(target,toTransform);
			for(float offTheta = 0; offTheta < 1; offTheta++){
				rotate2D(toTransform,out,theta);
				float mag = totalDistance(target,out);
				if (mag < bestMag){
					bestMag = mag;
					andTheta = theta;
					needsAFlip = doIt;
				}
				theta+=PApplet.PI;
			}
			//For second time, try flipping the x-axis of toTransform.
			//Note we do this both times so toTransform is back to where it started 
			for(int k = 0; k < toTransform.length; k++){
				toTransform[k][0]=-toTransform[k][0];
			}
		}
		//if the flip was better, flip to get there..
		if (needsAFlip==1){
			for(int k = 0; k < toTransform.length; k++){
				toTransform[k][0]=-toTransform[k][0];
			}
		}
		rotate2D(toTransform,out,andTheta);
		if (dbg != null) {
			dbg.printf(" %8.3f %-5d %8.3f",andTheta,needsAFlip,bestMag);
			dbg.println();
		}
	}
	//Rotates each a[i] by theta to compute out[i]
	public static void rotate2D(float[][] a, float[][] out, float theta){
		float ct = PApplet.cos(theta); 
		float st = PApplet.sin(theta);
		for(int k = 0; k < a.length; k++){
			out[k][0] = ct*a[k][0]+st*a[k][1];
			out[k][1] = -st*a[k][0]+ct*a[k][1];
		}
	}
	//Calculates a rotation to transform a set of 2D points into another.
	//
	//Assumes that both a and b have already been translated so their centroid is 0.
	//
	//This is the Kabsch algorithm, which computes the rotation that minimizes RMSD
	//https://en.wikipedia.org/wiki/Kabsch_algorithm
	public static float calcRotation2D(float[][] a, float[][] b){
		float[] H = new float[4];
		//Components of H are the cross-covariance matrix of a and b.
		for(int k = 0; k < a.length; k++){
			H[0]+=a[k][0]*b[k][0];
			H[1]+=a[k][0]*b[k][1];
			H[2]+=a[k][1]*b[k][0];
			H[3]+=a[k][1]*b[k][1];
		}
		//See
		//Witzgall, Christoph J., Javier Bernal, and James F. Lawrence. A Purely Algebraic Justification of the Kabsch-Umeyama Algorithm. No. Journal of Research (NIST JRES)-. 2019."
		//where it's shown that the justification behind Kabsch is that the rotation U that maximizes
		// tr(U H)   <-- (1)
		//in fact minimizes RMSD. The U maximizing (1) can be determined from the SVD of H, see the paper.
		//
		//Note that since U is a rotation matrix, it has the form:
		// [[cos(theta), -sin(theta)],
		//  [sin(theta),  cos(theta)]]
		//and so tr(U H) is
		//  (H[0] + H[3]) * cos(theta) + (H[1] - H[2]) * sin(theta)   <-- (2)
		//which is maximized when (cos(theta), sin(theta)) faces in the same direction as 
		//  (H[0] + H[3], H[1] - H[2]), i.e.:
		float theta = PApplet.atan2((H[1]-H[2]),(H[0]+H[3]));
		return theta;
	}
	//Computes Sum(distance(a[i], b[i]))
	public static float totalDistance(float[][] a, float[][] b){
		float toRet = 0;
		float dx, dy;
		for(int k = 0; k < a.length; k++){
			dx = a[k][0]-b[k][0];
			dy = a[k][1]-b[k][1];
			toRet += PApplet.sqrt(dx*dx+dy*dy);
		}
		return toRet;
	}
}
